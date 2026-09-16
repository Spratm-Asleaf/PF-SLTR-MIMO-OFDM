%{
    Online supplementary materials of the paper titled:
    A New Particle Filter for Target Tracking in MIMO OFDM Integrated Sensing and Communications

    @Author:    Shixiong Wang (s.wang@xjtu.edu.cn; s.wang@u.nus.edu; wsx.gugo@gmail.com)
    @Institute: School of Mathematics and Statistics, Xi'an Jiaotong University
    @Date:      17 November 2025, 27 April 2026, 4 September 2026
    @Home:      https://github.com/Spratm-Asleaf/PF-SLTR-MIMO-OFDM
%}


%% OFDM MIMO ISAC System Simulation
addpath('./Utils/');

%% Start
clear;
clc;
rng(2025);                  % For reproducibility; Year of 2025, is this a good reason to use 2025 here? Surely good!

%% Environment Setups
SNR = 0;                    % Signal-to-noise ratio in dB (NB: this experiment is not sensitive to SNR due to large MIMO)
Pn = 1;                     % Channel noise power
Pt = (10^(SNR/10))*Pn;      % Transmit power
C  = 299792458;             % Speed of light (m/s)

%% ISAC System Setups - Monostatic OFDM MIMO
% Base Station Position
p_s = [
    0
    0
];

% Operating Frequency
fc = 10e9;                  % Center carrier frequency (GHz); Doppler freq is propotional to this value; the larger, the more sensitive to Doppler
lambda = C/fc;              % Wavelength (m)

% MIMO Configuration
Nt = 64;                    % Number of transmit antennas
Nr = 64;                    % Number of receive antennas
a = @(Ang)StrVec(Nt, Ang);  % Tx steering vector
b = @(Ang)StrVec(Nr, Ang);  % Rx steering vector

% Pulse (Single) Configuration
Tc = 1e-6;                  % CP length: 1/8 length of OFDM symbol duration
To = 5*Tc;                  % OFDM symbol length in time (Without CP)
Delta_f = 1/To;             % OFDM subcarrier spacing
M  = 1;                     % Number of OFDM symbols per pulse
Tp = M * (To + Tc);         % Pulse duration  

Nc = 256;                   % Number of subcarrier; i.e., Number of samples in a OFDM symbol
B  = Nc * Delta_f;          % Bandwidth (Total spectra span: [-B/2, B/2])
fs = B;                     % Sampling rate = Bandwidth = 2 * Nyquist frequency = 2 * (B/2) = B
Ts = 1/fs;                  % Sampling interval

Ncp = ceil(Tc/Ts);          % Number of samples in a CP

Lss = M*(Nc + Ncp);         % Total samples in a pulse; M sub-pulses included in a pulse; each sub-pulse has (Nc + Ncp) samples
assert(Lss == ceil(Tp/Ts));

% Service Region Limits
% Service prioity: specify max range or specify max velocity?
ServicePrioity = 'Range';       % 'Range', 'Velocity'

switch ServicePrioity
    case 'Range'
        R_max = 2.1e3;                % Maximum-allowed range
        tau_max = GetTau(R_max, C);
        
        Tr = tau_max + Tp;          % Pulse repetition interval
        fr = 1/Tr;                  % Pulse repetition frequency (Hz)
        
        nu_max = fr/2;              % Maximum-allowed doppler frequency shift; cannot be changed because it is determined by coherent-processing mechanism
        v_max = GetVelocity(nu_max, lambda);     % Maximum-allowed velocity
        
    case 'Velocity'
        v_max = 85;
        nu_max = GetDoppler(v_max, lambda);

        fr = 2*nu_max;
        Tr = 1/fr;

        tau_max = Tr - Tp;
        R_max = GetR(tau_max, C);
        assert(R_max > 0);

    otherwise
        error('main :: Error in ServicePrioity');
end

% Time integration interval; i.e., Number of samples within a CPI Ti
L = floor(Tr/Ts);                  

% Maximum DoA: <= 60 deg: 
%   NB: Although half-wavelength spacing allows 90 deg viewfield, to ensure good DoA estimation performance, do not let array's electronical angle exceeds 60 deg
theta_max = deg2rad(60);

% Pulse Train Configuration
Np = 1;                 % Number of pulses in a pulse train; if Np = 1, no doppler processing is used
Ti = Np * Tr;           % Coherent Processing Interval (CPI)

% Show Radar Performance Indicators
disp(['Maximum Range       : ' num2str(R_max) ' m']);
disp(['Minimim Range       : ' num2str(GetR(Tp, C)) ' m']);      % Range blind zone in pulsed radar systems
disp(['Maximum Velocity    : ' num2str(v_max) ' m/s']);
disp(['Maximum Angle       : ' num2str(rad2deg(theta_max)) ' Deg']);
disp(['Delay    Resolution : ' num2str(Ts) ' s']);               % Delay resolution is 1/B = 1/fs = Ts
disp(['Range    Resolution : ' num2str(GetR(Ts, C)) ' m']);
disp(['Doppler  Resolution : ' num2str(1/Ti) ' Hz']);            % Doppler frequency resolution is 1/Ti
disp(['Velocity Resolution : ' num2str(GetVelocity(1/Ti, lambda)) ' m/s']);
disp(['Ghost Spacing       : ' num2str(GetR(To, C)) ' m']);
disp(['ISAC Period         : ' num2str(Ti) ' s']);               % An ISAC period (which functions sensing using OFDM) is a CPI "Ti"

%% Target Parameters
Tt   = 0.05;           % Real-time tracking interval is 25 ms; in practice, 50 ms is also a good choice 
assert(Tt >= Ti);      % A tracking window is no smaller than an ISAC period
disp(['Tracking Window     : ' num2str(Tt) ' s']);               % Tracking time window
TimeTrack = 0:Tt:100;  % Track the object for 100s
LenTimeTrack = length(TimeTrack);

% Cartesian Coordinate System
State0 = zeros(4, LenTimeTrack);    % True states: ground truth is marked with "0"

% Polar Coordinate System
R0 = zeros(1, LenTimeTrack);        % Range
Theta0 = zeros(1, LenTimeTrack);    % DoA
V0 = zeros(1, LenTimeTrack);        % Radial Speed

% State transition matrix
% Constant Velocity (CV) model
F_CV = [
    1   Tt
    0   1
];
F = blkdiag([F_CV, zeros(size(F_CV)); zeros(size(F_CV)), F_CV]);

% Noise driven matrix
G_CV = [
    Tt^2/2
    Tt
];
G = blkdiag([G_CV, zeros(size(G_CV)); zeros(size(G_CV)), G_CV]);

% Generate True Trajectory
TrajectoryMode = 2;                  % 1: Almost-Straightline, 2: Curved, 3: Circled

for t_track = 1:LenTimeTrack
    if t_track == 1
        % Initial true state
        State0(:, t_track) = [
            200                      % Position in x axis
            5                        % Velocity in x axis
            -250                     % Position in y axis
            8                        % Velocity in y axis
        ];

    else
        matQ_CV = eye(2);

        switch TrajectoryMode
            case 1
                State0(:, t_track) = F * State0(:, t_track - 1) + G * chol(matQ_CV, 'lower') * randn(2, 1);

            case 2
                Delta_Vx = 0.02*sin(0.1*TimeTrack(t_track));
                Delta_Vy = 0.02*cos(0.1*TimeTrack(t_track));
                State0(:, t_track) = F * State0(:, t_track - 1) + diag([0, 1, 0, 1]) * [0; Delta_Vx; 0; Delta_Vy] + G * chol(matQ_CV, 'lower') * randn(2, 1);

            case 3
                center = [5000; 5000];
                cTheta = linspace(0, pi, LenTimeTrack) + atan2(center(2) - State0(3, 1), center(1) - State0(1, 1));
                Delta_Vx = 0.005*cos(cTheta(t_track))*4;
                Delta_Vy = 0.005*sin(cTheta(t_track))*4;
                State0(:, t_track) = F * State0(:, t_track - 1) + diag([0, 1, 0, 1]) * [0; Delta_Vx; 0; Delta_Vy] + G * chol(matQ_CV, 'lower') * randn(2, 1);

            otherwise 
                error('main :: Error in TrajectoryMode');
        end
    end

    % True States in Cartesian Coordinate System
    p_obj = State0([1;3], t_track);         % Cartesian two-dimensional positions
    v_obj_cart = State0([2;4], t_track);    % Cartesian two-dimensional velocities (NB: not radial)

    % True States in Polar/Radar Coordinate System
    R0(t_track)      = norm(p_s - p_obj, 2);
    Theta0(t_track)  = atan2(p_obj(2) - p_s(2), p_obj(1) - p_s(1));     % positive for counterclockwise
    V0(t_track)      = v_obj_cart(1) * cos(Theta0(t_track)) + v_obj_cart(2) * sin(Theta0(t_track));       % positive for receding
end

% Validate Data: Warning: Do not remove the below "assert" statements unless you know what you are doing
assert(sum(R0 <= R_max) == length(R0));
Tau0 = GetTau(R0, C);       % Round-trip delay

assert(sum(abs(V0) <= v_max) == length(V0));
Nu0  = GetDoppler(V0, lambda);     % Doppler freq shift

assert(sum(abs(Theta0) <= theta_max) == length(Theta0));

%% Scatters Parameters
Scatters = [
    300  450  1000  800  1300  % Range
    20   -30  35    -8   15    % Theta
];
[~, Ns] = size(Scatters);      % Ns: Total number of significant scatters and other targets
% Ns = 0;

%% Plot States: Trajectory, Velocities, DoA
if true
    figure;
    plot(p_s(1), p_s(2), 'rs', 'linewidth', 1.5, 'markersize', 8)
    hold on;
    plot(State0(1, :), State0(3, :), 'b', 'linewidth', 2);
    % axis([-50 700 -700 50])
    set(gca, 'XAxisLocation', 'top')
    % axis equal;
    for j = 1:Ns
        plot(Scatters(1, j) * cos(deg2rad(Scatters(2, j))), Scatters(1, j) * sin(deg2rad(Scatters(2, j))), 'mo', 'linewidth', 1.5, 'markersize', 8);
    end
    xlabel('Cartesian Position x-Axis');
    ylabel('Cartesian Position y-Axis');
    
    figure;
    plot(TimeTrack, R0)
    xlabel('Time');
    ylabel('Range');
    
    figure;
    plot(TimeTrack, State0(2, :), 'r', TimeTrack, State0(4, :), 'b--')
    xlabel('Time');
    ylabel('Cartesian Velocities');
    legend('Vx', 'Vy');
    
    figure;
    plot(TimeTrack, V0)
    xlabel('Time');
    ylabel('Radial Velocity');
    
    figure;
    plot(TimeTrack, rad2deg(Theta0))
    xlabel('Time');
    ylabel('Degree of Arrival (DoA)');
    
end

% return;

%% Tracking
MonteCarlo = 1;

MSE_KF_ILTR   = zeros(MonteCarlo, LenTimeTrack);    
RunningTime_KF_ILTR    = zeros(MonteCarlo, LenTimeTrack);

MSE_PF_SLTR = zeros(MonteCarlo, LenTimeTrack);    
RunningTime_PF_SLTR  = zeros(MonteCarlo, LenTimeTrack);

for mc = 1:MonteCarlo
    % Estimated State
    State_hat_qAvrg_PF_SLTR = zeros(4, LenTimeTrack);   % Trimed (q-quantile) Posterior Average

    State_hat_KF_ILTR = zeros(4, LenTimeTrack);
    
    % Particles at each time step
    N_par = 200;                               % Number of Particles
    mu_prior = zeros(1, N_par) + 1/N_par;      % Prior Weights
    mu_likelihood = zeros(1, N_par) + 1/N_par; % Likelihood Weights
    State_Particles_SLTR = zeros(4, N_par);         % State Particles
    
    % Start
    for t_track = 1:LenTimeTrack
        if mod(t_track, 50) == 0
            disp(['Progress: ' num2str(100*t_track/LenTimeTrack) '%']);
        end
    
        %% Current-Time True Values
        tau0 = Tau0(t_track);
        nu0 = Nu0(t_track);
        theta0 = Theta0(t_track);
    
        %% Generate Measurements in a CPI
        % Tx: OFDM Waveform Generation: Time Domain
        S = cell(Np, 1);
        modOrder = 64;              % QAM Modulation order (e.g., 4 for QPSK; 16 for 16-QAM)
        for p = 1:Np
            S_p = zeros(Nt, Lss);
            for tx = 1:Nt
                data = randi([0 modOrder-1], Nc, M);
                block = qammod(data, modOrder, 'UnitAveragePower', true);       % OFDM resource block: Nc subcarriers and M symbols
                timeData = ifft(block, Nc, 1);                                  % IFFT to get time domain signal
                CPAdded  = [timeData((end-Ncp+1):end, :); timeData];            % Add CP
                S_p(tx, :) = (CPAdded(:))';                                     % Make it sequential 
            end
    
            % Zero Padding
            S{p} = [sqrt(Pt) * S_p, zeros(Nt, L - Lss)];                                    % Zero padding to simulate linear convolution
        end
        
        %% Signal Propogation: Time Domain
        % Scatters' Channel Effects: Scatters (including other targets)
        ScatterChannelEffects = cell(Np, 1);
        
        % Attenuation factor
        ScatterAttenuateFactor = (0.9 + 0.1*(2*rand(1, Ns) - 1)) .* exp(1j*2*pi*randn(1, Ns));
        for p = 1:Np
            ScatterChannelEffects{p} = zeros(Nr, L);
    
            for j = 1:Ns
                scatter_tau = GetTau(Scatters(1, j), C);
                scatter_theta = deg2rad(Scatters(2, j));
    
                % Scatters' Effects: Warning: Do not remove these parentheses unless you know what you are doing
                ScatterChannelEffects{p} = ScatterChannelEffects{p} + ... 
                                            ScatterAttenuateFactor(j) * ...
                                                ((b(scatter_theta) * a(scatter_theta)') * (GetTimeShiftedS(S{p}, floor(scatter_tau*fs + 0.5))));
            end
        end
    
        % Rx: Each receive pulse is Nt-by-(Nc*M)
        Y = cell(Np, 1);
    
        % Channel Attenuation Coefficient
        AttenuateFactor  = GeneratePathGain(R0(1), R0(t_track), 1);  % beta_{0, k}

        % Spatial Effect
        SpatialEffect   = b(theta0) * a(theta0)';
        
        for p = 1:Np
            % Doppler Effect
            DopplerEffect   = exp(1j*2*pi*nu0*(p-1)*Tr);
            % Time Delay Effect
            TimeShiftedTx = GetTimeShiftedS(S{p}, floor(tau0*fs));
            % Channel Noise
            NoiseEffect     = chol(Pn, 'lower') * (1/sqrt(2)) * (randn(Nr, L) + 1j * randn(Nr, L));
            
            % Overall Transmission: For each receive antenna n in [1:Nr], we have Nc*M = L samples
            Y{p} =     ... % Target Echo: Warning: Do not remove these parentheses unless you know what you are doing
                        AttenuateFactor * (SpatialEffect * (DopplerEffect * TimeShiftedTx))  ...   
                        ...
                        ... % Scatter Echoes
                        + ScatterChannelEffects{p} ...                                                 
                        ... 
                        ... % Channel Noise
                        + NoiseEffect;                                                                  
        end
    
        %% Sensing Rx Processing
        tic;

        % State Transition
        isWindowing = 0;
    
        if t_track == 1
            matR = diag([0.1, 0.01, deg2rad(0.1), deg2rad(0.01)]);
            v_cart = State0([2;4], t_track);
            v_polar = GetPolarVelocity(v_cart, theta0);
            State_Particles_SLTR = [GetR(tau0, C); v_polar(1); theta0; v_polar(2)/GetR(tau0, C)] + chol(matR, 'lower') * randn(4, N_par);

            State_hat_KF_ILTR(:, t_track) = [GetR(tau0, C); v_polar(1); theta0; v_polar(2)/GetR(tau0, C)] + chol(matR, 'lower') * randn(4, 1);
            P_KF_ILTR = eye(4, 4);

        else
            matQ = diag([20, deg2rad(10)]);
            State_Particles_SLTR = F * State_Particles_SLTR + G * chol(matQ, 'lower') * randn(2, N_par);

            State_hat_KF_ILTR(:, t_track) = F * State_hat_KF_ILTR(:, t_track-1);
            P_KF_ILTR = F*P_KF_ILTR*F' + G*diag([20, deg2rad(10)]/20)*G';
        end

        % Calculate Likelihoods
        Utility = zeros(1, N_par);
        for i = 1:N_par
            R      = State_Particles_SLTR(1, i);
            v      = State_Particles_SLTR(2, i);
            theta  = State_Particles_SLTR(3, i);
    
            Utility(i) = EvaluateDelayDopplerDoA_Pulsed(GetTau(R, C), GetDoppler(v, lambda), theta, Ts, Y, S, Lss, Np, Nr, L, b, a, Tr, isWindowing);
        end
    
        % Likelihood Weights
        Utility = (abs(Utility)).^2;
        mu_likelihood_sum = sum(Utility);
        mu_likelihood = Utility/mu_likelihood_sum;
    
        % Apply Generalized Bayes Rule
        xi = 1;
        mu_posterior = mu_prior .* (mu_likelihood.^xi);
        mu_posterior_sum = sum(mu_posterior);
        mu_posterior = mu_posterior/mu_posterior_sum;
    
        % Systematic Resampling
        if 1/(sum(mu_posterior.^2)) < N_par*0.5   % or 0.25, 0.3 are also good practice 
            [State_Particles_SLTR, mu_posterior] = Resampling(State_Particles_SLTR, mu_posterior);
        end
    
        mu_quantile = quantile(mu_posterior, 0.6);
        ind = find(mu_posterior >= mu_quantile);
        mu_posterior_part = mu_posterior(ind);
        mu_posterior_part = mu_posterior(ind)/sum(mu_posterior_part);
        State_hat_qAvrg_PF_SLTR(:, t_track) = State_Particles_SLTR(:, ind) * mu_posterior_part';
    
        % Prior at this time is the posterior for the next time
        mu_prior = mu_posterior;

        RunningTime_PF_SLTR(mc, t_track) = toc;

        Pos = [
            State_hat_qAvrg_PF_SLTR(1, t_track) .* cos(State_hat_qAvrg_PF_SLTR(3, t_track))
            State_hat_qAvrg_PF_SLTR(1, t_track) .* sin(State_hat_qAvrg_PF_SLTR(3, t_track))
        ];

        MSE_PF_SLTR(mc, t_track)  = norm(State0([1;3], t_track) - Pos, 'fro')^2;

        %% KF-ILTR
        H_KF_ILTR = [
            1 0 0 0
            0 1 0 0
            0 0 1 0
        ];

        tic;

        % Matched Filter for range-Dopper-DoA estimation
        if true
            GridLen = 30;
            for ml = 1:GridLen^3
                % To simulate how long it takes to evaluate the ambiguity-function value at all grid points in the "rang-Doppler-DoA" bins;
                %   For each dimesion, only GridLen = 30 points are evaluated for time-saving purpose
                EvaluateDelayDopplerDoA_Pulsed(GetTau(500, C), 0, 0, Ts, Y, S, Lss, Np, Nr, L, b, a, Tr, isWindowing);
            end
        end
        % Suppose we have found the maximum-likelihood estiamtes; We even suppose the measurements are unbiased, and with small measurement errors:
        % ---This assumption is very beneficial for KF-ILTR
        y = [GetR(tau0, C); GetVelocity(nu0, lambda); theta0] + chol(diag([1, 0.1, 0.01]), 'lower') * randn(3, 1);
        
        % Kalman updates
        K_KF_ILTR = P_KF_ILTR*H_KF_ILTR'*(H_KF_ILTR*P_KF_ILTR*H_KF_ILTR' + diag([1, 0.1, 0.01]))^-1;
        State_hat_KF_ILTR(:, t_track) = State_hat_KF_ILTR(:, t_track) + K_KF_ILTR*(y - H_KF_ILTR*State_hat_KF_ILTR(:, t_track));
        P_KF_ILTR = (eye(4) - K_KF_ILTR * H_KF_ILTR) * P_KF_ILTR * (eye(4) - K_KF_ILTR * H_KF_ILTR)' + K_KF_ILTR * diag([1, 0.1, 0.01]) * K_KF_ILTR';

        RunningTime_KF_ILTR(mc, t_track) = toc;

        Pos_ILTR = [
            State_hat_KF_ILTR(1, t_track) .* cos(State_hat_KF_ILTR(3, t_track))
            State_hat_KF_ILTR(1, t_track) .* sin(State_hat_KF_ILTR(3, t_track))
        ];
        MSE_KF_ILTR(mc, t_track)  = norm(State0([1;3], t_track) - Pos_ILTR, 'fro')^2/2;
    end
end

%% Plot
figure;
plot(State0(1, :), State0(3, :), 'r', 'linewidth', 2);
hold on;
Pos_SLTR = [
    State_hat_qAvrg_PF_SLTR(1, :) .* cos(State_hat_qAvrg_PF_SLTR(3, :))
    State_hat_qAvrg_PF_SLTR(1, :) .* sin(State_hat_qAvrg_PF_SLTR(3, :))
];
plot(Pos_SLTR(1, :), Pos_SLTR(2, :), 'b--', 'linewidth', 2);                               
mean(MSE_PF_SLTR, 2)
mean(RunningTime_PF_SLTR, 2)
% axis equal;
% axis([-50 700 -700 50])
set(gca, 'XAxisLocation', 'top')

figure;
plot(State0(1, :), State0(3, :), 'r', 'linewidth', 2);
hold on;
Pos_ILTR = [
    State_hat_KF_ILTR(1, :) .* cos(State_hat_KF_ILTR(3, :))
    State_hat_KF_ILTR(1, :) .* sin(State_hat_KF_ILTR(3, :))
];
plot(Pos_ILTR(1, :), Pos_ILTR(2, :), 'b--', 'linewidth', 2);                               
mean(MSE_KF_ILTR, 2)
mean(RunningTime_KF_ILTR, 2)
% axis equal;
% axis([-50 700 -700 50])
set(gca, 'XAxisLocation', 'top')