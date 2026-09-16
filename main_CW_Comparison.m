%{
    Online supplementary materials of the paper titled:
    A New Particle Filter for Target Tracking in MIMO OFDM Integrated Sensing and Communications

    @Author:    Shixiong Wang (s.wang@xjtu.edu.cn; s.wang@u.nus.edu; wsx.gugo@gmail.com)
    @Institute: School of Mathematics and Statistics, Xi'an Jiaotong University
    @Date:      17 November 2025, 27 April 2026, 4 September 2026
    @Home:      https://github.com/Spratm-Asleaf/PF-SLTR-MIMO-OFDM
%}


%% OFDM MIMO ISAC System Simulation - Continuous Wave Systems
% OFDM Communications Systems Reuse: full-duplex, no silent time between pulses are mandatory, although optional; here, one pulse means one OFDM burst
% All computational operations are in frequency domain, due to the perservation of circular convolution using CP
% Condition: Max Target Delay <= CP Length

addpath('./Utils/');

%% Start
clear;
clc;
rng(20250612);                  % For reproducibility; Year of 2025, is this a good reason to use 2025 here? Surely good!

%% Environment Setups
SNR = -10;                  % Signal-to-noise ratio in dB (NB: this experiment is not sensitive to SNR due to large MIMO)
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
Tc = 1e-6;
To = 5*Tc;                  % OFDM symbol length in time (Without CP)
Delta_f = 1/To;             % OFDM subcarrier spacing
M  = 1;                     % Number of OFDM symbols per pulse: Warning: Do Not Change This Value! (One OFDM symbol per pulse!)
Tp = M * (To + Tc);         % Pulse duration  

Nc = 256;                   % Number of subcarrier
B  = Nc * Delta_f;          % Bandwidth (Total spectra span: [-B/2, B/2])
fs = B;                     % Sampling rate = Bandwidth = 2 * Nyquist frequency = 2 * (B/2) = B
Ts = 1/fs;                  % Sampling interval

% Service Region Limits
tau_max = Tc;               % Maximum-allowed round-trip delay; CP duration
R_max = GetR(tau_max, C);

Tr = Tp;                    % Pulse repetition interval
fr = 1/Tr;                  % Pulse repetition frequency (Hz)

nu_max = fr/2;              % Maximum-allowed doppler frequency shift; cannot be changed because it is determined by coherent-processing mechanism
v_max = GetVelocity(nu_max, lambda);     % Maximum-allowed velocity

% Maximum DoA: <= 60 deg: 
%   NB: Although half-wavelength spacing allows 90 deg viewfield, to ensure good DoA estimation performance, do not let array's electronical angle exceeds 60 deg
theta_max = deg2rad(60);

% Pulse Train Configuration
Np = 1;                 % Number of pulses in a pulse train; if Np = 1, no doppler processing is used
Ti = Np * Tr;           % Coherent Processing Interval (CPI)

% Show Radar Performance Indicators
disp(['Maximum Range       : ' num2str(R_max) ' m']);
disp(['Minimim Range       : ' num2str(0) ' m']);
disp(['Maximum Velocity    : ' num2str(v_max) ' m/s']);
disp(['Maximum Angle       : ' num2str(rad2deg(theta_max)) ' Deg']);
disp(['Delay    Resolution : ' num2str(Ts) ' s']);               % Delay resolution is 1/B = 1/fs = Ts
disp(['Range    Resolution : ' num2str(GetR(Ts, C)) ' m']);
disp(['Doppler  Resolution : ' num2str(1/Ti) ' Hz']);            % Doppler frequency resolution is 1/Ti
disp(['Velocity Resolution : ' num2str(GetVelocity(1/Ti, lambda)) ' m/s']);
disp(['Ghost Spacing       : ' num2str(GetR(To, C)) ' m']);
disp(['ISAC Period         : ' num2str(Ti) ' s']);               % An ISAC period (which functions sensing using OFDM) is a CPI "Ti"

%% Target Parameters
Tt   = 0.025*2;          % Real-time tracking interval is 25 ms; in practice, 50 ms is also a good choice 
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
matQ_CV = eye(2);

for t_track = 1:LenTimeTrack
    assert(Tc <= 1.01e-6 && Tc >= 0.99e-6);          % Set Tc = 1e-6;
    if t_track == 1
        % Initial true state
        State0(:, t_track) = [
            50                      % Position in x axis
            2                        % Velocity in x axis
            -50                     % Position in y axis
            3                        % Velocity in y axis
        ];
    else
        center = [-100; -100];
        cTheta = linspace(0, pi, LenTimeTrack) + atan2(center(2) - State0(3, 1), center(1) - State0(1, 1));
        Delta_Vx = 0.005*cos(cTheta(t_track))/2;
        Delta_Vy = 0.005*sin(cTheta(t_track))/2;
        State0(:, t_track) = F * State0(:, t_track - 1) + diag([0, 1, 0, 1]) * [0; Delta_Vx; 0; Delta_Vy] + G * chol(matQ_CV, 'lower') * randn(2, 1);
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
    30   45   100   80   130   % Range
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
    % set(gca, 'XAxisLocation', 'top')
    % axis equal;
    for j = 1:Ns
        plot(Scatters(1, j) * cos(deg2rad(Scatters(2, j))), Scatters(1, j) * sin(deg2rad(Scatters(2, j))), 'mo', 'linewidth', 1.5, 'markersize', 8);
    end
    xlabel('Cartesian Position x-Axis');
    ylabel('Cartesian Position y-Axis');
    set(gca, 'fontsize', 16)
    
    figure;
    plot(TimeTrack, R0, 'b', 'linewidth', 2)
    xlabel('Time');
    ylabel('Range');
    set(gca, 'fontsize', 16)
    
    figure;
    plot(TimeTrack, State0(2, :), 'r', TimeTrack, State0(4, :), 'b--')
    xlabel('Time');
    ylabel('Cartesian Velocities');
    legend('V_x', 'V_y');
    set(gca, 'fontsize', 16)
    
    figure;
    plot(TimeTrack, V0, 'b', 'linewidth', 2)
    xlabel('Time');
    ylabel('Radial Velocity');
    set(gca, 'fontsize', 16)
    
    figure;
    plot(TimeTrack, rad2deg(Theta0), 'b', 'linewidth', 2)
    xlabel('Time');
    ylabel('DoA (Degree)');
    set(gca, 'fontsize', 16)
    
end

% return;

%% Tracking
MonteCarlo = 1;     % Use 1 is statistically sufficient because we have time-averaged performance over tracking time steps "k"

MSE_PF_SLTR   = zeros(MonteCarlo, LenTimeTrack);    
RunningTime_PF_SLTR    = zeros(MonteCarlo, LenTimeTrack);

MSE_KF_ILTR   = zeros(MonteCarlo, LenTimeTrack);    
RunningTime_KF_ILTR    = zeros(MonteCarlo, LenTimeTrack);

MSE_PF_SLTR_A   = zeros(MonteCarlo, LenTimeTrack);    
RunningTime_PF_SLTR_A    = zeros(MonteCarlo, LenTimeTrack);

for mc = 1:MonteCarlo
    disp(['mc = ' num2str(mc)]);

    % Estimated State
    State_hat_qAvrg_PF_SLTR = zeros(4, LenTimeTrack);   % Trimed (q-quantile) Posterior Average
    State_hat_KF_ILTR = zeros(4, LenTimeTrack);
    State_hat_qAvrg_PF_SLTR_A = zeros(6, LenTimeTrack);

    % Particles at each time step
    N_par = 200;                               % Number of Particles

    mu_prior_SLTR = zeros(1, N_par) + 1/N_par;      % Prior Weights
    mu_likelihood_SLTR = zeros(1, N_par) + 1/N_par; % Likelihood Weights
    State_Particles_SLTR = zeros(4, N_par);         % State Particles

    mu_prior_SLTR_A = zeros(1, N_par) + 1/N_par;
    mu_likelihood_SLTR_A = zeros(1, N_par) + 1/N_par; 
    State_Particles_SLTR_A = zeros(6, N_par); 
    
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
        % Tx: OFDM Waveform Generation; Frequency Domain
        Sf = cell(Np, 1);
        for p = 1:Np
            modOrder = 64;              % QAM Modulation order (e.g., 4 for QPSK; 16 for 16-QAM)
            data = randi([0 modOrder-1], Nt, Nc*M);
        
            % Each transmit pulse is Nt-by-(Nc*M)
            Sf{p} = sqrt(Pt) * qammod(data, modOrder, 'UnitAveragePower', true);
        end
        
        %% Signal Propogation; Frequency Domain
        % Scatters' Channel Effects: Scatters (including other targets)
        ScatterChannelEffects = cell(Np, 1);
    
        % Attenuation factor
        ScatterAttenuateFactor = (0.9 + 0.1*(2*rand(1, Ns) - 1)) .* exp(1j*2*pi*randn(1, Ns));
        for p = 1:Np
            ScatterChannelEffects{p} = zeros(Nr, Nc*M);
    
            for j = 1:Ns
                scatter_tau = GetTau(Scatters(1, j), C);
                scatter_theta = deg2rad(Scatters(2, j));
                
                % Scatters' Effects: Warning: Do not remove these parentheses unless you know what you are doing
                ScatterChannelEffects{p} = ScatterChannelEffects{p} + ... 
                                            ScatterAttenuateFactor(j) * ...
                                                ((b(scatter_theta) * a(scatter_theta)') * (repmat(exp(-1j*2*pi*scatter_tau*(0:Nc-1)*Delta_f), Nt, M) .* Sf{p}));
            end
        end
    
        % Rx: Each receive pulse is Nt-by-(Nc*M)
        Yf = cell(Np, 1);
    
        % Channel Attenuation Coefficient
        AttenuateFactor  = GeneratePathGain(R0(1), R0(t_track), 1);  % beta_{0, k}

        % Spatial Effect
        SpatialEffect   = b(theta0) * a(theta0)';
        
        for p = 1:Np
            % Doppler Effect
            DopplerEffect   = exp(1j*2*pi*nu0*(p-1)*Tr);
            % Time Delay Effect
            TimeDelayEffect = repmat(exp(-1j*2*pi*tau0*(0:Nc-1)*Delta_f), Nt, M);     % Every Frequency Component is Phase Shifted
            % Channel Noise
            NoiseEffect     = chol(Pn, 'lower') * (1/sqrt(2)) * (randn(Nr, Nc*M) + 1j*randn(Nr, Nc*M));
            
            % Overall Transmission: For each receive antenna n in [1:Nr], we have Nc*M = L samples
            Yf{p} =     ... % Target Echo: Warning: Do not remove these parentheses unless you know what you are doing
                        AttenuateFactor * (SpatialEffect * (DopplerEffect * (TimeDelayEffect .*  Sf{p})))  ...   
                        ...
                        ... % Scatter Echoes
                        + ScatterChannelEffects{p} ...                                                 
                        ... 
                        ... % Channel Noise
                        + NoiseEffect;                                                                  
        end
    
        %% Sensing Rx Processing

        % State Transition
        isWindowing = 0;
    
        tic;

        if t_track == 1
            matR = diag([0.1, 0.01, deg2rad(0.1), deg2rad(0.01)]);
            v_cart = State0([2;4], t_track);
            v_polar = GetPolarVelocity(v_cart, theta0);

            State_Particles_SLTR = [GetR(tau0, C); v_polar(1); theta0; v_polar(2)/GetR(tau0, C)] + chol(matR, 'lower') * randn(4, N_par);

            State_hat_KF_ILTR(:, t_track) = [GetR(tau0, C); v_polar(1); theta0; v_polar(2)/GetR(tau0, C)] + chol(matR, 'lower') * randn(4, 1);
            P_KF_ILTR = eye(4, 4);

            State_Particles_SLTR_A = [GetR(tau0, C); v_polar(1); theta0; v_polar(2)/GetR(tau0, C); real(AttenuateFactor); imag(AttenuateFactor)] + chol(blkdiag(matR, diag([0.01, 0.01])), 'lower') * randn(6, N_par);
        else
            State_Particles_SLTR = F * State_Particles_SLTR + G * chol(diag([20, deg2rad(10)]), 'lower') * randn(2, N_par);

            State_hat_KF_ILTR(:, t_track) = F * State_hat_KF_ILTR(:, t_track-1);
            P_KF_ILTR = F*P_KF_ILTR*F' + G*diag([20, deg2rad(10)]/20)*G';

            State_Particles_SLTR_A = blkdiag(F, eye(2)) * State_Particles_SLTR_A + blkdiag(G, Tt*eye(2)) * chol(diag([20, deg2rad(10), 0.05, 0.05]), 'lower') * randn(4, N_par);
        end

        ParticleTransitionTime = toc;

        %% PF-SLTR
        tic;

        % Calculate Likelihoods
        Utility_SLTR = zeros(1, N_par);
        for i = 1:N_par
            R      = State_Particles_SLTR(1, i);
            v      = State_Particles_SLTR(2, i);
            theta  = State_Particles_SLTR(3, i);
    
            Utility_SLTR(i) = EvaluateDelayDopplerDoA_CW(GetTau(R, C), GetDoppler(v, lambda), theta, Yf, Sf, Delta_f, b, a, Tr, Nc, M, isWindowing);
        end

        % Likelihood Weights
        Utility_SLTR = (abs(Utility_SLTR)).^2;
        mu_likelihood_sum_SLTR = sum(Utility_SLTR);
        mu_likelihood_SLTR = Utility_SLTR/mu_likelihood_sum_SLTR;
    
        % Apply Generalized Bayes Rule
        mu_posterior_SLTR = mu_prior_SLTR .* mu_likelihood_SLTR;
        mu_posterior_sum_SLTR = sum(mu_posterior_SLTR);
        mu_posterior_SLTR = mu_posterior_SLTR/mu_posterior_sum_SLTR;
    
        % Systematic Resampling
        if 1/(sum(mu_posterior_SLTR.^2)) < N_par*0.5   % or 0.25, 0.3 are also good practice 
            [State_Particles_SLTR, mu_posterior_SLTR] = Resampling(State_Particles_SLTR, mu_posterior_SLTR);
        end

        mu_quantile_SLTR = quantile(mu_posterior_SLTR, 0.6);
        ind_SLTR = find(mu_posterior_SLTR >= mu_quantile_SLTR);
        mu_posterior_part_SLTR = mu_posterior_SLTR(ind_SLTR);
        mu_posterior_part_SLTR = mu_posterior_SLTR(ind_SLTR)/sum(mu_posterior_part_SLTR);
        State_hat_qAvrg_PF_SLTR(:, t_track) = State_Particles_SLTR(:, ind_SLTR) * mu_posterior_part_SLTR';
    
        mu_prior_SLTR = mu_posterior_SLTR;

        RunningTime_PF_SLTR(mc, t_track) = toc + ParticleTransitionTime;

        Pos_SLTR = [
            State_hat_qAvrg_PF_SLTR(1, t_track) .* cos(State_hat_qAvrg_PF_SLTR(3, t_track))
            State_hat_qAvrg_PF_SLTR(1, t_track) .* sin(State_hat_qAvrg_PF_SLTR(3, t_track))
        ];
        MSE_PF_SLTR(mc, t_track)  = norm(State0([1;3], t_track) - Pos_SLTR, 'fro')^2/2;

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
                EvaluateDelayDopplerDoA_CW(GetTau(70, C), 0, 0, Yf, Sf, Delta_f, b, a, Tr, Nc, M, isWindowing);
            end
        end
        % Suppose we have found the maximum-likelihood estiamtes; We even suppose the measurements are unbiased, and with small measurement errors:
        % ---This assumption is very beneficial for KF-ILTR
        y = [GetR(tau0, C); GetVelocity(nu0, lambda); theta0] + chol(diag([1, 0.1, 0.01]), 'lower') * randn(3, 1);
        
        % Kalman updates
        K_KF_ILTR = P_KF_ILTR*H_KF_ILTR'*(H_KF_ILTR*P_KF_ILTR*H_KF_ILTR' + diag([1, 0.1, 0.01]))^-1;
        State_hat_KF_ILTR(:, t_track) = State_hat_KF_ILTR(:, t_track) + K_KF_ILTR*(y - H_KF_ILTR*State_hat_KF_ILTR(:, t_track));
        P_KF_ILTR = (eye(4) - K_KF_ILTR * H_KF_ILTR) * P_KF_ILTR * (eye(4) - K_KF_ILTR * H_KF_ILTR)' + K_KF_ILTR * diag([1, 0.1, 0.01]) * K_KF_ILTR';

        RunningTime_KF_ILTR(mc, t_track) = toc + ParticleTransitionTime;

        Pos_ILTR = [
            State_hat_KF_ILTR(1, t_track) .* cos(State_hat_KF_ILTR(3, t_track))
            State_hat_KF_ILTR(1, t_track) .* sin(State_hat_KF_ILTR(3, t_track))
        ];
        MSE_KF_ILTR(mc, t_track)  = norm(State0([1;3], t_track) - Pos_ILTR, 'fro')^2/2;

        %% PF-SLTR-A
        tic;

        % Calculate Likelihoods
        Utility_SLTR_A = zeros(1, N_par);
        for i = 1:N_par
            R      = State_Particles_SLTR_A(1, i);
            v      = State_Particles_SLTR_A(2, i);
            theta  = State_Particles_SLTR_A(3, i);
            real_beta = State_Particles_SLTR_A(5, i);
            imag_beta = State_Particles_SLTR_A(6, i);
    
            Utility_SLTR_A(i) = 1;

            % Channel Attenuation Coefficient
            AttenuateFactor_Hypothesized  = real_beta + 1j*imag_beta;
            % Spatial Effect
            SpatialEffect_Hypothesized   = b(theta) * a(theta)';
            
            for p = 1:Np
                % Doppler Effect
                DopplerEffect_Hypothesized   = exp(1j*2*pi*GetDoppler(v, lambda)*(p-1)*Tr);
                % Time Delay Effect
                TimeDelayEffect_Hypothesized = repmat(exp(-1j*2*pi*GetTau(R, C)*(0:Nc-1)*Delta_f), Nt, M);     % Every Frequency Component is Phase Shifted
                Error_Hypothesized = Yf{p} - AttenuateFactor_Hypothesized * (SpatialEffect_Hypothesized * (DopplerEffect_Hypothesized * (TimeDelayEffect_Hypothesized .*  Sf{p})));

                Utility_SLTR_A(i) = Utility_SLTR_A(i) * exp(-norm(Error_Hypothesized, 'fro')^2/(2*Pn));
            end
        end

        % Likelihood Weights
        Utility_SLTR_A = (abs(Utility_SLTR_A));
        if max(Utility_SLTR_A) < 1e-8
            Utility_SLTR_A = 1/N_par * ones(1, N_par);
        end
        mu_likelihood_sum_SLTR_A = sum(Utility_SLTR_A);
        mu_likelihood_SLTR_A = Utility_SLTR_A/mu_likelihood_sum_SLTR_A;
    
        % Apply Generalized Bayes Rule
        mu_posterior_SLTR_A = mu_prior_SLTR_A .* mu_likelihood_SLTR_A;
        mu_posterior_sum_SLTR_A = sum(mu_posterior_SLTR_A);
        mu_posterior_SLTR_A = mu_posterior_SLTR_A/mu_posterior_sum_SLTR_A;
    
        % Systematic Resampling
        if 1/(sum(mu_posterior_SLTR_A.^2)) < N_par*0.5   % or 0.25, 0.3 are also good practice 
            [State_Particles_SLTR_A, mu_posterior_SLTR_A] = Resampling(State_Particles_SLTR_A, mu_posterior_SLTR_A);
        end

        mu_quantile_SLTR_A = quantile(mu_posterior_SLTR_A, 0.6);
        ind_SLTR_A = find(mu_posterior_SLTR_A >= mu_quantile_SLTR_A);
        mu_posterior_part_SLTR_A = mu_posterior_SLTR(ind_SLTR_A);
        mu_posterior_part_SLTR_A = mu_posterior_SLTR(ind_SLTR_A)/sum(mu_posterior_part_SLTR_A);
        State_hat_qAvrg_PF_SLTR_A(:, t_track) = State_Particles_SLTR_A(:, ind_SLTR_A) * mu_posterior_part_SLTR_A';
    
        mu_prior_SLTR_A = mu_posterior_SLTR_A;

        RunningTime_PF_SLTR_A(mc, t_track) = toc + ParticleTransitionTime;

        Pos_SLTR_A = [
            State_hat_qAvrg_PF_SLTR_A(1, t_track) .* cos(State_hat_qAvrg_PF_SLTR_A(3, t_track))
            State_hat_qAvrg_PF_SLTR_A(1, t_track) .* sin(State_hat_qAvrg_PF_SLTR_A(3, t_track))
        ];
        MSE_PF_SLTR_A(mc, t_track)  = norm(State0([1;3], t_track) - Pos_SLTR_A, 'fro')^2/2;
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

figure;
plot(State0(1, :), State0(3, :), 'r', 'linewidth', 2);
hold on;
Pos_SLTR_A = [
    State_hat_qAvrg_PF_SLTR_A(1, :) .* cos(State_hat_qAvrg_PF_SLTR_A(3, :))
    State_hat_qAvrg_PF_SLTR_A(1, :) .* sin(State_hat_qAvrg_PF_SLTR_A(3, :))
];
plot(Pos_SLTR_A(1, :), Pos_SLTR_A(2, :), 'b--', 'linewidth', 2);                               
mean(MSE_PF_SLTR_A, 2)
mean(RunningTime_PF_SLTR_A, 2)
% axis equal;
% axis([-50 700 -700 50])
set(gca, 'XAxisLocation', 'top')




