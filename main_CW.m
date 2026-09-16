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
TrajectoryMode = 2;                  % 1: Almost-Straightline/Slightly-Curved, 2: Very-Circled
switch TrajectoryMode
    case 1  
        Tc = 2e-6;                  % CP length: 1/5 length of OFDM symbol duration
    case 2
        Tc = 1e-6;
    otherwise 
        error('main :: Error in TrajectoryMode');
end
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
Tt   = 0.05;          % Real-time tracking interval is 25 ms; in practice, 50 ms is also a good choice 
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
    switch TrajectoryMode
        case 1
            assert(Tc <= 2.01e-6 && Tc >= 1.99e-6);          % Set Tc = 2e-6;
            if t_track == 1
                % Initial true state
                State0(:, t_track) = [
                    10                      % Position in x axis
                    2                        % Velocity in x axis
                    -40                     % Position in y axis
                    2                        % Velocity in y axis
                ];
            else
                Delta_Vx = 0.02*sin(0.1*TimeTrack(t_track))/6;
                Delta_Vy = 0.02*cos(0.1*TimeTrack(t_track))/6;
                State0(:, t_track) = F * State0(:, t_track - 1) + diag([0, 1, 0, 1]) * [0; Delta_Vx; 0; Delta_Vy] + G * chol(matQ_CV, 'lower') * randn(2, 1);
            end

        case 2
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
            
        otherwise 
            error('main :: Error in TrajectoryMode');
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
MonteCarlo = 1;    % For a long trajectory, time-average and MonteCarlo-average make no significant differences
SquaredError   = zeros(MonteCarlo, LenTimeTrack);    
RunningTime    = zeros(MonteCarlo, LenTimeTrack);

PathGainRecord = zeros(1, LenTimeTrack);

for mc = 1:MonteCarlo
    % Estimated State
    State_hat_Avrg  = zeros(4, LenTimeTrack);   % Posterior Average
    State_hat_qAvrg = zeros(4, LenTimeTrack);   % Trimed (q-quantile) Posterior Average
    State_hat_ML    = zeros(4, LenTimeTrack);   % Max Likelihood   (ML)
    State_hat_MAP   = zeros(4, LenTimeTrack);   % Max A-Posteriori (MAP)
    
    % Particles at each time step
    N_par = 200;                               % Number of Particles
    mu_prior = zeros(1, N_par) + 1/N_par;      % Prior Weights
    mu_likelihood = zeros(1, N_par) + 1/N_par; % Likelihood Weights
    State_Particles = zeros(4, N_par);         % State Particles
    
    % Track mode
    TrackMode = 'Range-Theta';                 % 'Position-Velocity'; 'Range-Theta'
    
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
        PathGainRecord(t_track) = AttenuateFactor;

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
        tic;

        % State Particles Transition
        isWindowing = 0;
    
        switch TrackMode
            case 'Position-Velocity'
                if t_track == 1
                    matR = diag([0.1, 0.01, 0.1, 0.01]);
                    State_Particles = State0(:, t_track) + chol(matR, 'lower') * randn(4, N_par);
            
                    % Particle filter parameters initialization
                    switch TrajectoryMode
                        case 1
                            alpha = 2;                                 % Inflation factor of process noise (10 is a very good value)
                        case {2, 3}
                            alpha = 20;                                % Inflation factor of process noise (10 is a very good value)
                        otherwise
                            error('main :: Error in TrajectoryMode');
                    end
                else
                    State_Particles = F * State_Particles + G * chol(alpha * matQ_CV, 'lower') * randn(2, N_par);
                end
    
                % Calculate Likelihoods
                Utility = zeros(1, N_par);
                for i = 1:N_par
                    p_obj = State_Particles([1;3], i);
                    v_obj_cart = State_Particles([2;4], i);
            
                    R      = norm(p_s - p_obj, 2);
                    theta  = atan2(p_obj(2) - p_s(2), p_obj(1) - p_s(1));
                    v      = v_obj_cart(1) * cos(theta) + v_obj_cart(2) * sin(theta);
            
                    Utility(i) = EvaluateDelayDopplerDoA_CW(GetTau(R, C), GetDoppler(v, lambda), theta, Yf, Sf, Delta_f, b, a, Tr, Nc, M, isWindowing);
                end
    
            case 'Range-Theta'
                if t_track == 1
                    matR = diag([0.1, 0.01, deg2rad(0.1), deg2rad(0.01)]);
                    v_cart = State0([2;4], t_track);
                    v_polar = GetPolarVelocity(v_cart, theta0);
                    State_Particles = [GetR(tau0, C); v_polar(1); theta0; v_polar(2)/GetR(tau0, C)] + chol(matR, 'lower') * randn(4, N_par);
                else
                    matQ = diag([20, deg2rad(10)]);
                    State_Particles = F * State_Particles + G * chol(matQ, 'lower') * randn(2, N_par);
                end
    
                % Calculate Likelihoods
                Utility = zeros(1, N_par);
                for i = 1:N_par
                    R      = State_Particles(1, i);
                    v      = State_Particles(2, i);
                    theta  = State_Particles(3, i);
            
                    Utility(i) = EvaluateDelayDopplerDoA_CW(GetTau(R, C), GetDoppler(v, lambda), theta, Yf, Sf, Delta_f, b, a, Tr, Nc, M, isWindowing);
                end
    
            otherwise 
                error('main :: Error in TrackMode');
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
            [State_Particles, mu_posterior] = Resampling(State_Particles, mu_posterior);
        end
    
        % FusionMode = 2;
        % switch FusionMode
            % case 1
                State_hat_Avrg(:, t_track) = State_Particles * mu_posterior';
    
            % case 2
                mu_quantile = quantile(mu_posterior, 0.6);
                ind = find(mu_posterior >= mu_quantile);
                mu_posterior_part = mu_posterior(ind);
                mu_posterior_part = mu_posterior(ind)/sum(mu_posterior_part);
                State_hat_qAvrg(:, t_track) = State_Particles(:, ind) * mu_posterior_part';
    
            % case 3
                [~, ind] = max(mu_likelihood);
                State_hat_ML(:, t_track) = State_Particles(:, ind);
            % case 4
                [~, ind] = max(mu_posterior);
                State_hat_MAP(:, t_track) = State_Particles(:, ind);
            % otherwise
                % error('main :: Error in Fusion Mode');
        % end
    
        % Prior at this time is the posterior for the next time
        mu_prior = mu_posterior;

        RunningTime(mc, t_track) = toc;

        switch TrackMode
            case 'Position-Velocity'
                Pos = State_hat_qAvrg([1;3], t_track);

            case 'Range-Theta'
                Pos = [
                    State_hat_qAvrg(1, t_track) .* cos(State_hat_qAvrg(3, t_track))
                    State_hat_qAvrg(1, t_track) .* sin(State_hat_qAvrg(3, t_track))
                ];
    
            otherwise 
                error('main :: Error in TrackMode');
        end
        SquaredError(mc, t_track)  = norm(State0([1;3], t_track) - Pos, 'fro')^2/2;
    end
end

%% Display Position Tracking Results: MSEs and Plots
disp('=============================');
disp(['Running Times (ms):   ' num2str(mean(RunningTime(1:300), 2))]);
disp('=============================');
disp(['Tracking MSE:         ' num2str(mean(SquaredError(:)))]);
disp('=============================');
switch TrackMode
    case 'Position-Velocity'
        figure;
        plot(State0(1, :), State0(3, :), 'r', 'linewidth', 2);
        hold on;
        plot(State_hat_Avrg(1, :), State_hat_Avrg(3, :), 'b--', 'linewidth', 2);
        disp(['Full PF      : ' num2str(norm(State0([1;3], :) - State_hat_Avrg([1;3], :), 'fro')^2/(2*LenTimeTrack))]);
        % title(['Num Particles: ' num2str(N_par)]);
        % axis equal;
        % axis([-50 700 -700 50])
        set(gca, 'XAxisLocation', 'top')
        title('Full PF');

        figure;
        plot(State0(1, :), State0(3, :), 'r', 'linewidth', 2);
        hold on;
        plot(State_hat_qAvrg(1, :), State_hat_qAvrg(3, :), 'b--', 'linewidth', 2);                                    
        disp(['q-Quantile PF: ' num2str(norm(State0([1;3], :) - State_hat_qAvrg([1;3], :), 'fro')^2/(2*LenTimeTrack))]);
        % axis equal;
        % axis([-50 700 -700 50])
        set(gca, 'XAxisLocation', 'top')
        title('q-Quantile PF');

        figure;
        plot(State0(1, :), State0(3, :), 'r', 'linewidth', 2);
        hold on;
        plot(State_hat_ML(1, :), State_hat_ML(3, :), 'b--', 'linewidth', 2);                                    
        disp(['ML PF        : ' num2str(norm(State0([1;3], :) - State_hat_ML([1;3], :), 'fro')^2/(2*LenTimeTrack))]);
        % axis equal;
        % axis([-50 700 -700 50])
        set(gca, 'XAxisLocation', 'top')
        title('ML PF');

        figure;
        plot(State0(1, :), State0(3, :), 'r', 'linewidth', 2);
        hold on;
        plot(State_hat_MAP(1, :), State_hat_MAP(3, :), 'b--', 'linewidth', 2);
        disp(['MAP PF       : ' num2str(norm(State0([1;3], :) - State_hat_MAP([1;3], :), 'fro')^2/(2*LenTimeTrack))]);
        % axis equal;
        % axis([-50 700 -700 50])
        set(gca, 'XAxisLocation', 'top')
        title('MAP PF');

    case 'Range-Theta'
        figure;
        plot(State0(1, :), State0(3, :), 'r', 'linewidth', 2);
        hold on;
        Pos = [
            State_hat_Avrg(1, :) .* cos(State_hat_Avrg(3, :))
            State_hat_Avrg(1, :) .* sin(State_hat_Avrg(3, :))
        ];
        plot(Pos(1, :), Pos(2, :), 'b--', 'linewidth', 2);
        disp(['Full PF      : ' num2str(norm(State0([1;3], :) - Pos, 'fro')^2/(2*LenTimeTrack))]);
        % title(['Num Particles: ' num2str(N_par)]);
        % axis equal;
        % axis([-50 700 -700 50])
        set(gca, 'XAxisLocation', 'top')
        title('Full PF');

        figure;
        plot(State0(1, :), State0(3, :), 'r', 'linewidth', 2);
        hold on;
        Pos = [
            State_hat_qAvrg(1, :) .* cos(State_hat_qAvrg(3, :))
            State_hat_qAvrg(1, :) .* sin(State_hat_qAvrg(3, :))
        ];
        plot(Pos(1, :), Pos(2, :), 'b--', 'linewidth', 2);                               
        disp(['q-Quantile PF: ' num2str(norm(State0([1;3], :) - Pos, 'fro')^2/(2*LenTimeTrack))]);
        % axis equal;
        % axis([-50 700 -700 50])
        set(gca, 'XAxisLocation', 'top')
        title('q-Quantile PF');

        figure;
        plot(State0(1, :), State0(3, :), 'r', 'linewidth', 2);
        hold on;
        Pos = [
            State_hat_ML(1, :) .* cos(State_hat_ML(3, :))
            State_hat_ML(1, :) .* sin(State_hat_ML(3, :))
        ];
        plot(Pos(1, :), Pos(2, :), 'b--', 'linewidth', 2);                             
        disp(['ML PF        : ' num2str(norm(State0([1;3], :) - Pos, 'fro')^2/(2*LenTimeTrack))]);
        % axis equal;
        % axis([-50 700 -700 50])
        set(gca, 'XAxisLocation', 'top')
        title('ML PF');

        figure;
        plot(State0(1, :), State0(3, :), 'r', 'linewidth', 2);
        hold on;
        Pos = [
            State_hat_MAP(1, :) .* cos(State_hat_MAP(3, :))
            State_hat_MAP(1, :) .* sin(State_hat_MAP(3, :))
        ];
        plot(Pos(1, :), Pos(2, :), 'b--', 'linewidth', 2);
        disp(['MAP PF       : ' num2str(norm(State0([1;3], :) - Pos, 'fro')^2/(2*LenTimeTrack))]);
        % axis equal;
        % axis([-50 700 -700 50])
        set(gca, 'XAxisLocation', 'top')
        title('MAP PF');

    otherwise 
            error('main :: Error in TrackMode');
end

%%
% if SNR >= 10
%     SquaredError = sort(SquaredError);
%     MSEs = mean(SquaredError(1:1885))
% end