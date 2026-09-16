%{
    Online supplementary materials of the paper titled:
    A New Particle Filter for Target Tracking in MIMO OFDM Integrated Sensing and Communications

    @Author:    Shixiong Wang (s.wang@xjtu.edu.cn; s.wang@u.nus.edu; wsx.gugo@gmail.com)
    @Institute: School of Mathematics and Statistics, Xi'an Jiaotong University
    @Date:      17 November 2025, 27 April 2026, 4 September 2026
    @Home:      https://github.com/Spratm-Asleaf/PF-SLTR-MIMO-OFDM
%}

function [particles_new, weights_new] = Resampling(particles, weights)
% Systematic resampling for particle filter
%
% Inputs:
%   particles : [D x N] matrix of N particles, D is the state dimension
%   weights   : [1 x N] vector of normalized weights (sum to 1)
%
% Outputs:
%   particles_new : [D x N] matrix of resampled particles
%   weights_new   : [1 x N] vector of uniform weights (1/N)

    % Ensure weights are normalized
    % weights = weights / sum(weights);
    N = length(weights);
    
    % Step 1: Compute cumulative sum of weights
    cumulative_sum = cumsum(weights);
    
    % Step 2: Generate systematic points
    u0 = rand / N;                          % random offset in [0, 1/N]
    positions = u0 + (0:N-1)/N;             % N evenly spaced points

    % Step 3: Systematic resampling
    indexes = zeros(1, N);
    i = 1; 
    j = 1;
    while i <= N
        if positions(i) < cumulative_sum(j)
            indexes(i) = j;
            i = i + 1;
        else
            j = j + 1;
        end
    end

    % Step 4: Resample particles and assign uniform weights
    particles_new = particles(:, indexes);
    weights_new = ones(1, N) / N;
end
