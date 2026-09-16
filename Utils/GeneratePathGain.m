%{
    Online supplementary materials of the paper titled:
    A New Particle Filter for Target Tracking in MIMO OFDM Integrated Sensing and Communications

    @Author:    Shixiong Wang (s.wang@xjtu.edu.cn; s.wang@u.nus.edu; wsx.gugo@gmail.com)
    @Institute: School of Mathematics and Statistics, Xi'an Jiaotong University
    @Date:      17 November 2025, 27 April 2026, 4 September 2026
    @Home:      https://github.com/Spratm-Asleaf/PF-SLTR-MIMO-OFDM
%}

function beta0_k = GeneratePathGain(R1, Rk, mean)
    r_phase = rand;

    sigma = 0.001 + exprnd(mean);   % normalized RCS fluctuation with unit mean; 0.001 added for numerical stability

    % if sigma < 0.1          % 0.1 censored for numerical stability
    %     sigma = 0.1;
    % end

    % Easy Model 
    % beta0_k = (1 + 0.2*(2*rand - 1)) * exp(1j*2*pi*r_phase);

    % Another Model: Which is more interetable
    beta0_k = (R1/Rk)^2 * sqrt(sigma) * exp(1j*2*pi*r_phase);

    % Yet another Model: Which is physically more understandable: But somehow equivalent to the above
    % −j*4π*(R0k​/λ)
    % beta0_k = (R1/Rk)^2 * sqrt(sigma) * exp(-1j*4*pi*Rk/0.0300);
end