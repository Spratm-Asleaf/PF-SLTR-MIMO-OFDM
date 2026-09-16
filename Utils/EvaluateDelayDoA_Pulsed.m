%{
    Online supplementary materials of the paper titled:
    A New Particle Filter for Target Tracking in MIMO OFDM Integrated Sensing and Communications

    @Author:    Shixiong Wang (s.wang@xjtu.edu.cn; s.wang@u.nus.edu; wsx.gugo@gmail.com)
    @Institute: School of Mathematics and Statistics, Xi'an Jiaotong University
    @Date:      17 November 2025, 27 April 2026, 4 September 2026
    @Home:      https://github.com/Spratm-Asleaf/PF-SLTR-MIMO-OFDM
%}

function ret = EvaluateDelayDoA_Pulsed(tau, theta, Ts, Y, S, Lss, Nr, L, StrVecRx, StrVecTx, isWindowing)
    k = floor(tau/Ts) + 1;
    HypothesizedRx = (StrVecRx(theta) * StrVecTx(theta)') * S(:, 1:Lss);

    buffer = 5;
    HypothesizedRx = [zeros(Nr, buffer), HypothesizedRx, zeros(Nr, buffer)];

    lower = max(1, k-buffer);
    upper = min(k+Lss-1+buffer, L);
    int_region = lower:upper;       % Suppose true target delay is within this window; good practice is buffer = 3, 5
    ret = trace(Y(:, int_region)*HypothesizedRx');
end

