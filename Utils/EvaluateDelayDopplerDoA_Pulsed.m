%{
    Online supplementary materials of the paper titled:
    A New Particle Filter for Target Tracking in MIMO OFDM Integrated Sensing and Communications

    @Author:    Shixiong Wang (s.wang@xjtu.edu.cn; s.wang@u.nus.edu; wsx.gugo@gmail.com)
    @Institute: School of Mathematics and Statistics, Xi'an Jiaotong University
    @Date:      17 November 2025, 27 April 2026, 4 September 2026
    @Home:      https://github.com/Spratm-Asleaf/PF-SLTR-MIMO-OFDM
%}

function ret = EvaluateDelayDopplerDoA_Pulsed(tau, nu, theta, Ts, Y, S, Lss, Np, Nr, L, StrVecRx, StrVecTx, Tr, isWindowing)
    % [Np, ~] = size(Y);  % Just pass parameters, instead of infer them, to speed up in large-scale algorithmic loops

    ret = 0;
    for p = 1:Np
        ret = ret + EvaluateDelayDoA_Pulsed(tau, theta, Ts, Y{p}, S{p}(:, 1:Lss), Lss, Nr, L, StrVecRx, StrVecTx, isWindowing) * exp(-1j*2*pi*nu*(p-1)*Tr);
    end
end

 