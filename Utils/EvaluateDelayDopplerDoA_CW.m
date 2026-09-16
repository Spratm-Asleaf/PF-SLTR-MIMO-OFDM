%{
    Online supplementary materials of the paper titled:
    A New Particle Filter for Target Tracking in MIMO OFDM Integrated Sensing and Communications

    @Author:    Shixiong Wang (s.wang@xjtu.edu.cn; s.wang@u.nus.edu; wsx.gugo@gmail.com)
    @Institute: School of Mathematics and Statistics, Xi'an Jiaotong University
    @Date:      17 November 2025, 27 April 2026, 4 September 2026
    @Home:      https://github.com/Spratm-Asleaf/PF-SLTR-MIMO-OFDM
%}

function ret = EvaluateDelayDopplerDoA_CW(tau, nu, theta, Yf, Sf, Delta_f, StrVecRx, StrVecTx, Tr, Nc, M, isWindowing)
    [Np, ~] = size(Yf);

    ret = 0;
    for p = 1:Np
        ret = ret + EvaluateDelayDoA_CW(tau, theta, Yf{p}, Sf{p}, Delta_f, StrVecRx, StrVecTx, Nc, M, isWindowing) * exp(-1j*2*pi*nu*(p-1)*Tr);
    end
end

 