%{
    Online supplementary materials of the paper titled:
    A New Particle Filter for Target Tracking in MIMO OFDM Integrated Sensing and Communications

    @Author:    Shixiong Wang (s.wang@xjtu.edu.cn; s.wang@u.nus.edu; wsx.gugo@gmail.com)
    @Institute: School of Mathematics and Statistics, Xi'an Jiaotong University
    @Date:      17 November 2025, 27 April 2026, 4 September 2026
    @Home:      https://github.com/Spratm-Asleaf/PF-SLTR-MIMO-OFDM
%}

function ret = EvaluateDelayDoA_CW(tau, theta, Yf, Sf, Delta_f, StrVecRx, StrVecTx, Nc, M, isWindowing)
    [Nr, ~] = size(Yf);
    [Nt, ~] = size(Sf);

    if isWindowing
        % Using windowing technique to suppress sidelobes; In this case, widening main lobe is a sweet case
        ret =  (StrVecRx(theta) .* Windowing(Nr))' * (Yf * (repmat(exp(-1j*2*pi*tau*(0:Nc-1)*Delta_f), Nt, M) .* repmat((Windowing(Nc))', Nt, M) .* Sf)') * (StrVecTx(theta));
    else
        ret =  (StrVecRx(theta)                 )' * (Yf * (repmat(exp(-1j*2*pi*tau*(0:Nc-1)*Delta_f), Nt, M)                                    .* Sf)') * (StrVecTx(theta));
    end
end
