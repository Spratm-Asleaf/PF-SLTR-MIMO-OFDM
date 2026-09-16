%{
    Online supplementary materials of the paper titled:
    A New Particle Filter for Target Tracking in MIMO OFDM Integrated Sensing and Communications

    @Author:    Shixiong Wang (s.wang@xjtu.edu.cn; s.wang@u.nus.edu; wsx.gugo@gmail.com)
    @Institute: School of Mathematics and Statistics, Xi'an Jiaotong University
    @Date:      17 November 2025, 27 April 2026, 4 September 2026
    @Home:      https://github.com/Spratm-Asleaf/PF-SLTR-MIMO-OFDM
%}

function s = GetTimeShiftedS(s, n)
% Circular Time Shift of the Transmitted Signal "s"
% In radar signal processing, do not try to call this function unless you know what you are doing

    [N, ~] = size(s);
    s = [zeros(N, n), s(1:N, 1:end-n)];
end

