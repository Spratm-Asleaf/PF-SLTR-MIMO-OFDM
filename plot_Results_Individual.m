%{
    Online supplementary materials of the paper titled:
    A New Particle Filter for Target Tracking in MIMO OFDM Integrated Sensing and Communications

    @Author:    Shixiong Wang (s.wang@xjtu.edu.cn; s.wang@u.nus.edu; wsx.gugo@gmail.com)
    @Institute: School of Mathematics and Statistics, Xi'an Jiaotong University
    @Date:      17 November 2025, 27 April 2026, 4 September 2026
    @Home:      https://github.com/Spratm-Asleaf/PF-SLTR-MIMO-OFDM
%}


%% Plot Trajectory
plot(State0(1, :), State0(3, :), 'r', 'linewidth', 2);
hold on;
Pos = [
    State_hat_qAvrg(1, :) .* cos(State_hat_qAvrg(3, :))
    State_hat_qAvrg(1, :) .* sin(State_hat_qAvrg(3, :))
];
plot(Pos(1, :), Pos(2, :), 'b--', 'linewidth', 2);          
hold on;
for j = 1:Ns
    plot(Scatters(1, j) * cos(deg2rad(Scatters(2, j))), Scatters(1, j) * sin(deg2rad(Scatters(2, j))), 'mo', 'linewidth', 1.5, 'markersize', 8);
end
norm(State0([1;3], :) - Pos, 'fro')^2/(2*LenTimeTrack)
% axis equal;
% axis([50 150 -50 60]);
% set(gca, 'XAxisLocation', 'top')
set(gca, 'FontSize', 16);
xlabel('Cartesian Position x-Axis');
ylabel('Cartesian Position y-Axis');
legend({'True Trajectory', 'Estimated Trajectory'}, 'FontSize', 16);

%% Plot Range
figure;
plot(TimeTrack, R0, 'r', 'linewidth', 2);
hold on;
plot(TimeTrack, State_hat_qAvrg(1, :), 'b--', 'linewidth', 2);
xlabel('Time (s)');
ylabel('Range (m)');
set(gca, 'FontSize', 16);
legend({'True Range', 'Estimated Range'}, 'FontSize', 16);

%% Plot Radial Velocity
figure;
plot(TimeTrack, V0, 'r', 'linewidth', 2);
hold on;
plot(TimeTrack, State_hat_qAvrg(2, :), 'b--', 'linewidth', 2);
xlabel('Time (s)');
ylabel('Radial Velocity (m/s)');
set(gca, 'FontSize', 16);
legend({'True Radial Velocity', 'Estimated Radial Velocity'}, 'FontSize', 16);

%% Plot Angle
figure;
plot(TimeTrack, rad2deg(Theta0), 'r', 'linewidth', 2);
hold on;
plot(TimeTrack, rad2deg(State_hat_qAvrg(3, :)), 'b--', 'linewidth', 2);
xlabel('Time (s)');
ylabel('DoA (degree)');
set(gca, 'FontSize', 16);
legend({'True Angle', 'Estimated Angle'}, 'FontSize', 16);