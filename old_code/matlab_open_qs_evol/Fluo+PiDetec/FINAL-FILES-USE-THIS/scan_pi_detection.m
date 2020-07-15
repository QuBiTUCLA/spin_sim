close all;
clear all;



detection_duration_arr = [10:10:500]*1e-9;
waiting_time_arr = [0:1:50]*1e-9;

% the initial density matrix, has to be one row vector
N = 7;
rho_init = zeros(1, N*N);


no = 0;
ind = 0;
for tau = detection_duration_arr
	results = [];
	ind = ind + 1;
	for kappa = waiting_time_arr
		no = no + 1;
		no/(length(detection_duration_arr)*length(waiting_time_arr))

		results(end+1).detection_duration = tau;
		results(end).waiting_time = kappa;

		rho_init = zeros(1, N*N);
		rho_init(1) = 1; % zero_g
		[T, Y] = pi_do_prop(tau, kappa, rho_init);

		results(end).init0.T = T;
		results(end).init0.Y = Y;
	
		rho_init = zeros(1, N*N);
		rho_init(17) = 1; % one_g
		[T, Y] = pi_do_prop(tau, kappa, rho_init);

		results(end).init1.T = T;
		results(end).init1.Y = Y;

	end;
	save(['simu_results_' num2str(ind) '.mat'], 'results');
end;


asd




pop_zero_g = real(Y(:, 1));
pop_mone_g = real(Y(:, 9));
pop_one_g = real(Y(:, 17));
pop_meta = real(Y(:, 25));
pop_zero_e = real(Y(:, 33));
pop_mone_e = real(Y(:, 41));
pop_one_e = real(Y(:, 49));

figure(1);
subplot(2,1,1);
plot(T, pop_zero_g, 'r', T, pop_one_g, 'k', T, pop_mone_g, 'b', ...
T, pop_meta, 'g', ...
T, pop_zero_e, 'r--', T, pop_one_e, 'k--', T, pop_mone_e, 'b--');

hold on;
plot(T, arrayfun(@(t) Rabi1(t)/max(Rabi1(t)), T) - 0.2, 'g-', 'LineWidth', 10);
plot(T, arrayfun(@(t) Rabim1(t)/max(Rabim1(t)), T) - 0.2, 'r-', 'LineWidth', 10);
plot(T, arrayfun(@(t) laser_rabi(t)/max(laser_rabi(t)), T) - 0.2, 'k-', 'LineWidth', 10);
hold off;

legend('|0>', '|1>', '|-1>', '|M>', '|0>_e', '|1>_e', '|-1>_e');

tr = pop_zero_g + pop_one_g + pop_mone_g + pop_meta + pop_zero_e + pop_one_e + pop_mone_e;
subplot(2,1,2);
plot(T, tr - 1);

% remember to clear global variables


