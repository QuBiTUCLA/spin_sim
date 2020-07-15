function [T, Y] = pi_do_prop(detection_duration, waiting_time, rho_init)

% laser
global sat;
global laser_rabi;

% decay
global carrier_decay;
global cross_decay;
global zero_e_to_meta;
global ones_e_to_meta;
global meta_to_zero_g;
 
% microwave
global Rabi1;
global Rabim1;
 
global D1;
global Dm1;

global N;
% number of levels
N = 7;




% time to start integration
t_begin = 0;



% scanning parameters
%detection_duration = 50e-9;
%waiting_time = 5*15e-9;



% time when laser is switched on
laser_start_point = 0e-9;
% laser saturation parameters
sat = 1;

% microwave; rabi frequencies as observed in the lab
microwave_startpoint = detection_duration + waiting_time;
microwave_freq = 20e6;
microwave_duration = 1./microwave_freq/2; % corresponding pi time

laser_start_point2 = microwave_startpoint+microwave_duration;

% time to stop integration
t_end = 2*detection_duration + waiting_time + microwave_duration;

% the factor of 0.5 comes from the fact that the Hamiltonian in the RWA approximation with a cos(w*t) microwave excitation yields a factor of 0.5
Rabi1 = @(t) 0.5 * 2*pi*microwave_freq * ( (t>=microwave_startpoint) && (t<=(microwave_startpoint+microwave_duration)) );
Rabim1 = @(t) 0.5 * 2*pi*0e6 * t;

% detunings of microwave, positive is red detuned
D1 = 2*pi*0e6;
Dm1 = 2*pi*0e6;


% the initial density matrix, has to be one row vector
%rho_init = zeros(1, N*N);
%rho_init(1) = 0; % zero_g
%rho_init(17) = 1; % one_g

%rho_init(9) = 0; % minus_one_g




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general do not change parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

laser_rabi = @(t) sat*77e6 * ( (t>=laser_start_point && t<=(laser_start_point+detection_duration)) ...
+ (t>=(laser_start_point2) && t<=(laser_start_point2+detection_duration)) );

% decay
carrier_decay = 77e6;
cross_decay = 1.5e6;
zero_e_to_meta = 0e6;
ones_e_to_meta = 30e6;
meta_to_zero_g = 3.3e6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%figure(666);
%t_arr = t_begin:(t_end-t_begin)/1000:t_end;
%Rabi1_arr = [];
%Rabim1_arr = [];
%for k = t_arr
%	Rabi1_arr(end+1) = Rabi1(k);
%	Rabim1_arr(end+1) = Rabim1(k);
%end;
%plot(t_arr, Rabi1_arr, t_arr, Rabim1_arr);
%drawnow;

% set the tolerances for the integration
options = odeset('RelTol',1e-12,'AbsTol',ones(1,N*N)*1e-12);

% this step calls the integration ode45 function
[T,Y] = ode45(@propagate49, [t_begin t_end], rho_init, options);





%pop_zero_g = real(Y(:, 1));
%pop_mone_g = real(Y(:, 9));
%pop_one_g = real(Y(:, 17));
%pop_meta = real(Y(:, 25));
%pop_zero_e = real(Y(:, 33));
%pop_mone_e = real(Y(:, 41));
%pop_one_e = real(Y(:, 49));
%
%figure(1);
%subplot(2,1,1);
%plot(T, pop_zero_g, 'r', T, pop_one_g, 'k', T, pop_mone_g, 'b', ...
%T, pop_meta, 'g', ...
%T, pop_zero_e, 'r--', T, pop_one_e, 'k--', T, pop_mone_e, 'b--');
%
%hold on;
%plot(T, arrayfun(@(t) Rabi1(t)/max(Rabi1(t)), T) - 0.2, 'g-', 'LineWidth', 10);
%plot(T, arrayfun(@(t) Rabim1(t)/max(Rabim1(t)), T) - 0.2, 'r-', 'LineWidth', 10);
%plot(T, arrayfun(@(t) laser_rabi(t)/max(laser_rabi(t)), T) - 0.2, 'k-', 'LineWidth', 10);
%hold off;
%
%legend('|0>', '|1>', '|-1>', '|M>', '|0>_e', '|1>_e', '|-1>_e');
%
%tr = pop_zero_g + pop_one_g + pop_mone_g + pop_meta + pop_zero_e + pop_one_e + pop_mone_e;
%subplot(2,1,2);
%plot(T, tr - 1);
%
%% remember to clear global variables
%

