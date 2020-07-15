function U_I_0 = direct_int_timeprop_2lvl(atomparameters, laserparameters, detunings_I_1, detunings_I_2)

load physics_constants;

global mu_1;
global mu_2;
global det1;
global det2;
global e0_1;
global e0_2;
global shift_pulse;
global tp_1;
global tp_2;
global my_hbar;

det1 		= detunings_I_1;
det2 		= detunings_I_2;
mu_2 		= atomparameters.mu_2;
e0_2 		= laserparameters.e0_2;
tp_2 		= laserparameters.tp_2;
mu_1 		= atomparameters.mu_1;
e0_1 		= laserparameters.e0_1;
tp_1 		= laserparameters.tp_1;
shift_pulse = laserparameters.shift_pulse;
my_hbar 	= hbar;

int_rlimit 	= laserparameters.right_integration_limit; % limit where the integration is stopped is defined

N 		= atomparameters.N*atomparameters.no_of_trap_levels;

% set the tolerances for the integration
options = odeset('RelTol',laserparameters.RelTol,'AbsTol',ones(1,N*N)*laserparameters.AbsTol);

U_I_C = cell(N,N);
U_I_0 = [];


[T,Y] = ode45(@calc_direct_timeprop,[0 int_rlimit],[1  0  0  1],options);
T_C = T.';
U_I_C{1,1} = [Y(:,1).'];
U_I_C{2,1} = [Y(:,2).'];
U_I_C{1,2} = [Y(:,3).'];
U_I_C{2,2} = [Y(:,4).'];

help = zeros(N,N);
for k = 1:N
	for l = 1:N
		help(k,l) = U_I_C{k,l}(length(U_I_C{k,l}));
	end;
end;

U_I_0 = help;

clear mu_1;
clear mu_2;
clear det1;
clear det2;
clear e0_1;
clear e0_2;
clear shift_pulse;
clear tp_1;
clear tp_2;





% integration function


function dy = calc_direct_timeprop(t,y)

global mu_1;
global mu_2;
global det1;
global det2;
global e0_1;
global e0_2;
global shift_pulse;
global tp_1;
global tp_2;
global my_hbar;

% Optical Bloch Equations

dy = zeros(4,1);    % a column vector

dy(1) = ( 0                    + e0_1 * sech(1.763*( t-shift_pulse )/tp_1) * i * exp(-i*det1(1,2)*t) * mu_1(1,2) * y(2));
dy(2) = ( e0_1 * sech(1.763*( t-shift_pulse )/tp_1) * i * exp( i*det1(1,2)*t) * mu_1(1,2)' * y(1) + 0                   );
dy(3) = ( 0                    + e0_1 * sech(1.763*( t-shift_pulse )/tp_1) * i * exp(-i*det1(1,2)*t) * mu_1(1,2) * y(4));
dy(4) = ( e0_1 * sech(1.763*( t-shift_pulse )/tp_1) * i * exp( i*det1(1,2)*t) * mu_1(1,2)' * y(3) + 0                   );




