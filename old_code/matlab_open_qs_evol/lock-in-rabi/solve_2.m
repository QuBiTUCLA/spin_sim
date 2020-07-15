function pop = solve_2()


N = 2;

% set the tolerances for the integration
options = odeset('RelTol',1e-12,'AbsTol',ones(1,N)*1e-12);

U_I_C = cell(N,N);
U_I_0 = [];


int_rlimit = 30;

[T,Y] = ode45(@calc_direct_timeprop,[0 int_rlimit],[1  0],options);
T_C = T.';

pop = abs(Y(:,1).').^2;


%U_I_C{1,1} = [Y(:,1).'];
%U_I_C{2,1} = [Y(:,2).'];
%U_I_C{1,2} = [Y(:,3).'];
%U_I_C{2,2} = [Y(:,4).'];
%
%help = zeros(N,N);
%for k = 1:N
%	for l = 1:N
%		help(k,l) = U_I_C{k,l}(length(U_I_C{k,l}));
%	end;
%end;
%
%U_I_0 = help;
%




% integration function


function dy = calc_direct_timeprop(t,y)



w1 = 0;
e0_1 = 1;
mu1 = 1;
eta = 1.6;
Omega = 2.5;

% Optical Bloch Equations

dy = zeros(2,1);    % a column vector

dy(1) = e0_1 * i * exp(-i*w1*t - i*eta*sin(Omega*t)) * mu1 * y(2);
dy(2) = e0_1 * i * exp( i*w1*t + i*eta*sin(Omega*t)) * mu1 * y(1);




