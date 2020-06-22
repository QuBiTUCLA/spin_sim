%function [numerical, my_func] = numericalexactBOE(Rabi_par, fit_par,angle_par)

function [numerical, my_func] = numericalexactBOE(fit_par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Rabi g1 g2 anglerot Delta DabEq aux nb_aux noise;

do_rabi =0;
do_rotary = 1;

Rabi = 2*pi*1e6*20;
%Rabi = Rabi_par;
g1 = 0;%1/1e-3; %1/100e-3; %Right now T1 infinite
g2 = 0; %1/1e-6;

%For rotary
anglerot =5*pi;
%anglerot =angle_par;
Delta=2*pi*1e6*1;     %2.1e6;        %3/100*Rabi; % By definition, w0 - w, but here parametrized by % of Rabi;

DabEq =0; %all pop in ground state
rhoggEq = (-DabEq + 1)/2;

u_init = 0;
v_init = 0;
w_init = -1;

nbcycle =64;    % ceil(Rabi/anglerot/g2); 
%nbcycle is Rotary cycles
% total time: nbcycle*2*anglerot/Rabi
timeunit =anglerot/Rabi;
period = 2*timeunit;

N = 80; % nb of points per cycle/2  ie nb of points per cycle is 2*N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%original
%options = odeset('RelTol',1e-12,'AbsTol',ones(1,3)*1e-12);
%modified to see if improves convergence
options = odeset('RelTol',1e-3,'AbsTol',ones(1,3)*1e-6);

%
%Numerical no noise

Wnn = [];
Unn = [];
Vnn = [];

WRotnn = [];
URotnn = [];
VRotnn = [];

u_new = u_init;
v_new = v_init;
w_new = w_init;


%% Usual Rabi
if do_rabi
for q = 0:1:nbcycle-1
    
[T,Y] = ode45(@Propagate_phase0nn,[q*period (2*q+1)/2*period],[u_new v_new w_new],options);

t0 = q*period:period/2/N:(2*q+1)/2*period;
y0=interp1(T,Y,t0);

if q ~= 0
   Wnn = Wnn(1:1:end-1);
   Unn = Unn(1:1:end-1);
   Vnn = Vnn(1:1:end-1);
end
Wnn = [Wnn y0(:,3).'];
Unn = [Unn y0(:,1).'];
Vnn = [Vnn y0(:,2).'];

u_new = y0(end,1);
v_new = y0(end,2);
w_new = y0(end,3);

[TT,YY] = ode45(@Propagate_phase0nn,[(2*q+1)/2*period (q+1)*period],[u_new v_new w_new],options);

tt0 = (2*q+1)/2*period:period/2/N:(q+1)*period;

yy0=interp1(TT,YY,tt0);

Wnn = Wnn(1:1:end-1);
Unn = Unn(1:1:end-1);
Vnn = Vnn(1:1:end-1);

Wnn = [Wnn yy0(:,3).'];
Unn = [Unn yy0(:,1).'];
Vnn = [Vnn yy0(:,2).'];

u_new = yy0(end,1);
v_new = yy0(end,2);
w_new = yy0(end,3);

end
end



%% Rotary
if do_rotary
u_new = 0;
v_new = 0;
w_new = -1;
for q = 0:1:nbcycle-1
    
[T,Y] = ode45(@Propagate_phase0nn,[q*period (2*q+1)/2*period],[u_new v_new w_new],options);

t0 = q*period:period/2/N:(2*q+1)/2*period;
y0=interp1(T,Y,t0);

if q ~= 0
   WRotnn = WRotnn(1:1:end-1);
   URotnn = URotnn(1:1:end-1);
   VRotnn = VRotnn(1:1:end-1);
end
WRotnn = [WRotnn y0(:,3).'];
URotnn = [URotnn y0(:,1).'];
VRotnn = [VRotnn y0(:,2).'];

u_new = y0(end,1);
v_new = y0(end,2);
w_new = y0(end,3);

[TT,YY] = ode45(@Propagate_phase180nn,[(2*q+1)/2*period (q+1)*period],[u_new v_new w_new],options);

tt0 = (2*q+1)/2*period:period/2/N:(q+1)*period;

yy0=interp1(TT,YY,tt0);

WRotnn = WRotnn(1:1:end-1);
URotnn = URotnn(1:1:end-1);
VRotnn = VRotnn(1:1:end-1);

WRotnn = [WRotnn yy0(:,3).'];
URotnn = [URotnn yy0(:,1).'];
VRotnn = [VRotnn yy0(:,2).'];

u_new = yy0(end,1);
v_new = yy0(end,2);
w_new = yy0(end,3);
 
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FITS

%%%%% in time
tunit =  (2*anglerot*nbcycle/Rabi)/(nbcycle*N*2); %1e-9; %make same number per cycle as in the simulation
tfinal = 2*anglerot*nbcycle/Rabi;
t = 0:tunit:tfinal;

%Paola's round formula
Proundrot =  0.5*(1 + cos(anglerot/2)^2) + 0.5*sin(anglerot/2)^2.*cos(4*Delta/Rabi*round(Rabi*t/2/anglerot)).*cos(t*Rabi); %ORIGINAL
roundrot =  0.5*(1 + cos(anglerot/2)^2) + 0.5*sin(anglerot/2)^2.*(cos(4*Delta/Rabi*floor(Rabi*t/2/anglerot + 1/2)).*cos(t*Rabi/2).^2); %test

%addi =  0.5*sin(anglerot/2)^2.*(-cos(4*Delta/Rabi*floor(Rabi*t/2/anglerot + 32)).*sin(t*Rabi/2).^2) + 0.5; %works for Rabi/Delta = 20;
%addi =  0.5*sin(anglerot/2)^2.*(cos(4*Delta/Rabi*floor(-Rabi*t/2/anglerot + 16)).*sin(t*Rabi/2).^2) + 0.5; %TRY NOT WORKING Rabi/Delta = 20;
%addi =  0.5*sin(anglerot/2)^2.*(cos(4*Delta/Rabi*floor(-Rabi*t/2/anglerot + 24)).*sin(t*Rabi/2).^2) + 0.5; %works for Rabi/Delta = 30;
%%%%%%%%%%addi =  0.5*sin(anglerot/2)^2.*(cos(4*Delta/Rabi*floor(Rabi*t/2/anglerot - 31)).*sin(t*Rabi/2).^2) + 0.5; %works for Rabi/Delta = 40;
%addi =  0.5*sin(anglerot/2)^2.*(cos(4*Delta/Rabi*floor(-Rabi*t/2/anglerot + 32)).*sin(t*Rabi/2).^2) + 0.5; %NEW TRY Rabi/Delta = 40;
%addi =  0.5*sin(anglerot/2)^2.*(cos(4*Delta/Rabi*floor(-Rabi*t/2/anglerot +40)).*sin(t*Rabi/2).^2) + 0.5; %works for Rabi/Delta = 50;


%addi =  0.5*sin(anglerot/2)^2.*(cos(4*Delta/Rabi*floor(-Rabi*t/2/anglerot + fit_par)).*sin(t*Rabi/2).^2) + 0.5; %TRY NOT WORKING Rabi/Delta = 20;
%addi =  0.5*sin(anglerot/2)^2.*(-cos(4*Delta/Rabi*floor(Rabi*t/2/anglerot + fit_par)).*sin(t*Rabi/2).^2) + 0.5; %TRY NOT WORKING Rabi/Delta = 20;



%original
addi =  0.5*sin(anglerot/2)^2.*(-cos(4*Delta/Rabi*floor(Rabi*t/2/anglerot) + 4*Delta/Rabi*fit_par).*sin(t*Rabi/2).^2) + 0.5; %TRY NOT WORKING Rabi/Delta = 20;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting

%figure(100);
hold on
if do_rabi
%%plot(0:1/20:nbcycle,Zavg,'r');  %numerical with noise
plot([0:1/2/N:nbcycle]*2*anglerot/Rabi,(1-Wnn)/2,'r'); %numerical no noise
plot(t,rabitest,'m');
end

if do_rotary
%plot(0:1/20:nbcycle,Zrotavg,'b');


%plot([0:1/2/N:nbcycle]*2*anglerot/Rabi,(1-WRotnn)/2,'r'); %USED plot in time



%plot([0:1/20:nbcycle],(1-WRotnn)/2,'b'); %USED plot in cycles
%models
%plot([0:1/20:nbcycle]*2*anglerot/Rabi,norelaxrot,'+k');

%plot(t,norelaxrot3,'g+');
%plot(t,relaxrot3,'g+');


%plot(t,Proundrot,'g')
%plot(t,roundrot + addi - 0.5,'b*')
%plot(t,signal,'g')
diffe = (1-WRotnn)/2 - roundrot + 0.5; % difference num and roundrot
%plot(t/(2*anglerot/Rabi),diffe, 'r')
%plot(t/(2*anglerot/Rabi),addi,'k*')
%plot(t/(2*anglerot/Rabi),diffe-addi,'g')

numerical = diffe;
my_func = addi;

%plot([0:1/20:nbcycle]*2*anglerot/Rabi,prb,'k'); %USED plot in time
%plot(t/(2*anglerot/Rabi),norelaxrottest,'m');

%test plot time of maximum signal, approximately 
%n=0:1:2;
%plot(n*pi*anglerot/2/Delta/sin(anglerot/2),zeros(1,length(n)),'s')
%coseno = cos(t*2*Delta/anglerot*sin(anglerot/2))/2 + 0.5;
%plot(t,coseno,'r')
end

hold off
%title('pop in 0, red is rabi, blue is rotary, line is num, points is mm')
%title('Red: Numerics, Green: Analy'); %, Green: Without mod (original form)')
%title('Red numerics, blue new formula')

% %mm
% plot(0:1/2:(length(Zmmavg)-1)/2,real(Zmmavg),'sr')
% plot(0:1/2:(length(Zmmrotavg)-1)/2,real(Zmmrotavg),'sb')
% 
% %mm no noise
% plot(0:1/2:(length(Zmmnn)-1)/2,real(Zmmnn),'dk')
% plot(0:1/2:(length(Zmmrotnn)-1)/2,real(Zmmrotnn),'dk')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% integration function
function dy = Propagate_phase0(t,y)

global Rabi g1 g2 Delta DabEq aux nb_aux noise;

Eps =noise(aux,nb_aux);

% Optical Bloch Equations

dy = zeros(3,1);    % a column vector

%y(1) u
%y(2) v
%y(3) w

dy(1) =  +y(2)*(Delta+Eps) -y(1)*g2;
dy(2) =  -y(1)*(Delta+Eps) -y(3)*Rabi - y(2)*g2;
dy(3) = -(y(3) - DabEq)*g1 + Rabi*y(2);

function dy = Propagate_phase180(t,y)

global Rabi g1 g2 Delta DabEq aux nb_aux noise;

Eps =noise(aux,nb_aux);

% Optical Bloch Equations

dy = zeros(3,1);    % a column vector

%y(1) u
%y(2) v
%y(3) w

dy(1) =  +y(2)*(Delta+Eps) -y(1)*g2;
dy(2) =  -y(1)*(Delta+Eps) +y(3)*Rabi - y(2)*g2;
dy(3) = -(y(3) - DabEq)*g1 - Rabi*y(2);

% integration function
function dy = Propagate_phase0nn(t,y)

global Rabi g1 g2 Delta DabEq aux nb_aux noise;

% Optical Bloch Equations

dy = zeros(3,1);    % a column vector

%y(1) u
%y(2) v
%y(3) w

dy(1) =  +y(2)*(Delta) -y(1)*g2;
dy(2) =  -y(1)*(Delta) -y(3)*Rabi - y(2)*g2;
dy(3) = -(y(3) - DabEq)*g1 + Rabi*y(2);

function dy = Propagate_phase180nn(t,y)

global Rabi g1 g2 Delta DabEq aux nb_aux noise;

% Optical Bloch Equations

dy = zeros(3,1);    % a column vector

%y(1) u
%y(2) v
%y(3) w

dy(1) =  +y(2)*(Delta) -y(1)*g2;
dy(2) =  -y(1)*(Delta) +y(3)*Rabi - y(2)*g2;
dy(3) = -(y(3) - DabEq)*g1 - Rabi*y(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

