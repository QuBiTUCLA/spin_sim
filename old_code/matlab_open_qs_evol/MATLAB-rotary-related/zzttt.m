function zzttt()
%2pi 0.0322
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Rabi g1 g2 anglerot Delta DabEq aux nb_aux noise;

do_rabi =0;
do_rotary = 1;

Rabi = 2*pi*40e6;        %10e6; % bare Rabi in circular freq
g1 = 0;%1/1e-3; %1/100e-3; %Right now T1 infinite
g2 = 0; %1/1e-6;

%For rotary
anglerot =1/2*pi;
Delta=2*pi*1e6*2;     %2.1e6;        %3/100*Rabi; % By definition, w0 - w, but here parametrized by % of Rabi;

DabEq =0; %all pop in ground state
rhoggEq = (-DabEq + 1)/2;

u_init = 0;
v_init = 0;
w_init = -1;

nbcycle =50;    % ceil(Rabi/anglerot/g2); 
%nbcycle is Rotary cycles
% total time: nbcycle*2*anglerot/Rabi
timeunit =anglerot/Rabi;
period = 2*timeunit;

N = 20; % nb of points per cycle/2  ie nb of points per cycle is 2*N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NUMERICAL

%original
%options = odeset('RelTol',1e-12,'AbsTol',ones(1,3)*1e-12);
%modified to see if improves convergence
options = odeset('RelTol',1e-6,'AbsTol',ones(1,3)*1e-6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
 close all
figure(100);
hold on
if do_rabi
plot([0:1/2/N:nbcycle]*2*anglerot/Rabi,(1-Wnn)/2,'r'); %numerical no noise
plot(t,rabitest,'m');
end

if do_rotary
 
plot([0:1/2/N:nbcycle]*2*anglerot/Rabi,(1-WRotnn)/2,'r'); %USED plot in time

 tunit =  (2*anglerot*nbcycle/Rabi)/(nbcycle*N*2); %1e-9; %make same number per cycle as in the simulation
 tfinal = 2*anglerot*nbcycle/Rabi;
 t = 0:tunit:tfinal;
 %form = (cos(anglerot * Delta *t*(2*pi*0.0322))).*(cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))))/2 + 0.5; %working 
 form = 0.5*sin(anglerot/2)^2*cos(2*Delta*t/anglerot*sin(anglerot/2)).*(cos(2*anglerot/sin(anglerot/2)*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5)))) + 0.5 + 0.5*cos(anglerot/2)^2; %working 
 plot([0:1/2/N:nbcycle],(1-WRotnn)/2,'r',[0:1/2/N:nbcycle],form,'b'); %USED plot in time
% 
% plot([0:1/2/N:nbcycle],(1-WRotnn)/2 - form,'g')
% 
% resgreen = Delta/Rabi/2*(sin(anglerot*Delta*t*(2*pi*0.0322)).* cos(2*anglerot*abs(pi/4+t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))));
% plot([0:1/2/N:nbcycle],resgreen,'m') %,[0:1/2/N:nbcycle],(1-WRotnn)/2 - form - resgreen, 'k')

% timearray = t;
% NFFT = 2^nextpow2(length(timearray));
% Fs = 1/abs(timearray(1)-timearray(2));                    % Sampling frequency
% fftx = Fs/2*linspace(0,1,NFFT/2+1);
% 
% ffty = fft((1-WRotnn)/2 - form-resgreen-mean((1-WRotnn)/2 - form-resgreen),NFFT)/length(timearray);
% ffty = ffty(1:NFFT/2+1);
% 
% figure(345)
% plot(fftx/1e6,abs(ffty), 'g')
l;l;l;l

timearray = t;
NFFT = 2^nextpow2(length(timearray));
Fs = 1/abs(timearray(1)-timearray(2));                    % Sampling frequency
fftx = Fs/2*linspace(0,1,NFFT/2+1);

ffty = fft((1-WRotnn)/2-mean((1-WRotnn)/2),NFFT)/length(timearray);
ffty = ffty(1:NFFT/2+1);

figure(345)
plot(fftx/1e6,abs(ffty), 'g')

end

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
