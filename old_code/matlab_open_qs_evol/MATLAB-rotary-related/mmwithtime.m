function mmwithtime()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Rabi g1 g2 anglerot Delta DabEq aux nb_aux noise;

do_rotary = 1;

Rabi = 2*pi*10e6;        %10e6; % bare Rabi in circular freq
g1 = 0;%1/1e-3; %1/100e-3; %Right now T1 infinite
g2 = 0; %1/1e-6;

DabEq =0; %all pop in ground state
rhoggEq = (-DabEq + 1)/2;

u_init = 0;
v_init = 0;
w_init = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TotalSign = [];
angleinit = 1/4; %in pi units
anglestep = 1/4;
anglefin = 1;
k = 1;
for anglerot = pi*angleinit:pi*anglestep:pi*anglefin
% total time: nbcycle*2*anglerot/Rabi
timeunit =anglerot/Rabi;
period = 2*timeunit;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SignalDeltaNNN = [];
DeltaInit = 1e6;
DeltaFin = Rabi/2/pi;
DeltaStep = 0.1e6;
for Delta = 2*pi*DeltaInit:2*pi*DeltaStep:2*pi*DeltaFin
    
    time = pi*Rabi/4/Delta/sin(anglerot/2);
    nbcycle = ceil(time);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No stochastic noise

% Matrix multiplication

Zmmnn = [];
ummnn = [];
vmmnn = [];

Zmmrotnn = [];
ummrotnn = [];
vmmrotnn = [];

%with dissipation, column vectors
rho = [(1-w_init)/2,(u_init + 1i*v_init)/2,(u_init - 1i*v_init)/2,1-(1-w_init)/2]';
Zmmnn(1) = rho(1);
ummnn(1) = 2*real(rho(2));
vmmnn(1) = -2*imag(rho(2));
%% Usual Rabi
for nb_aux = 1:1:2*nbcycle
    
    % dissipation
    M = [0,-1i*Rabi/2,1i*Rabi/2,g1;-1i*Rabi/2,1i*(Delta) - g2,0,1i*Rabi/2;1i*Rabi/2, 0,-1i*(Delta)-g2,-1i*Rabi/2;0,1i*Rabi/2,-1i*Rabi/2,-g1];
    
    rho = expm(timeunit*M)*rho;
  
    Zmmnn(nb_aux+1) = rho(1);
    ummnn(nb_aux+1) = 2*real(rho(2));
    vmmnn(nb_aux+1) = -2*imag(rho(2));
    
end

h = 1;
rho = [(1-w_init)/2,(u_init + 1i*v_init)/2,(u_init - 1i*v_init)/2,1-(1-w_init)/2]';
Zmmrotnn(1) = rho(1);
ummrotnn(1) = 2*real(rho(2));
vmmrotnn(1) = -2*imag(rho(2));

%% Rotary
for nb_aux = 1:1:2*nbcycle
        
    % dissipation
    M = [0,-1i*Rabi*h/2,1i*Rabi*h/2,g1;-1i*Rabi*h/2,1i*(Delta) - g2,0,1i*Rabi*h/2;1i*Rabi*h/2, 0,-1i*(Delta)-g2,-1i*Rabi*h/2;0,1i*Rabi*h/2,-1i*Rabi*h/2,-g1];
    
    rho = expm(timeunit*M)*rho;
    Zmmrotnn(nb_aux+1) = rho(1);
    ummrotnn(nb_aux+1) = 2*real(rho(2));
    vmmrotnn(nb_aux+1) = -2*imag(rho(2));
   
    h = -h;
    
end

%mm no noise
figure(100)
hold on
plot(0:1/2:(length(Zmmnn)-1)/2,real(Zmmnn),'b')
plot(0:1/2:(length(Zmmrotnn)-1)/2,real(Zmmrotnn),'r')
hold off

end

end
