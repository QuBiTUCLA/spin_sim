function cleancodeGENERALIZEDtestCSTNOISEinEACHAVG() 
%built on cleancodeGENERALIZED(), adds the feature that one can have N
%points at each half-cycle of matrix-mult evolution

%Static noise

%plotting in time
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Rabi g1 g2 anglerot DabEq;
global Delta

do_Rabi = 1;
do_Rot = 1;

Rabi = 2*pi*20e6; % bare Rabi in circular freq

noise_strengthD = 0/100*Rabi;
noise_strengthR = 5/100*Rabi;

g1 = 0; %Right now T1 infinite
g2 = 0; % (noise_strengthD)^(1)*(sqrt(2)*sin(anglerot/2)/anglerot)^(1);%*(timeunit_noise_is_cst_Delta/N);        %1/(100e-6);

%For rotary
angle = 5; %in pi units; change here
anglerot =angle*pi;
Delta=2*pi*(1)*1e6;   %2.1e6; % By definition, w0 - w, but here parametrized by % of Rabi;

DabEq =0; %all pop in ground state
rhoggEq = (-DabEq + 1)/2;

u_init = 0;
v_init = 0;
w_init = -1;

nbcycle =20; %ceil(Rabi/anglerot/g2); 

%N is the nb of pts per half cycles, ie, 2N is the nb of pts per cycle
N =4;
timeunit = anglerot/Rabi/N;

% Do noise in Rabi frequency?
do_stanoiseR = 1;

% Do noise in Delta (bath noise)?
do_stanoiseD = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noise params

%how many times to average over noise, both in Rabi and in Delta
avg_noise_nb_of_times = 500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%static noise in Delta (bath noise)
DistrD = 'Normal';
AD = 0; %A param for normal is mean
BD = noise_strengthD; %B param for normal is std dev

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%static noise in Rabi (excitation noise)
DistrR = 'Normal';
AR = 0; %A param for normal is mean
BR = noise_strengthR; %B param for normal is std dev
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To calculate fidelity decay
Auxpop = {};
AuxU = {};
AuxV = {}; 
AuxpopRot = {}; 
AuxURot = {};
AuxVRot = {}; 

% vectors to average over realizations of noise
Zmmavg = zeros(1,2*nbcycle*N+1);
Zmmrotavg = zeros(1,2*nbcycle*N+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (do_stanoiseR)||(do_stanoiseD)

    g11 = 0; %g1;
    g22 = 0; %g2;

for aux = 1:1:avg_noise_nb_of_times
    
    aux/avg_noise_nb_of_times
    noiseD = random(DistrD,AD,BD);
    noiseR = random(DistrR,AR,BR);

% Matrix multiplication

Zmm = [];
umm = [];
vmm = [];

Zmmrot = [];
ummrot = [];
vmmrot = [];


%with dissipation, column vectors
rho = [(1-w_init)/2,(u_init + 1i*v_init)/2,(u_init - 1i*v_init)/2,1-(1-w_init)/2]';
Zmm(1) = rho(1);
umm(1) = 2*real(rho(2));
vmm(1) = -2*imag(rho(2));
%% Usual Rabi
if do_Rabi
    k = 0;
for nb_aux = 1/N:1/N:(2*nbcycle)
    
    k = k + 1;
    
    EpsD = noiseD;
    EpsR = noiseR;
    
    % dissipation
    M = [0,-1i*(Rabi+EpsR)/2,1i*(Rabi+EpsR)/2,g11;-1i*(Rabi+EpsR)/2,1i*(Delta + EpsD) - g22,0,1i*(Rabi+EpsR)/2;1i*(Rabi+EpsR)/2, 0,-1i*(Delta+EpsD)-g22,-1i*(Rabi+EpsR)/2;0,1i*(Rabi+EpsR)/2,-1i*(Rabi+EpsR)/2,-g11];
    
    rho = expm(timeunit*M)*rho;
  
    Zmm(k+1) = rho(1);
    umm(k+1) = 2*real(rho(2));
    vmm(k+1) = -2*imag(rho(2));
    
end
end

h = 1;
rho = [(1-w_init)/2,(u_init + 1i*v_init)/2,(u_init - 1i*v_init)/2,1-(1-w_init)/2]';
Zmmrot(1) = rho(1);
ummrot(1) = 2*real(rho(2));
vmmrot(1) = -2*imag(rho(2));


if do_Rot
%% Rotary
kk = 0;
for nb_aux = 1/N:1/N:(2*nbcycle)
    
    kk = kk + 1;
    
    EpsD = noiseD;
    EpsR = noiseR;
        
    % dissipation
    M = [0,-1i*(Rabi+EpsR)*h/2,1i*(Rabi+EpsR)*h/2,g11;-1i*(Rabi+EpsR)*h/2,1i*(Delta+EpsD) - g22,0,1i*(Rabi+EpsR)*h/2;1i*(Rabi+EpsR)*h/2, 0,-1i*(Delta+EpsD)-g22,-1i*(Rabi+EpsR)*h/2;0,1i*(Rabi+EpsR)*h/2,-1i*(Rabi+EpsR)*h/2,-g11];
    
    rho = expm(timeunit*M)*rho;
    Zmmrot(kk+1) = rho(1);
    ummrot(kk+1) = 2*real(rho(2));
    vmmrot(kk+1) = -2*imag(rho(2));
   
    if rem(kk,N) == 0
    h = -h;
    end
    
end
end

Zmmrotavg = Zmmrotavg + Zmmrot;
Zmmavg = Zmmavg + Zmm;

Auxpop{aux} = Zmm;
AuxU{aux} = umm;
AuxV{aux} = vmm;

AuxpopRot{aux} = Zmmrot;
AuxURot{aux} = ummrot;
AuxVRot{aux} = vmmrot;

end
end

%averaging
Zmmavg = Zmmavg/avg_noise_nb_of_times;
Zmmrotavg = Zmmrotavg/avg_noise_nb_of_times;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No stochastic noise

% Matrix multiplication

Zmm = [];
umm = [];
vmm = [];

Zmmrot = [];
ummrot = [];
vmmrot = [];

if do_Rabi
%with dissipation, column vectors
rho = [(1-w_init)/2,(u_init + 1i*v_init)/2,(u_init - 1i*v_init)/2,1-(1-w_init)/2]';
Zmm(1) = rho(1);
umm(1) = 2*real(rho(2));
vmm(1) = -2*imag(rho(2));
%% Usual Rabi
k = 0;
for nb_aux = 1/N:1/N:(2*nbcycle)
    
    k = k + 1;
    
    % dissipation
    M = [0,-1i*Rabi/2,1i*Rabi/2,g1;-1i*Rabi/2,1i*(Delta) - g2,0,1i*Rabi/2;1i*Rabi/2, 0,-1i*(Delta)-g2,-1i*Rabi/2;0,1i*Rabi/2,-1i*Rabi/2,-g1];
    
    rho = expm(timeunit*M)*rho;
  
    Zmm(k+1) = rho(1);
    umm(k+1) = 2*real(rho(2));
    vmm(k+1) = -2*imag(rho(2));
    
end
end

if do_Rot
h = 1;
rho = [(1-w_init)/2,(u_init + 1i*v_init)/2,(u_init - 1i*v_init)/2,1-(1-w_init)/2]';
Zmmrot(1) = rho(1);
ummrot(1) = 2*real(rho(2));
vmmrot(1) = -2*imag(rho(2));

%% Rotary
kk = 0;
for nb_aux = 1/N:1/N:2*nbcycle
    
    kk = kk + 1;
        
    % dissipation
    M = [0,-1i*Rabi*h/2,1i*Rabi*h/2,g1;-1i*Rabi*h/2,1i*(Delta) - g2,0,1i*Rabi*h/2;1i*Rabi*h/2, 0,-1i*(Delta)-g2,-1i*Rabi*h/2;0,1i*Rabi*h/2,-1i*Rabi*h/2,-g1];
    
    rho = expm(timeunit*M)*rho;
    Zmmrot(kk+1) = rho(1);
    ummrot(kk+1) = 2*real(rho(2));
    vmmrot(kk+1) = -2*imag(rho(2));
   
    if rem(kk,N) == 0
    h = -h;
    end
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_stanoiseR || do_stanoiseD 

%Fidelity decay
%density matrices
    M_Ideal = [];
    MRot_Ideal = [];
    
    %fidelities
    Fid = [];
    FidRot = [];
    Fidelity = [];
    FidelityRot = [];
    if do_Rabi
        for aux = 1:1:length(Zmm)
        M_Ideal(:,:,aux) = [Zmm(aux),(umm(aux)-1i*vmm(aux))/2;(umm(aux)+1i*vmm(aux))/2,1-Zmm(aux)];
        end
    end
    
    if do_Rot
        for aux = 1:1:length(Zmmrot)
        MRot_Ideal(:,:,aux) = [Zmmrot(aux),(ummrot(aux)-1i*vmmrot(aux))/2;(ummrot(aux)+1i*vmmrot(aux))/2,1-Zmmrot(aux)];
        end
    end

    for hlp = 1:1:avg_noise_nb_of_times
       
        for aux=1:1:length(Auxpop{hlp})
            if do_Rabi
            M = [Auxpop{hlp}(aux),(AuxU{hlp}(aux)-1i*AuxV{hlp}(aux))/2;(AuxU{hlp}(aux)+1i*AuxV{hlp}(aux))/2,1-Auxpop{hlp}(aux)];
            end
            if do_Rot
            MRot =  [AuxpopRot{hlp}(aux),(AuxURot{hlp}(aux)-1i*AuxVRot{hlp}(aux))/2;(AuxURot{hlp}(aux)+1i*AuxVRot{hlp}(aux))/2,1-AuxpopRot{hlp}(aux)];
            end
            
            %calculate fidelity comparing (ideal) with (non-ideal)
            %real needed here to correct for calculation mistakes; errors
            %give imag part o(1e-16)
            if do_Rabi
            Fid(hlp,aux) = real(trace(M_Ideal(:,:,aux)*M));
            end
            if do_Rot
            FidRot(hlp,aux) = real(trace(MRot_Ideal(:,:,aux)*MRot));
            end
        end
        
    end
    
     %avg over noise
     if do_Rabi
        for aux = 1:1:length(Zmm)
            Fidelity(:,aux) = mean(Fid(:,aux));
        end
     end
     
    if do_Rot
            for aux = 1:1:length(Zmmrot)
                FidelityRot(:,aux) = mean(FidRot(:,aux));
            end 
    end
    
    %commenting for the time being
%     figure(21); 
%     hold on
%      if do_Rabi
%      plot([0:1/2:(length(Zmmrot)-1)/2]/N*2*anglerot/Rabi,Fidelity,'k--');
%      end
%      hold on
%      if do_Rot
%      plot([0:1/2:(length(Zmmrot)-1)/2]/N*2*anglerot/Rabi,FidelityRot,'b--');
%      end
%      hold off
%      title(sprintf('Fidelity in time, red is rabi, blue is rotary'));   

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting for rotary

%PRB modified results for rotary
tunit =  (2*anglerot*nbcycle/Rabi)/(nbcycle*N*2); %1e-9; %make same number per cycle as in the simulation
tfinal = 2*anglerot*nbcycle/Rabi;
t = 0:tunit:tfinal;
c = (3*g1 + 2*g2)/4;
C = sqrt(1 + c^2/Rabi^2);
theta = atan(Rabi/c);
%prb = 0.5*(sin(anglerot/2)^2*C*sin(t*Rabi*pi/anglerot + theta).*cos(2*Delta*t/anglerot*sin(anglerot/2)).*exp(-(c*t)*(3/2)*pi/anglerot*sin(anglerot/2))) + 0.5*(1 + cos(anglerot/2)^2*exp(-(c*t)*anglerot/pi*sin(anglerot/2)));
prb = 0.5*(sin(anglerot/2)^2*C*sin(t*Rabi*pi/anglerot + theta).*cos(2*Delta*t/anglerot*sin(anglerot/2)).*exp(-t*g2/2*(3/2))) + 0.5*(cos(anglerot/2)^2*exp(-t*g2/2*(2/3)) + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting

figure(11); 

hold on
if do_Rot
% %to save rotary data
% tempon3 = [0:1/2:(length(Zmmrot)-1)/2]/N*2*anglerot/Rabi;
% signaln3 = real(Zmmrotavg);
% save('noisy3.mat','tempon3','signaln3')

% mm result rotary without noise
%plot([0:1/2:(length(Zmmrot)-1)/2]/N*2*anglerot/Rabi,real(Zmmrot),'b') %plot in time

%subplot(2,1,1)
plot([0:1/2/N:(length(Zmmrot)-1)/2/N]*2*anglerot/Rabi,real(Zmmrot),'-.bs') %plot in cycles
hold on

%prb formula
%plot(0:1/2/N:(length(Zmmrot)-1)/2/N,prb,'g') %plot in cycles
%hold on

%plot exponential
%T2 = (1/noise_strengthD)^(1);%*(anglerot/sqrt(2)/sin(anglerot/2))^(1);
T2 = (1/noise_strengthD)^(0.93);%*(anglerot/sqrt(2)/sin(anglerot/2))^(1);%(45/64), magic exponent

%plot(0:1/2/N:(length(Zmmrot)-1)/2/N,0.5*(1+exp(-(t./T2).^2)),'-.ko') %plot in cycles
hold on

%matrix multiplication result rotary
if  do_stanoiseR || do_stanoiseD
%plot([0:1/2:(length(Zmmrot)-1)/2]/N*2*anglerot/Rabi,real(Zmmrotavg),'b')
%%plot in time
plot(0:1/2/N*2*anglerot/Rabi:(length(Zmmrotavg)-1)/2/N*2*anglerot/Rabi,real(Zmmrotavg),'r') %plot in cycles
hold on
end

end


if do_Rabi
%subplot(2,1,2)
% mm result rabi without noise
%plot([0:1/2:(length(Zmm)-1)/2]/N*2*anglerot/Rabi,real(Zmm),'r--') %plot in time
%plot(0:1/2/N:(length(Zmm)-1)/2/N,real(Zmm),'b--') %plot in cycles
hold on
%matrix multiplication result rabi
if do_stanoiseR || do_stanoiseD 
%plot([0:1/2:(length(Zmmavg)-1)/2]/N*2*anglerot/Rabi,real(Zmmavg),'r') %plot in time
%plot([0:1/2/N:(length(Zmmavg)-1)/2/N]*2*anglerot/Rabi,real(Zmmavg),'b') %plot in cycles
end

end
hold off

%title('Rotary: pop in 0, blue is rot matrix mult, black is rot PRB approx, red is rabi matrix mult')
title('blue is rot matrix mult, red is rabi matrix mult')
xlabel(sprintf('(2*anglerot)/Rabi cycles, anglerot = %0.3f * pi', angle))

if  do_stanoiseR || do_stanoiseD

figure(210); 
    if do_Rabi
    plot([0:1/2:(length(Zmmrot)-1)/2]/N*2*anglerot/Rabi,Fidelity,'r');
    hold on
    plot([0:1/2:(length(Zmmrot)-1)/2]/N*2*anglerot/Rabi,real(Zmm).*real(Zmmavg),'g')
    end
     hold on
     if do_Rot
     plot([0:1/2:(length(Zmmrot)-1)/2]/N*2*anglerot/Rabi,FidelityRot,'b');
     hold on
     plot([0:1/2:(length(Zmmrot)-1)/2]/N*2*anglerot/Rabi,real(Zmmrot).*real(Zmmrotavg),'k')
     end
     
     title(sprintf('Fidelity in time, red is rabi, blue is rotary'));   
end