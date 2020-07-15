function cleancodeGENERALIZEDtest() %does Rabi, Rasmsey and Rot-theta
%built on cleancodeGENERALIZED(), adds the feature that one can have N
%points at each half-cycle of matrix-mult evolution

%Static noise

%Stochastic noise 
%characterized by 0 mean and autocorrelation function 
% = stdev^2*exp(-timeunit/tauC) = variance*exp(-timeunit/tauC)

%in the Detuning frequency Det Delta + dDelta(t), plus in the Rabi

%plotting in time
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Rabi g1 g2 anglerot DabEq;
global Delta

do_Rabi = 1;
do_Rot = 1;
do_Ram = 1;

Rabi = 2*pi*20e6; % bare Rabi in circular freq
g1 = 0; %Right now T1 infinite
g2 = 0; %(5/100*Rabi)^(3/4)*(sqrt(2)*sin(anglerot/2)/anglerot)^(1/2);%*(timeunit_noise_is_cst_Delta/N);        %1/(100e-6);

%For rotary
angle =1/2; %in pi units; change here
anglerot =angle*pi;
Delta=0; %*2*pi*(1)*2e6;   %2.1e6; % By definition, w0 - w, but here parametrized by % of Rabi;

DabEq =0; %all pop in ground state
rhoggEq = (-DabEq + 1)/2;

u_init = 0;
v_init = 0;
w_init = -1;

nbcycle =100; %ceil(Rabi/anglerot/g2); 

%N is the nb of pts per half cycles, ie, 2N is the nb of pts per cycle
N =4;
timeunit = anglerot/Rabi/N;

% Do noise in Rabi frequency?
do_stonoiseR = 0;
do_stanoiseR = 0;

% Do noise in Delta (bath noise)?
do_stonoiseD = 1;
do_stanoiseD = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noise params

%how many times to average over noise, both in Rabi and in Delta
avg_noise_nb_of_times = 100;

% one half cycle has N points; number below is from 1 to ...
timeunit_noise_is_cst_Delta =1; %time duration of piecewise-constant noise in timeunits
timeunit_noise_is_cst_Rabi = 2;  %time duration of piecewise-constant noise in temeunits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOISE IN DELTA
% stoch noise in Delta params: autocorrelation function = variance*exp(-timeunit/tauC)
varianceD = (25/100*Rabi)^2;
tauCD = 200e-9; % 5*pi/Rabi; %   pi/2/Rabi;    %200*timeunit;
stonoiseD=sample_stoch_noise(varianceD,tauCD,timeunit,2*nbcycle*N,avg_noise_nb_of_times);
%way to call it: stonoise(nb_of_noise_avg_param,cycle_time_param)

%static noise in Delta (bath noise)
Distr = 'Normal';
AD = 0; %A param for normal is mean
BD = 0;%(5/100*Rabi); %B param for normal is std dev
stanoiseD = random(Distr,AD,BD,[avg_noise_nb_of_times,2*nbcycle*N]);
% correction for the length of the piecewise constant noise

%code below extends stanoiseD but it's ok bc extended array is not used
for hel = 0:1:ceil(2*N*nbcycle/timeunit_noise_is_cst_Delta)-1
    for hel2 = 1:1:timeunit_noise_is_cst_Delta-1
stanoiseD(:,1 + hel*timeunit_noise_is_cst_Delta + hel2) = stanoiseD(:,1 + hel*timeunit_noise_is_cst_Delta);
    end
end
stanoiseD = stanoiseD(:,1:1:2*nbcycle*N); %cutting extended stanoiseR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOISE IN RABI
% stoch noise in Rabi params: autocorrelation function = variance*exp(-timeunit/tauC)
varianceR =0;%(5/100*Rabi)^2;
tauCR = 0;%200e-9;    %pi/200/Rabi; %70e-9;                   %1*timeunit; %2*timeunit gives echo time
stonoiseR=sample_stoch_noise(varianceR,tauCR,timeunit,2*nbcycle*N,avg_noise_nb_of_times);
%way to call it: stonoise(nb_of_noise_avg_param,cycle_time_param)

%static noise in Rabi (excitation noise)
Distr = 'Normal';
AR = 0; %A param for normal is mean
BR = 0;% (5/100*Rabi); %B param for normal is std dev
stanoiseR = random(Distr,AR,BR,[avg_noise_nb_of_times,2*nbcycle*N]);
%code below extends stanoiseR but it's ok bc extended array is not used
for hel = 0:1:ceil(2*N*nbcycle/timeunit_noise_is_cst_Rabi)-1
    for hel2 = 1:1:timeunit_noise_is_cst_Rabi-1
stanoiseR(:,1 + hel*timeunit_noise_is_cst_Rabi + hel2) = stanoiseR(:,1 + hel*timeunit_noise_is_cst_Rabi);
    end
end
stanoiseR = stanoiseR(:,1:1:2*nbcycle*N); %cutting extended stanoiseR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SUM
noiseD = zeros(size(stonoiseD));
noiseR = zeros(size(stonoiseD));

if do_stonoiseD
    noiseD = stonoiseD;    
end

if do_stanoiseD
    noiseD = noiseD + stanoiseD;    
end


if do_stonoiseR
    noiseR = stonoiseR;    
end

if do_stanoiseR
    noiseR = noiseR + stanoiseR;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To calculate fidelity decay
Auxpop = {};
AuxU = {};
AuxV = {}; 
AuxpopRot = {}; 
AuxURot = {};
AuxVRot = {}; 
AuxpopRam = {}; 
AuxURam = {};
AuxVRam = {}; 

% vectors to average over realizations of noise
Zmmavg = zeros(1,2*nbcycle*N+1);
Zmmrotavg = zeros(1,2*nbcycle*N+1);
Zmmramavg = zeros(1,2*nbcycle*N+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (do_stonoiseR)||(do_stanoiseR)||(do_stonoiseD)||(do_stanoiseD)

    g11 = 0; %g1;
    g22 = 0; %g2;

for aux = 1:1:avg_noise_nb_of_times
    
    aux/avg_noise_nb_of_times

% Matrix multiplication

Zmm = [];
umm = [];
vmm = [];

Zmmrot = [];
ummrot = [];
vmmrot = [];

Zmmram = [];
ummram = [];
vmmram = [];


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
    
    EpsD = noiseD(aux,k);
    EpsR = noiseR(aux,k);
    
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
    
    EpsD = noiseD(aux,kk);
    EpsR = noiseR(aux,kk);
        
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


rho = [1/2,1i/2,-1i/2,1/2]'; %[1,0,0,0] rotated pi/2 around x
Zmmram(1) = rho(1);
ummram(1) = 2*real(rho(2));
vmmram(1) = -2*imag(rho(2));

if do_Ram
%%Ramsey
kk = 0;
for nb_aux = 1/N:1/N:(2*nbcycle)
    
    kk = kk + 1;
    
    EpsD = noiseD(aux,kk);
    EpsR = noiseR(aux,kk);
    
    M = [0,-1i*(EpsR)/2,1i*(EpsR)/2,g11;-1i*(EpsR)/2,1i*(Delta + EpsD) - g22,0,1i*(EpsR)/2;1i*(EpsR)/2, 0,-1i*(Delta+EpsD)-g22,-1i*(EpsR)/2;0,1i*(EpsR)/2,-1i*(EpsR)/2,-g11];
    
    rho = expm(timeunit*M)*rho;
    rhof = 0.5*[rho(1) + 1i*(rho(2) - rho(3)) + rho(4),1i*rho(1) + rho(2) + rho(3) - 1i*rho(4) ,-1i*rho(1) + rho(2) + rho(3) + 1i*rho(4) ,rho(1) - 1i*rho(2) + 1i*rho(3) + rho(4)]'; %final rho rotated pi/2 around x
    Zmmram(kk+1) = rhof(1);
    ummram(kk+1) = 2*real(rhof(2));
    vmmram(kk+1) = -2*imag(rhof(2));
end
end

Zmmramavg = Zmmramavg + Zmmram;
Zmmrotavg = Zmmrotavg + Zmmrot;
Zmmavg = Zmmavg + Zmm;

Auxpop{aux} = Zmm;
AuxU{aux} = umm;
AuxV{aux} = vmm;

AuxpopRot{aux} = Zmmrot;
AuxURot{aux} = ummrot;
AuxVRot{aux} = vmmrot;

AuxpopRam{aux} = Zmmram;
AuxURam{aux} = ummram;
AuxVRam{aux} = vmmram;

end

end

%averaging
Zmmavg = Zmmavg/avg_noise_nb_of_times;
Zmmrotavg = Zmmrotavg/avg_noise_nb_of_times;
Zmmramavg = Zmmramavg/avg_noise_nb_of_times;

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

Zmmram = [];
ummram = [];
vmmram = [];

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


if do_Ram
%%Ramsey

rho = [1/2,1i/2,-1i/2,1/2]'; %[1,0,0,0] rotated pi/2 around x
Zmmram(1) = rho(1);
ummram(1) = 2*real(rho(2));
vmmram(1) = -2*imag(rho(2));

kk = 0;
for nb_aux = 1/N:1/N:(2*nbcycle)
    
    kk = kk + 1;
    
    M = [0,0,0,g1;0,1i*(Delta) - g2,0,0;0, 0,-1i*(Delta)-g2,0;0,0,0,-g1];
    
    rho = expm(timeunit*M)*rho;
    rhof = 0.5*[rho(1) + 1i*(rho(2) - rho(3)) + rho(4),1i*rho(1) + rho(2) + rho(3) - 1i*rho(4) ,-1i*rho(1) + rho(2) + rho(3) + 1i*rho(4) ,rho(1) - 1i*rho(2) + 1i*rho(3) + rho(4)]'; %final rho rotated pi/2 around x
    Zmmram(kk+1) = rhof(1);
    ummram(kk+1) = 2*real(rhof(2));
    vmmram(kk+1) = -2*imag(rhof(2));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_stonoiseR || do_stanoiseR || do_stonoiseD || do_stanoiseD 

%Fidelity decay
%density matrices
    M_Ideal = [];
    MRot_Ideal = [];
    MRam_Ideal = [];
    
    %fidelities
    Fid = [];
    FidRot = [];
    FidRam = [];
    Fidelity = [];
    FidelityRot = [];
    FidelityRam = [];
%     if do_Rabi
%         for aux = 1:1:length(Zmm)
%         M_Ideal(:,:,aux) = [Zmm(aux),(umm(aux)-1i*vmm(aux))/2;(umm(aux)+1i*vmm(aux))/2,1-Zmm(aux)];
%         end
%     end
%     
%     if do_Rot
%         for aux = 1:1:length(Zmmrot)
%         MRot_Ideal(:,:,aux) = [Zmmrot(aux),(ummrot(aux)-1i*vmmrot(aux))/2;(ummrot(aux)+1i*vmmrot(aux))/2,1-Zmmrot(aux)];
%         end
%     end
%     
%     if do_Ram
%         for aux = 1:1:length(Zmmram)
%         MRam_Ideal(:,:,aux) = [Zmmram(aux),(ummram(aux)-1i*vmmram(aux))/2;(ummram(aux)+1i*vmmram(aux))/2,1-Zmmram(aux)];
%         end
%     end
% 
%     for hlp = 1:1:avg_noise_nb_of_times
%        
%         if do_Rabi
%             for aux=1:1:length(Auxpop{hlp})
%             M = [Auxpop{hlp}(aux),(AuxU{hlp}(aux)-1i*AuxV{hlp}(aux))/2;(AuxU{hlp}(aux)+1i*AuxV{hlp}(aux))/2,1-Auxpop{hlp}(aux)];
%             end
%         end            
%             
%         if do_Rot
%             for aux=1:1:length(AuxpopRot{hlp})
%             MRot =  [AuxpopRot{hlp}(aux),(AuxURot{hlp}(aux)-1i*AuxVRot{hlp}(aux))/2;(AuxURot{hlp}(aux)+1i*AuxVRot{hlp}(aux))/2,1-AuxpopRot{hlp}(aux)];
%             end
%         end
%         
%         if do_Ram
%             for aux=1:1:length(AuxpopRam{hlp})
%             MRam =  [AuxpopRam{hlp}(aux),(AuxURam{hlp}(aux)-1i*AuxVRam{hlp}(aux))/2;(AuxURam{hlp}(aux)+1i*AuxVRam{hlp}(aux))/2,1-AuxpopRam{hlp}(aux)];
%             end
%         end
%             
%             %calculate fidelity comparing (ideal) with (non-ideal)
%             %real needed here to correct for calculation mistakes; errors
%             %give imag part o(1e-16)
%             if do_Rabi
%                 for aux=1:1:length(Auxpop{hlp})
%             Fid(hlp,aux) = real(trace(M_Ideal(:,:,aux)*M));
%                 end
%             end
%             if do_Rot
%                 for aux=1:1:length(AuxpopRot{hlp})
%             FidRot(hlp,aux) = real(trace(MRot_Ideal(:,:,aux)*MRot));
%                 end
%             end
%             if do_Ram
%                  for aux=1:1:length(AuxpopRam{hlp})
%             FidRam(hlp,aux) = real(trace(MRam_Ideal(:,:,aux)*MRam));
%                  end
%             end
%        
%         
%     end
%     
%      %avg over noise
%      if do_Rabi
%         for aux = 1:1:length(Zmm)
%             Fidelity(:,aux) = mean(Fid(:,aux));
%         end
%      end
%      
%     if do_Rot
%             for aux = 1:1:length(Zmmrot)
%                 FidelityRot(:,aux) = mean(FidRot(:,aux));
%             end 
%     end
%     
%     if do_Ram
%             for aux = 1:1:length(Zmmram)
%                 FidelityRam(:,aux) = mean(FidRam(:,aux));
%             end 
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
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
    if do_Ram
         for aux = 1:1:length(Zmmram)
         MRam_Ideal(:,:,aux) = [Zmmram(aux),(ummram(aux)-1i*vmmram(aux))/2;(ummram(aux)+1i*vmmram(aux))/2,1-Zmmram(aux)];
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
            
             if do_Ram
             MRam =  [AuxpopRam{hlp}(aux),(AuxURam{hlp}(aux)-1i*AuxVRam{hlp}(aux))/2;(AuxURam{hlp}(aux)+1i*AuxVRam{hlp}(aux))/2,1-AuxpopRam{hlp}(aux)];
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
            
         if do_Ram
         FidRam(hlp,aux) = real(trace(MRam_Ideal(:,:,aux)*MRam));
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
    
    if do_Ram
            for aux = 1:1:length(Zmmram)
                 FidelityRam(:,aux) = mean(FidRam(:,aux));
             end 
     end
    
    
    %%%%%%%%%%%%%%%%%%%%555
    
%    commenting for the time being
%     figure(21); 
%     hold on
%     % if do_Rabi
%     % plot([0:1/2:(length(Zmmrot)-1)/2]/N*2*anglerot/Rabi,Fidelity,'k--');
%     % end
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
%plot(0:1/2/N:(length(Zmmrot)-1)/2/N,real(Zmmrot),'-.bs') %plot in cycles
hold on

%prb formula
%plot(0:1/2/N:(length(Zmmrot)-1)/2/N,prb,'g') %plot in cycles
%hold on

%plot exponential
%T2 = (1/(5/100*Rabi))^(sqrt(2)/2); %(45/64);%(anglerot/sqrt(2)/sin(anglerot/2))^(1/4);
%plot(0:1/2/N:(length(Zmmrot)-1)/2/N,0.5*(1+exp(-(t./T2).^1)),'-.ko') %plot in cycles
%hold on

%plot Paola's formula
%Ps = exp(-varianceD*tauCD*(t + tauCD*(exp(-t/tauCD -1))).*4*sin(anglerot/2)^2/(mod(anglerot,2*pi))^2); %original corrected with mod
%Ps = exp(-varianceD*tauCD*(t + tauCD*(exp(-t/tauCD -1))).*4*sin(anglerot/2)^2/(mod(anglerot,2*pi) + floor(anglerot/(2*pi)))^2); %WRONG
Ps = exp(-varianceD*tauCD*(t + tauCD*(exp(-t/tauCD) -1)).*4*sin(anglerot/2)^2/(anglerot)^2); %P's original
%Psdois = exp(-varianceD*tauCD*(t + tauCD*(exp(-t/tauCD) -1)).*2*pi*sin(anglerot/2)^2/(anglerot)^2); %P's modified - good for very high tau_c * sigma
%Ps = exp(-varianceD*tauCD*(t + tauCD*(exp(-t/tauCD - 1))).*4/(mod(anglerot,2*pi))^2*sin(anglerot/2)^2);%test
%plot([0:1/2/N:(length(Zmmrot)-1)/2/N]*2*anglerot/Rabi,0.5*(1+cos(anglerot/2)^2+sin(anglerot/2)^2*Ps),'-.ko') %plot in cycles
%plot(0:1/2/N:(length(Zmmrot)-1)/2/N,0.5*(1+Ps),'-.bo') %plot in cycles
hold on
%plot(0:1/2/N:(length(Zmmrot)-1)/2/N,0.5*(1++cos(anglerot/2)^2+sin(anglerot/2)^2*Psdois),'-.go') %plot in cycles

%norelaxrot3 =  0.5*(1 + cos(anglerot/2)^2) + 0.5*sin(anglerot/2)^2.*cos((t*pi/anglerot)*(2*Delta/pi)*sin(anglerot/2)).*cos((t*pi/mod(anglerot,2*pi))*Rabi).*Ps; %correct, with mod
%plot(0:1/2/N:(length(Zmmrot)-1)/2/N,norelaxrot3,'-.ko') %plot in cycles

%%%%%%% NOISE IN RABI FREQ
tauc = tauCR;
x = t;
sigma = sqrt(varianceR);
factor = 1/2;
%infidrot = @(theta) 0.5*(1 + exp(-factor*2*tauc^2*sigma^2*(x./tauc + x*Rabi/theta*(exp(-theta/Rabi/tauc) -1) - tanh(theta/2/Rabi/tauc)^2*(x*Rabi/theta*(exp(-theta/Rabi/tauc)+1) + exp(-x./tauc) - 1))));
%plot([0:1/2/N:(length(Zmmrot)-1)/2/N]*2*anglerot/Rabi,infidrot(anglerot),'-.k') %plot in time



plot([0:1/2/N:(length(Zmmrot)-1)/2/N]*2*anglerot/Rabi,0.5*(1+cos(anglerot/2)^2+sin(anglerot/2)^2*Ps),'-.ko') %plot in time
plot([0:1/2:(length(Zmmrot)-1)/2]/N*2*anglerot/Rabi,FidelityRot,'b');
%infid = (varianceR/Rabi^2)*t.^2*Delta^2/8*((2 + anglerot^2*cos(anglerot) - 2*anglerot*sin(anglerot))/anglerot^2);
%plot([0:1/2/N:(length(Zmmrot)-1)/2/N]*2*anglerot/Rabi,1-infid,'-m') %plot in time
hold on

figure(16)
Psram = exp(-varianceR*tauCR*(t + tauCR*(exp(-t/tauCR) -1))); %P's original
plot([0:1/2/N:(length(Zmmrot)-1)/2/N]*2*anglerot/Rabi,0.5*(1+ Psram),'-.m') %plot in time
hold on
plot([0:1/2:(length(Zmmrot)-1)/2]/N*2*anglerot/Rabi,Fidelity,'-m');

lklk

%matrix multiplication result rotary
if  do_stonoiseR || do_stanoiseR || do_stonoiseD || do_stanoiseD
%plot([0:1/2:(length(Zmmrot)-1)/2]/N*2*anglerot/Rabi,real(Zmmrotavg),'b')
%%plot in time
%plot([0:1/2/N:(length(Zmmrotavg)-1)/2/N],real(Zmmrotavg),'r') %plot in cycles
plot([0:1/2/N:(length(Zmmrotavg)-1)/2/N]*2*anglerot/Rabi,real(Zmmrotavg),'r') %plot in time
end

end

if do_Ram
plot(0:1/2/N:(length(Zmmram)-1)/2/N,real(Zmmram),'-.ms') %plot in cycles
hold on   
%matrix multiplication result rotary
if  do_stonoiseR || do_stanoiseR || do_stonoiseD || do_stanoiseD
%plot([0:1/2:(length(Zmmrot)-1)/2]/N*2*anglerot/Rabi,real(Zmmrotavg),'b')
%%plot in time
plot(0:1/2/N:(length(Zmmramavg)-1)/2/N,real(Zmmramavg),'g') %plot in cycles
end
end


if do_Rabi
%subplot(2,1,2)
% mm result rabi without noise
%plot([0:1/2:(length(Zmm)-1)/2]/N*2*anglerot/Rabi,real(Zmm),'r--') %plot in time
%plot(0:1/2/N:(length(Zmm)-1)/2/N,real(Zmm),'r--') %plot in cycles
hold on
%matrix multiplication result rabi
if  do_stonoiseR || do_stanoiseR || do_stonoiseD || do_stanoiseR 
%plot([0:1/2:(length(Zmmavg)-1)/2]/N*2*anglerot/Rabi,real(Zmmavg),'r') %plot in time
%plot([0:1/2/N:(length(Zmmavg)-1)/2/N],real(Zmmavg),'b') %plot in cycles
end

end
hold off

%title('Rotary: pop in 0, blue is rot matrix mult, black is rot PRB approx, red is rabi matrix mult')
title('blue is rot matrix mult, red is rabi matrix mult')
xlabel(sprintf('(2*anglerot)/Rabi cycles, anglerot = %0.3f * pi', angle))

if  do_stonoiseR || do_stanoiseR || do_stonoiseD || do_stanoiseD

figure(210); 
%     if do_Rabi
%     plot([0:1/2:(length(Zmmrot)-1)/2]/N*2*anglerot/Rabi,Fidelity,'r');
%     hold on
%     plot([0:1/2:(length(Zmmrot)-1)/2]/N*2*anglerot/Rabi,real(Zmm).*real(Zmmavg),'g')
%     end
     hold on
     if do_Rot
     plot([0:1/2:(length(Zmmrot)-1)/2]/N*2*anglerot/Rabi,FidelityRot,'b');
     hold on
     plot([0:1/2:(length(Zmmrot)-1)/2]/N*2*anglerot/Rabi,real(Zmmrot).*real(Zmmrotavg),'k')
     end
     
     title(sprintf('Fidelity in time, red is rabi, blue is rotary'));   
end