function DiffDeltas()

clear all

%working; model incorporates detuning Delta and pulse errors (1+Eps)*Rabi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Rabi g1 g2 Delta DabEq Distr A B;

Rabi = 2*pi*10e6; % bare Rabi in Hz
g1 = 0;% 1/100e-3; %100e-3;
g2 = 1/(5e-6/2);%1/1e-6;%1/2e-6; %20e-6;

DabEq = -1; %all pop in ground state
rhoggEq = (-DabEq + 1)/2;

%(Pulse error + Instability Rabi) Eps distributed following:

Distr = 'Normal';
A = 0; %A param for normal is mean
B = 0.01;%0.01;%0.5; %B param for normal is std dev

%how many times to average over noise
avg_noise_nb_of_times =10;

positive_or_negative_det = -1;
%-1 corresponds to w above w0
%+1 corresponds to w below w0


%By definition, Delta = w0 - w, but here parametrized as a % of the
%Rabi frequency
perc = 0; %change here
Delta = positive_or_negative_det*perc*Rabi/100;

plot_nb = 1;
PopulationRotForDer = [];
PopulationForDer = [];

anglerot_init = 1/3; %in pi units
anglerot_step = 1/3; %in pi units
anglerot_final = 5/3; %in pi units

%number of +X/-X ***2pi*** cycles
nbcycle = ceil(Rabi/pi/g2); 

figure(60);
hold on
colors = ['r' 'g' 'b' 'c' 'm' 'y' 'k'];
sign = ['*', '-'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for anglerotnow = anglerot_init:anglerot_step:anglerot_final
    
    timeunit =anglerot_init*pi/Rabi; %timeunit now is the smallest, initial anglerot
   
    popT = {};
    uT = {};
    vT = {};
    
    poprotT = {};
    urotT = {};
    vrotT = {};
    
    for aux = 1:1:avg_noise_nb_of_times
        
        % Matrix multiplication
        
        pop = [];
        u = [];
        v = [];
        
        poprot = [];
        urot = [];
        vrot = [];
        
        u_init = 0;
        v_init = 0;
        w_init = -1;
        
        %with dissipation, column vectors
        rho = [(1-w_init)/2,(u_init + 1i*v_init)/2,(u_init - 1i*v_init)/2,1-(1-w_init)/2]';
        pop(1) = rho(1);
        u(1) = 2*real(rho(2));
        v(1) = -2*imag(rho(2));
        %% Rabi with noise and detuning
        for nb_aux = 1:1:2*nbcycle/anglerot_init
            
            if mod(nb_aux,2) ~= 0
            Eps = random(Distr,A,B);
            end
            
            % dissipation
            M = [0,-1i*Rabi*(1+Eps)/2,1i*Rabi*(1+Eps)/2,g1;-1i*Rabi*(1+Eps)/2,1i*Delta - g2,0,1i*Rabi*(1+Eps)/2;1i*Rabi*(1+Eps)/2, 0,-1i*Delta-g2,-1i*Rabi*(1+Eps)/2;0,1i*Rabi*(1+Eps)/2,-1i*Rabi*(1+Eps)/2,-g1];
            
            rho = expm(timeunit*M)*rho;
            
            pop(end+1) = rho(1);
            u(end+1) = 2*real(rho(2));
            v(end+1) = -2*imag(rho(2));
            
        end
        
        h = 1;
        rho = [(1-w_init)/2,(u_init + 1i*v_init)/2,(u_init - 1i*v_init)/2,1-(1-w_init)/2]';
        poprot(1) = rho(1);
        urot(1) = 2*real(rho(2));
        vrot(1) = -2*imag(rho(2));
        
        %% Rotary with noise and detuning
        for nb_aux = 1:1:2*nbcycle/anglerot_init
            
            if mod(nb_aux,2) ~= 0 %If beginning of cycle
            Eps = random(Distr,A,B);
            end
            
            % dissipation
            M = [0,-1i*Rabi*h*(1+Eps)/2,1i*Rabi*h*(1+Eps)/2,g1;-1i*Rabi*h*(1+Eps)/2,1i*Delta - g2,0,1i*Rabi*h*(1+Eps)/2;1i*Rabi*h*(1+Eps)/2, 0,-1i*Delta-g2,-1i*Rabi*h*(1+Eps)/2;0,1i*Rabi*h*(1+Eps)/2,-1i*Rabi*h*(1+Eps)/2,-g1];
            
            rho = expm(timeunit*M)*rho;
            
            poprot(end+1) = rho(1);
            urot(end+1) = 2*real(rho(2));
            vrot(end+1) = -2*imag(rho(2));
            
           
            if mod(nb_aux,anglerotnow/anglerot_init) == 0
            h = -h;    
            end
            
           
        end
        
        popT{aux} = pop;
        uT{aux} = u;
        vT{aux} = v;
        
        poprotT{aux} = poprot;
        urotT{aux} = urot;
        vrotT{aux} = vrot;
        
    end
        
        %% Ideal Rabi with detuning, no error
        
        popI = [];
        uI = [];
        vI = [];
        
        %with dissipation, column vectors
        rho = [(1-w_init)/2,(u_init + 1i*v_init)/2,(u_init - 1i*v_init)/2,1-(1-w_init)/2]';
        popI(1) = rho(1);
        uI(1) = 2*real(rho(2));
        vI(1) = -2*imag(rho(2));
        %% Rabi with noise and detuning
        for nb_aux = 1:1:2*nbcycle/anglerot_init
            
            % dissipation
            M = [0,-1i*Rabi/2,1i*Rabi/2,g1;-1i*Rabi/2,1i*Delta - g2,0,1i*Rabi/2;1i*Rabi/2, 0,-1i*Delta-g2,-1i*Rabi/2;0,1i*Rabi/2,-1i*Rabi/2,-g1];
            
            rho = expm(timeunit*M)*rho;
            
            popI(end+1) = rho(1);
            uI(end+1) = 2*real(rho(2));
            vI(end+1) = -2*imag(rho(2));
            
        end
        %% Ideal Rotary with detuning, no error
        poprotI = [];
        urotI = [];
        vrotI = [];
        
        h = 1;
        rho = [(1-w_init)/2,(u_init + 1i*v_init)/2,(u_init - 1i*v_init)/2,1-(1-w_init)/2]';
        poprotI(1) = rho(1);
        urotI(1) = 2*real(rho(2));
        vrotI(1) = -2*imag(rho(2));
        for nb_aux = 1:1:2*nbcycle/anglerot_init
           
            
            % dissipation
            M = [0,-1i*Rabi*h/2,1i*Rabi*h/2,g1;-1i*Rabi*h/2,1i*Delta - g2,0,1i*Rabi*h/2;1i*Rabi*h/2, 0,-1i*Delta-g2,-1i*Rabi*h/2;0,1i*Rabi*h/2,-1i*Rabi*h/2,-g1];
            
            rho = expm(timeunit*M)*rho;
            
            poprotI(end+1) = rho(1);
            urotI(end+1) = 2*real(rho(2));
            vrotI(end+1) = -2*imag(rho(2));
            
            if mod(nb_aux,anglerotnow/anglerot_init) == 0
            h = -h;
            end
        end
        
 %saving here the pops for ideal case no noise
 PopulationForDer(plot_nb,:) = popI;
 PopulationRotForDer(plot_nb,:) = poprotI;
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting

 %density matrices
    M_Ideal = [];
    MRot_Ideal = [];
    
    %fidelities
    Fid = [];
    FidRot = [];
    Fidelity = [];
    FidelityRot = [];
    
    for aux = 1:1:length(popI)
        %%%%check here for u, v pm
        M_Ideal(:,:,aux) = [popI(aux),(uI(aux)-1i*vI(aux))/2;(uI(aux)+1i*vI(aux))/2,1-popI(aux)];
        MRot_Ideal(:,:,aux) = [poprotI(aux),(urotI(aux)+1i*vrotI(aux))/2;(urotI(aux)-1i*vrotI(aux))/2,1-poprotI(aux)];
    end

    for hlp = 1:1:avg_noise_nb_of_times
        Auxpop = popT{hlp};
        AuxU = uT{hlp};
        AuxV = vT{hlp};
        
        AuxpopRot = poprotT{hlp};
        AuxURot = urotT{hlp};
        AuxVRot = vrotT{hlp};
        
        for aux=1:1:length(popI)
            %check here for u,v pm
            M = [Auxpop(aux),(AuxU(aux)-1i*AuxV(aux))/2;(AuxU(aux)+1i*AuxV(aux))/2,1-Auxpop(aux)];
            MRot = [AuxpopRot(aux),(AuxURot(aux)+1i*AuxVRot(aux))/2;(AuxURot(aux)-1i*AuxVRot(aux))/2,1-AuxpopRot(aux)];
            
            %calculate fidelity comparing (ideal - with detuning) with (non-ideal eith errors - with
            %detuning)
            %real needed here to correct for calculation mistakes; errors
            %give imag part o(1e-16)
            Fid(hlp,aux) = real(trace(M_Ideal(:,:,aux)*M));
            FidRot(hlp,aux) = real(trace(MRot_Ideal(:,:,aux)*MRot));
        end
        
    end
    
     %avg over noise
    for aux = 1:1:length(popI)
        
        Fidelity(:,aux) = mean(Fid(:,aux));
        FidelityRot(:,aux) = mean(FidRot(:,aux));
        
    end
    
%     figure(plot_nb);
%     plot(0:anglerot_init/2:nbcycle,Fidelity,'r');
%      hold on
%      plot(0:anglerot_init/2:nbcycle,FidelityRot,'b');
%      hold off

%title(sprintf('Fidelity in time, red is rabi, blue is rotary, %0.3f pi anglerot, flag %d',anglerotnow,positive_or_negative_det));

plot(0:anglerot_init/2:nbcycle,Fidelity,['.' colors(plot_nb)]);
plot(0:anglerot_init/2:nbcycle,FidelityRot,['*' colors(plot_nb)]);
    
    FF(plot_nb,:) = Fidelity;
    FFR(plot_nb,:) = FidelityRot;
    
    plot_nb = plot_nb + 1;
    
    
end
hold off

% Plotting surface

[A,B] = meshgrid(0:1/2:(size(FF,2)-1)/2,anglerot_init:anglerot_step:anglerot_final);

figure(1982)
h1 = surf(A,B,FF);
set(h1,'LineStyle','none')
set(h1,'FaceColor',[1 0 0],'FaceAlpha',0.5);
hold on;
h3 = surf(A,B,FFR);
set(h3,'LineStyle','none')
set(h3,'FaceColor',[0 0 1],'FaceAlpha',0.5);
title('Fidelity decay with time and detuning - rabi is red, rotary is blue')
xlabel('Nb Cycles')
ylabel('Detuning in % Rabi')
zlabel('Fidelity')