function ChainHamilt(numberOfSpins)
%% Spin chain
%Parameters
B0 = 150*1e-4;
gammaCarbon = 10.705e6; %carbon13 known value, Hz/T
mu0 = 4*pi*1e-7; %H.m-1
hbar = 1.0545715964207855e-34; %rad
CCdistance = 2e-10;

%Create Hamiltonian
HnatData = Hamiltonian();

%Create applied static field
BStatic = Field([0 0 B0]); % field in Gauss
BStatic.Name = 'BStatic';
HnatData.addField(BStatic);

%Spins
for ctSpin = 1:numberOfSpins
    carbonList(ctSpin) = Spin(1/2);
    carbonList(ctSpin).Name = strcat('Carbon',num2str(ctSpin));
    HnatData.addSpin(carbonList(ctSpin));
end

%Interactions
for ctSpin = 1:numberOfSpins
    % Nuclear Zeeman (Hz/G)
    HnatData.createInteraction(carbonList(ctSpin).Name,'BStatic','Zeeman',-gammaCarbon*eye(3));
    
    if(ctSpin < numberOfSpins)
        C1 = [ctSpin 0 1]*CCdistance;
        C2 = [ctSpin+1 0 1]*CCdistance;

        %distance carbon1 - carbon2
        VC1C2 = C1-C2;

        %Full dipolar tensor in Hz
        dipCstCC = -mu0*hbar/2*gammaCarbon^2/(norm(VC1C2)^3);

        CarbonCarbonTensor = dipCstCC*(3*(VC1C2.')*VC1C2/(norm(VC1C2)^2)-eye(3));
            
        HnatData.createInteraction(carbonList(ctSpin).Name,carbonList(ctSpin+1).Name,'Hyperfine', CarbonCarbonTensor);
    end
end

%Create Hamiltonian matrix
HnatData.createFullHamiltonian();
HnatData.RotatingMatrix = HnatData.Matrix; %time step small enough that secular is false
save('Data\HnatData','HnatData');
end