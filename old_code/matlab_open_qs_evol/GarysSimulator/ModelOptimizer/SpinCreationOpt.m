%% Constants
mu0 = 4*pi*1e-7; %H.m-1
hbar = 1.0545715964207855e-34; %rad 
gammaN14 = 3.0766e6; %N14 gyromagnetic value in Hz (known value), e-4 pour T->G
gammaE = -2.8024953e10; %Hz/T
gammaCarbon = 10.705e6; %carbon13 known value, Hz/T
ZFS_NV = 2.87e9; %Hz

%% Create Spins
HnatModel = Hamiltonian();

%Create applied static field
BStatic = Field([0 0 params.B0]); % field in Gauss
BStatic.Name = 'BStatic';
HnatModel.addField(BStatic);

%Create NV and numberOfCarbon carbon 13 
NV = Spin(params.NVSpinNumber);
NV.Name = 'NV';
HnatModel.addSpin(NV);

Nitrogen = Spin(1);
Nitrogen.Name = 'Nitrogen';
HnatModel.addSpin(Nitrogen);

carbonList = cell(params.numberOfCarbon,1);
for ctC = 1:params.numberOfCarbon
    carbonList{ctC} = Spin(1/2);
    carbonList{ctC}.Name = strcat('Carbon',num2str(ctC));
    HnatModel.addSpins(carbonList{ctC});
end

%% Create Interactions
% Zero Field Splitting (Hz)
HnatModel.createInteraction('NV',[],'ZFS',diag([0 0 ZFS_NV]));

% Electron Zeeman (Hz/G)
HnatModel.createInteraction('NV','BStatic','Zeeman',-gammaE*eye(3));

% Nuclear Quadrupole (Hz)
NitrogenQuad = diag([0 0 -5e6]*2/3); %For N14 (I=1>1/2), measured
HnatModel.createInteraction('Nitrogen',[],'ZFS',NitrogenQuad);

% Hyperfine (Hz)
NVHyperfineTensor = diag([2.1 2.1 2.3]*1e6); %Measured, see Gali
HnatModel.createInteraction('NV','Nitrogen','Hyperfine', NVHyperfineTensor);

% Nuclear Zeeman (Hz/G) 
HnatModel.createInteraction('Nitrogen','BStatic','Zeeman',-gammaN14*eye(3));

% Carbon interactions (Tensors to optimize)
% carbonIntTensors = cell(round(params.numberOfCarbon*((params.numberOfCarbon+3)/2)),1);
k=1;

%Carbon-NV/Nitrogen/Static Field
if(params.numberOfCarbon>0)
    for i=1:params.numberOfCarbon
        % Nuclear Zeeman (Hz/G)
        HnatModel.createInteraction(carbonList{i}.Name,'BStatic','Zeeman',-gammaCarbon*eye(3));
        
        %TENSORS MUST BE SYMMETRIC
        % NV Electron - C Dipolar
        HnatModel.createInteraction('NV',carbonList{i}.Name,'Hyperfine', params.NVC*ones(3));
        carbonIntTensors{k} = HnatModel.findObject(strcat('NV',carbonList{i}.Name,'Hyperfine'));
        averageCoupling(k) = params.NVC;
        k=k+1;
                
        if(~params.symmetryFlag)
            % Nitrogen - C Dipolar
            HnatModel.createInteraction('Nitrogen',carbonList{i}.Name,'Hyperfine', params.NiC*ones(3));
            carbonIntTensors{k} = HnatModel.findObject(strcat('Nitrogen',carbonList{i}.Name,'Hyperfine'));
            averageCoupling(k) = params.NiC;
            k=k+1;
        end
    end
end

if(~params.symmetryFlag)
    %Carbon-carbon
    if(params.numberOfCarbon>1)
        for i=2:params.numberOfCarbon
            for j=1:i-1
                HnatModel.createInteraction(carbonList{i}.Name,carbonList{j}.Name,'Hyperfine', params.CC*ones(3));
                carbonIntTensors{k} = HnatModel.findObject(strcat(carbonList{i}.Name,carbonList{j}.Name,'Hyperfine')); %#ok<*SAGROW>
                averageCoupling(k) = params.CC;
                k=k+1;
            end
        end
    end
end
