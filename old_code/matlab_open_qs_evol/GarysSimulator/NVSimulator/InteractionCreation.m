%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Setup all interactions %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create NV internal and field interactions
% SEE SETTINGS FOR VALUES

% Zero Field Splitting (Hz)
Hnat.createInteraction('NV',[],'ZFS',diag([0 0 ZFS_NV]));

% Electron Zeeman (Hz/G)
Hnat.createInteraction('NV','BStatic','Zeeman',-gammaE*eye(3));

% Nuclear Quadrupole (Hz)
NitrogenQuad = diag([0 0 -5e6]*2/3); %For N14 (I=1>1/2), measured
Hnat.createInteraction('Nitrogen',[],'ZFS',NitrogenQuad);

%Hyperfine (Hz)
NVHyperfineTensor = diag([2.1 2.1 2.3]*1e6); %Measured, see Gali
Hnat.createInteraction('NV','Nitrogen','Hyperfine', NVHyperfineTensor);

% Nuclear Zeeman (Hz/G) 
Hnat.createInteraction('Nitrogen','BStatic','Zeeman',-gammaN14*eye(3));

%% Create Carbon interactions

%Calculate the hyperfine data
%carbonDipoleTensor(i), NitrogenCarbonTensor(i), CarbonCarbonTensor(i,j)
%values
CarbonTensorCreation; 

if(size(carbonListData,1)>0)
    for i=1:size(carbonListData,1)
        % Nuclear Zeeman (Hz/G)
        Hnat.createInteraction(carbonList(i).Name,'BStatic','Zeeman',-gammaCarbon*eye(3));

        % NV Electron - C Dipolar with Contact interaction
        Hnat.createInteraction('NV',carbonList(i).Name,'Hyperfine', NVCarbonTensor{i});

        % Nitrogen - C Dipolar
        Hnat.createInteraction('Nitrogen',carbonList(i).Name,'Hyperfine', NitrogenCarbonTensor{i});
    end
end

%Carbon-carbon
if(size(carbonListData,1)>1)
    for i=2:size(carbonListData,1)
        for j=1:i-1
            Hnat.createInteraction(carbonList(i).Name,carbonList(j).Name,'Hyperfine', CarbonCarbonTensor{i,j});
        end
    end
end


