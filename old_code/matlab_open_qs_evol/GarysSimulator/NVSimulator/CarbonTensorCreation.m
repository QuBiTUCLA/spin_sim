%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Calculate all dipole constant %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NV - Carbon dipole
%Use hyperfine table from Gali

NVCarbonTensor = cell(size(carbonListData,1));
for i=1:size(carbonListData,1)
    currentCarbon = carbonListData(i,:);
    
    %Contact constant
    a_C = currentCarbon(7)*eye(3);

    %Diagonal dipole matrix 
    D_C = [currentCarbon(8) 0                0
           0                currentCarbon(9) 0
           0                0                currentCarbon(10)];

    %Principal axis to NV axis basis change matrix
    V_C = [currentCarbon(14) currentCarbon(15) currentCarbon(16)
           currentCarbon(17) currentCarbon(18) currentCarbon(19)
           currentCarbon(20) currentCarbon(21) currentCarbon(22)];

    NVCarbonTensor{i} = (V_C*(a_C + D_C)/V_C)*1e6; %Values given in MHz
    
    %Symmetrization
    NVCarbonTensor{i} = (NVCarbonTensor{i} + NVCarbonTensor{i}.')/2;
end

%% Nitrogen - Carbon dipole

NitrogenCarbonTensor = cell(size(carbonListData,1));
for i=1:size(carbonListData,1)
    C1 = carbonListPosition(i,:); 
    N1 = [0 0 -1.676*1e-10];
    
    %distance carbon from nitrogen
    V12 = C1-N1;
    
    %Full dipolar tensor in Hz
    dipCstNC = -mu0*hbar/2*gammaN14*gammaCarbon/(norm(V12)^3);
    
    NitrogenCarbonTensor{i} = dipCstNC*(3*(V12.')*V12/(norm(V12)^2)-eye(3));
end

%% Carbon - Carbon dipole

CarbonCarbonTensor = cell(size(carbonListPosition,1));
if(size(carbonListPosition,1)>1)
    for i=2:size(carbonListPosition,1)
        C1 = carbonListPosition(i,:); 
        for j=1:i-1
            C2 = carbonListPosition(j,:);

            %distance carbon1 - carbon2
            VC1C2 = C1-C2;

            %Full dipolar tensor in Hz
            dipCstCC = -mu0*hbar/2*gammaCarbon^2/(norm(VC1C2)^3);

            CarbonCarbonTensor{i,j} = dipCstCC*(3*(VC1C2.')*VC1C2/(norm(VC1C2)^2)-eye(3));
        end
    end
end

%% NV - Carbon dipole checking (calculated)
% NVCarbonTensor2 = cell(size(carbonListData,1));
% for i=1:size(carbonListData,1)
%     C1 = carbonListPosition(i,:); 
%     
%     %distance carbon from NV
%     V12 = C1;
%     
%     %Full dipolar tensor in Hz
%     dipCstNVC = mu0*hbar/2*gammaE*gammaCarbon/(norm(V12)^3);
%     
%     NVCarbonTensor2{i} = -dipCstNVC*(3*(V12.')*V12/(norm(V12)^2)-eye(3));
% end