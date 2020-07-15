%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Settings for Model Optimizer %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%Static field
params.B0 = 400*1e-4; %100G ~ 300 MHz for NV, ~ 0.1 MHz for Carbon
params.NVSpinNumber = 3/4; %3/4 strictly for rotating frame

%Propagator
params.timeStep = 1e-8;
params.timeLength = 1000;

%Number of Carbon to optimize
params.numberOfCarbon = 2;
params.carbonDistRadius = [2 3];
params.trotter = 5; %power of 2

%Load Data (To replace real sample to test this program)
% LauncherHamiltonian(params.NVSpinNumber,params.numberOfCarbon,params.carbonDistRadius,params.B0);
load('HnatData.mat') %Hamiltonian from simulator, variable HnatData
params.Resonance = HnatData.Resonance;

%SCGA parameters
params.initPopulation = 15; % > number of parameters, should be at least 10
params.maxGenerations = 30; %Max number of generations
params.goodCost = 0.02; %global stop (best < goodcost)
params.costError = 1e-6; %Local (max var of cost in local simplex)

%Simplification
params.symmetryFlag = true; %Only NV-C if true, otherwise also C-C, N-C

%Signal noise
params.noiseRatio = 0.00; %signal + [-noise/2 , +noise/2], signal=[0 1]

%Number of pulse to try to optimize
params.numberOfTrial = 1;

%Random pulse params
params.initPulseScale = 1e7;
params.LPFilter = 5e8; %low pass filter on random pulse

%Interaction tensors parameters: used for initial values, the
%amplitude will be normalized for each type of interaction. For example,
%interaction between electron and nuclei are very different.
params.valueRange = 4;% ~[-val +val]*averageCoupling (NVC,NiC,...)
                      % corresponds to initial search space, ~*2 in the end
params.NVC = 1e6; %NV-Carbon
params.NiC = 1e1; %Nitrogen-Carbon
params.CC = 1e1; %Carbon-carbon

%% Controled states, Input state
if(params.NVSpinNumber == 3/4)
    NV = Spin(1/2);
else
    NV = Spin(1);
end
Nitrogen = Spin(1);
Carbon = Spin(1/2);

%Control fields
%NV
params.ControlMats{1} = HnatData.expandOperator('NV',NV.Sx);

% %Bath
% gNV2Nratio = -3.0766e6/2.8024953e10; %Gyromagnetic ratio: from NV to N14 type
% gNV2Cratio = -10.705e6/2.8024953e10; %Gyromagnetic ratio: from NV to C13 type
% RF2MWratio = 10; %Experimental setting
% 
% HControlBath = gNV2Nratio*HnatData.expandOperator('Nitrogen',Nitrogen.Sx);
% for ctC = 1:params.numberOfCarbon
%     HControlBath = HControlBath + gNV2Cratio*...
%         HnatData.expandOperator(strcat('Carbon',num2str(ctC)),Carbon.Sx);
% end
% params.ControlMats{2} = RF2MWratio*HControlBath;

%Initial density matrix
rhoin = HnatData.expandOperator('NV',diag([0 1]));
params.rhoin = rhoin/trace(rhoin);
clear rhoin

