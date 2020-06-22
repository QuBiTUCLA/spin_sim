%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Define all parameters %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reset the params structure
params = [];

%Hamiltonian
params.NVSpinNumber = 3/4; %Spin of NV ; 1 for normal, 3/4 for 1 reduced to 1/2
params.numberOfCarbon = 0; %Number of Carbon
params.carbonDistRadius = [2.472 2.472]; %Carbon distributed randomly inside this radius
params.trotter = 5; %Trotter number for exponential approximation (power of 2)
params.ampDisc = 3000; %Discretization of the amplitude of the fields,dt (-amp:amp)
%NEED TO MAKE AMPDISC DIFFERENT FOR MANY FIELDS

%Load Data
params.B0 = 150*1e-4;
LauncherHamiltonian(3/4,0,[0 2],150*1e-4);
%LauncherHamiltonian(params.NVSpinNumber,params.numberOfCarbon,...
 %   params.carbonDistRadius,params.B0);
load('HnatData.mat') %Hamiltonian from simulator, variable HnatData

%MethodsGradient
params.optMethod = 'Conjugate'; %Optimization methods: 'Gradient', 'Conjugate', 'BFGS'
params.optType = 'Unitary'; %Optimization types :'Unitary', 'Transfer'
params.gradientOrder = 1; %Improvement is not assured for increased order

%Pulse params
params.pulseLength = 454;
params.minTimeStep = 1e-9;
params.maxTimeStep = 1e-9;
%To have non zero time gradient, dt will be random between min and max at
%initial step.

%Maximum power
params.maxPower = 1e7; %5e6-1e7 maximum experimental power in Hz, must be positive
params.powerGradPenalty = 0; %around 1
params.timeGradPenalty = 0;

%Some parameters for the random guess 
params.randomPulseFlag = 1; %Disable random, if 0 sometimes block algorithm, if 1 works will often be best
params.initPulseScale = 0.99; %in percent [0 1] in the scale from -1 to 1
params.randPower = 1; %How random is the initial pulse [0,1], ~frequency

%Annealing parameters for global convergence
params.annealingSteps = 10; %Number of randomization, can make algorithm quite long

%Golden search step search
params.maxStepSize = 0.2; % Smaller is better but slower, if too small may not work
params.minStepInterval = 1e-3; %Minimum section interval in goldStep

%Desired fidelity for the unitary (this is the trace squared fidelity
params.fidelity = 0.001; %INFIDELITY
params.improvechk = 0; %Tolerance for improving DIRECTION VARIATION in percent
params.fidelityImp = 0;%Improvement in fidelity

%LOG
params.logCount = 1; %Print values on command window every logCount

%% Controled states, Input and wanted output
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
% 
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
%     
%Input state
rhoin = HnatData.expandOperator('NV',diag([0 1]));
% rhoin =kron(diag([0 1]),diag([0 1 0]));
params.rhoin = rhoin/trace(rhoin);
clear rhoin

% %Output states (many goal possible)
% rhogoal = HnatData.expandOperator('NV',diag([1 0]));
% params.opGoal = rhogoal/norm(rhogoal,'fro');
% clear rhogoal
% 
% %Input state
% rhoin = kron(kron(diag([1 0]),eye(3)/3),diag([0 1])); %NV/Nitrogen/Carbon
% params.rhoin = rhoin/trace(rhoin);
% clear rhoin

% % Operator goal
gate = Gates();
params.opGoal = gate.Rotation(HnatData,'NV',0,pi,0,0);
% params.opGoal = gate.CNOT(HnatData,'NV','Carbon1');

%% Plots
%Figure to be shown (projected state evolution)
params.watchSpin{1} = 'NV';
% params.watchSpin{2} = 'Nitrogen';
% params.watchSpin{2} = 'Carbon1';



