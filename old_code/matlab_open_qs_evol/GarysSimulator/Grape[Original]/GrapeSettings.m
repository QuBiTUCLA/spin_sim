%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Define all parameters %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Reset the params structure
params = [];

%Hamiltonian
params.NVSpinNumber = 3/4; %Spin of NV ; 1 for normal, 3/4 for 1 reduced to 1/2
params.numberOfCarbon = 0; %Number of Carbon
params.carbonDistRadius = [0,0]; %Carbon distributed randomly inside this radius
params.Trotter = 5; %Trotter number for exponential approximation

%Load Data
LauncherHamiltonian(params.NVSpinNumber,params.numberOfCarbon,params.carbonDistRadius,110e-4);
load('HnatData.mat') %Hamiltonian from simulator, variable HnatData

%MethodsGradient
params.optMethod = 'Conjugate'; %Optimization methods: 'Gradient', 'Conjugate', 'BFGS'
params.optStep = 'Golden'; %Golden for exact Golden Method, Wolfe for inexact Wolfe-Powell Method
params.optType = 'Transfer'; %Optimization types :'Unitary', 'Transfer'
params.gradientOrder = 1; %Improvement is not assured for increased order

%Length of pulses: sweeping min:step:max
%timeStep-Length of each time step
params.maxpLength = 300;
params.pLengthStep =0;
params.minpLength = 300; %Must be larger than 2
params.timeStep = 1e-9; %of the order of oscillations (1e-9 minimum experimental pulse generator 
%+ RWA resonance freq max)

%Maximum power
params.maxPower = 1e7; %5e6-1e7 maximum experimental power in Hz, must be positive
params.powerGradPenalty = 0; %around 1, amplitude of power penalty in gradient

%Some parameters for the random guess 
%initPulseScale-Max intensity of the random guess for RF pulse
%randPower-How random is the initial pulse [0,1]
params.randomPulseFlag = 1; %Disable random, if 0 sometimes block algorithm, if 1 works will often be best
params.initPulseScale = params.maxPower/5; %Only for random
params.randPower = 1;

%Annealing parameters for global convergence
params.annealingSteps = 6; %Number of randomization, can make algorithm quite long

%initStepSize-Initial stepsize
%minstepInterval-Minimum section interval (if we're not moving anywhere we should stop searching)
params.initStepSize = 0.4; % 1 is good, only for Golden method, not used in Wolfe
params.minStepInterval = 1e-3; %In percent of initStepSize

%Desired fidelity for the unitary (this is the trace squared fidelity
params.fidelity = 0.999;
params.improvechk = 5e-4; %Tolerance for improving DIRECTION VARIATION in percent
params.lengthImprov = -1; %This condition stops the length loop if the fidelity
%starts decreasing more then this value (negative, -1 to disable)

%% Controled states, Input and wanted output
NV = Spin(params.NVSpinNumber);
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
    
%Input state
%rhoin = HnatData.expandOperator('NV',diag([0 1]));
%params.rhoin = rhoin/trace(rhoin);
%clear rhoin

% %Output states (many goal possible)
% rhogoal = HnatData.expandOperator('NV',diag([1 0]));
 rhogoal = kron(diag([1 0]),[0,0,0;0,1,0;0,0,0]);
 params.opGoal = rhogoal/norm(rhogoal,'fro');
 %clear rhogoal

% %Input state
 rhoin = kron(diag([0 1]),[0,0,0;0,1,0;0,0,0]); %NV/Nitrogen/Carbon
 %rhoin = HnatData.expandOperator('NV',diag([0 1]));
 params.rhoin = rhoin/trace(rhoin);
 clear rhoin

% % Operator goal
%gate = Gates();
%params.opGoal = gate.Rotation2(HnatData,'NV',0,pi,0,0);
% params.opGoal = gate.CNOT(HnatData,'NV','Carbon1');

%% Plots
%Figure to be shown (projected state evolution)
params.watchSpin{1} = 'NV';
params.watchSpin{2} = 'Nitrogen';
% params.watchSpin{2} = 'Carbon1';



