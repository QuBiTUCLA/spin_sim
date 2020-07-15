%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Define all parameters %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reset the params structure
params = [];

%Hamiltonian
params.numberOfCarbon = 3; %Number of Carbon
params.trotter = 5; %Trotter number for exponential approximation
params.ampDisc = 2000; %Discretization of the amplitude of the fields (-amp:amp), min 100

%Load Data
ChainHamilt(params.numberOfCarbon);
load('HnatData.mat') %Hamiltonian from simulator, variable HnatData

%MethodsGradient
params.optMethod = 'Conjugate'; %Optimization methods: 'Gradient', 'Conjugate', 'BFGS'
params.optType = 'Unitary'; %Optimization types :'Unitary', 'Transfer'
params.gradientOrder = 1; %Improvement is not assured for increased order

%Pulse params
params.pulseLength = 15000;
params.minTimeStep = 8e-7;
params.maxTimeStep = 8e-7;

%Maximum power
params.maxPower = 8e3; %5e6-1e7 maximum experimental power in Hz, must be positive
params.powerGradPenalty = 0;
params.timeGradPenalty = 0;

%Some parameters for the random guess 
%randPower-How random is the initial pulse [0,1]
params.randomPulseFlag = 1; %Disable random, if 0 sometimes block algorithm, if 1 works will often be best
params.initPulseScale = 1; %in percent [0 1] in the scale from -1 to 1
params.randPower = 1;

%Annealing parameters for global convergence
params.annealingSteps = 15; %Number of randomization, can make algorithm quite long

%initStepSize-Initial stepsize
%minstepInterval-Minimum section interval (if we're not moving anywhere we should stop searching)
params.maxStepSize = 0.1; % 1 is good
params.minStepInterval = 5e-3;

%Desired fidelity for the unitary (this is the trace squared fidelity
params.fidelity = 0.001;
params.improvechk = 0; %Tolerance for improving DIRECTION VARIATION in percent
params.fidelityImp = 5e-5;%Improvement in fidelity (usually bad test)

params.logCount = 1;

%% Controled states, Input and wanted output
Carbon = Spin(1/2);

%Control fields
HControlBath = 0;
for ctC = 1:params.numberOfCarbon
    HControlBath = HControlBath + ...
        HnatData.expandOperator(strcat('Carbon',num2str(ctC)),Carbon.Sx);
end
params.ControlMats{1} = HControlBath;

%Input state
rhoin = diag([1 0]);
for ctC = 1:params.numberOfCarbon-1
    rhoin = kron(rhoin,eye(2)/2);
end
params.rhoin = rhoin/norm(rhoin,'fro');
clear rhoin

% % Operator goal
gate = Gates();
params.opGoal = gate.SWAP(HnatData,'Carbon1',strcat('Carbon',num2str(params.numberOfCarbon)));

%% Plots
%Figure to be shown (projected state evolution)
params.watchSpin{1} = 'Carbon1';
params.watchSpin{2} = strcat('Carbon',num2str(params.numberOfCarbon));