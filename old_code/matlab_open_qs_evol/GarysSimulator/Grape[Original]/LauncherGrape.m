%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Launcher: GRAPE methods for NV %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Gary Wolfowicz, March 2011 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Adapted from program by Ryan Colm %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% profile on
%% Initialization
clear all;

%Load Settings
GrapeSettings; %Prepare all settings, get 'params'

%% Initializing propagator

%Prepare propagator
propagator = PropagatorObj();
propagator.TimeStep = params.timeStep;

propagator.Hnat = HnatData.RotatingMatrix;

for ctMat=1:length(params.ControlMats)
    propagator.ControlMats{ctMat} = params.ControlMats{ctMat};
end
propagator.rhoIn = params.rhoin;

%Precalculation for matrix exponential calculation
%If Trotter number is a power of 2, calculations are faster
propagator.precalculation(params.Trotter)

%% Starting Optimization
%Starting Goal
goalTotalTime = 0;
goalStepCount = 0;

bestFidelity.Intensity = 0; %For first loop

%Pulse Length
currentLength = params.minpLength;

%Optimization of pulse length loop
lengthLoopFlag = 1;
while(lengthLoopFlag) %bestFidelity{goalCount}.Intensity < params.fidelity) 
    
disp(sprintf('\nPulse length: %d', currentLength));
propagator.TimeLength = currentLength;

%Create pulse and fidelity
currentPulse = PulseObj([length(propagator.ControlMats) currentLength]);
if(params.randomPulseFlag)
    currentPulse.Pulse = currentPulse.makeRandomPulse(params.randPower,params.initPulseScale);
else
    currentPulse.Pulse = params.initPulseScale*ones(length(propagator.ControlMats), currentLength);
end
currentFidelity = Fidelity([length(propagator.ControlMats) currentLength]);
% currentFidelity = TimeOptFidelity([length(propagator.ControlMats) currentLength]);
currentFidelity.opGoal = params.opGoal;
currentFidelity.GradientOrder = params.gradientOrder;

%Direction optimization object
optDir = OptDirection();

%Initialize some parameters for the main loop
stuckFlag = 0;
stepCount = 0;
nextStepFlag = true;
annealingCount = 0;
annealStepSize = params.initStepSize;

%Calculate gradient (First step only here due to Wolfe calculating gradient and not Golden)
propagator.ControlFields = currentPulse;
[optDir.newGradient,~] = currentFidelity.makeGradient(params.optType,propagator);

while(nextStepFlag)
    stepCount = stepCount + 1;
    tic

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% MAIN PART %%%%%%%%%%%%%%%%%%%%%%%%%%   
    %Power penalty
    maxGrad = max(abs(optDir.newGradient(:)));
    optDir.newGradient = optDir.newGradient - params.powerGradPenalty*(currentPulse.Pulse/params.maxPower)*maxGrad/10;
    
    %Find optimal direction
    directionImprov = optDir.Next(params.optMethod);

    %Find optimal step size
    optDir.deltaPulse = currentPulse.Pulse; %For BFGS
    switch(params.optStep)
        case 'Golden'
            [currentPulse, currentFidelity, currentStepSize] = GoldStepOpt(params, propagator, currentPulse, currentFidelity, optDir);
            
            %Calculate gradient for next step
            propagator.ControlFields = currentPulse;
            [optDir.newGradient,~] = currentFidelity.makeGradient(params.optType,propagator);
        case 'Wolfe'
            [currentPulse, currentFidelity, currentStepSize, optDir] = WolfeStepOpt(params, propagator, currentPulse, currentFidelity, optDir);
        
        otherwise
            ERROR('No step size line search method with this name');
    end
    optDir.deltaPulse = currentPulse.Pulse - optDir.deltaPulse; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Memorize the best pulse
    fidelImprov = currentFidelity.Intensity - bestFidelity.Intensity;
    if(fidelImprov > 0)
        bestPulse = currentPulse;
        bestFidelity = currentFidelity;
    end
    
    %Check if algorithm is stuck, flag=1 for gradient, flag=2 for BFGS and
    %Conjugate
    if(currentStepSize == 0)
        stuckFlag = stuckFlag + 1;
    else
        stuckFlag = 0;
    end

    %Prepare next step conditions
    currentTime=toc;
    goalTotalTime = goalTotalTime + currentTime;
    goalStepCount = goalStepCount + 1;
    disp(sprintf('Fidelity: %0.5g | Step Duration: %0.3g |Step size: %0.5g',currentFidelity.Intensity,currentTime,currentStepSize)); %#ok<*DSPS>

    %Condition for next loop, be careful with the order below of conditions
    nextStepFlag =  directionImprov > params.improvechk;
    if(strcmp(params.optMethod,'Gradient'))
        nextStepFlag = nextStepFlag && stuckFlag < 1;
    else
        nextStepFlag = nextStepFlag && stuckFlag < 2;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Annealing/Randomization for global search
    if(nextStepFlag == false && annealingCount < params.annealingSteps)
        annealingCount = annealingCount + 1;
        disp(sprintf('\nAnnealing step: %d', annealingCount));
        nextStepFlag = true;
        
        %Randomization with decay when close to fidelity = 1
        currentPulse.Pulse = currentPulse.Pulse+ params.maxPower*...
            (1-currentFidelity.Intensity)^4.5*(2*rand(size(currentPulse.Pulse))-1);
        params.initStepSize = annealStepSize/annealingCount;
        
        %Reinitialization of direction and fidelity objects
        optDir = OptDirection();
        currentFidelity.Intensity = 0;
        currentFidelity.OldIntensity = 0;
        
        %First next step
        propagator.ControlFields = currentPulse;
        [optDir.newGradient,~] = currentFidelity.makeGradient(params.optType,propagator);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Fidelity achieved
    nextStepFlag = nextStepFlag && bestFidelity.Intensity < params.fidelity;
    
end %End current pulse optimization
clear currentPulse

currentLength = currentLength + params.pLengthStep;

lengthLoopFlag = currentLength <= params.maxpLength;
lengthLoopFlag = lengthLoopFlag && (currentFidelity.Intensity - bestFidelity.Intensity) > params.lengthImprov;
lengthLoopFlag = lengthLoopFlag && bestFidelity.Intensity <= params.fidelity;
end %End pulse length

disp(sprintf('----> Total computing duration: %0.5g | Step count: %0.5g',goalTotalTime,stepCount));

%% Measurements of spin states
clear Measurement
clf

figure(1)
propagator.ControlFields = bestPulse;
propagator.TimeLength = bestPulse.Size(2);
Measurement = PlotSpins();
for ctSpin=1:length(params.watchSpin)
    Measurement.WatchSpin(params.watchSpin{ctSpin},HnatData);
end
Measurement.Propagate(propagator);
Measurement.Plot();

figure(4)
plot(bestPulse.Pulse','*-')

%% Tests on the final result
%Test of Trotter Expansion + sparsing
propagator.Hnat = HnatData.RotatingMatrix;
bestEnd = propagator.fullPropagation('Unitary');
bestEndReal = propagator.fullPropagation2('Unitary');
TrotterError = log(abs(1-trace(bestEndReal'*bestEnd/length(bestEnd)))^2);
disp(sprintf('\nError of real propagator vs Trotter: %2.5g', TrotterError));

% %Test of spin reducing 1 to 1/2
% %BAD TEST
% %Must use the right carbon radius so that it does not change
% LauncherHamiltonian(1,params.numberOfCarbon,params.carbonDistRadius);
% load('HnatData.mat');
% NV = Spin(1);
% propagator.Hnat = HnatData.RotatingMatrix;
% propagator.ControlMats{1} = HnatData.expandOperator('NV',NV.Sx);
% propagator.precalculation(params.Trotter);
% bestEndOne = propagator.fullPropagation('Unitary');
% SpinReducingError = log(abs(1-trace(bestEndOne(1:length(bestEnd),1:length(bestEnd))'*bestEnd/length(bestEnd)))^2);
% disp(sprintf('\nError of NV 1/2 vs NV 1: %2.5g', SpinReducingError));

% profile off