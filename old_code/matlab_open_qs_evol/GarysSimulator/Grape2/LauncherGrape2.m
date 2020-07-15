%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Launcher: GRAPE methods for NV %%%%%%%%%%%%%%%%
%%%%%%%% Use Time-optimal and amplitude discretization %%%%%%%%%
%%%%%%%%%%%%%%%%%% Gary Wolfowicz, June 2011 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IMPORTANT: THIS SIMULATOR MINIMIZE THE COST FUNCTION SO THE FIDELITY GOES
%TO ZERO = INFIDELITY. THE DISPLAY VALUE DOES NOT TAKE INTO ACCOUNT THE
%PENALTY FUNCTIONS SO IT MIGHT GO DOWN DURING THE COMPUTATION
mex -c mtimesx.c
profile on
%% Initialization
clear all;

disp('Initialization...');

%Load Settings
GrapeSettings2; %Prepare all settings, get 'params'
% ChainSettings %Settings for spin chain

%Checking parameters
if(params.maxTimeStep < params.minTimeStep)
    temp = params.maxTimeStep;
    params.maxTimeStep = params.minTimeStep;
    params.minTimeStep = temp;
    clear temp;
end
if(params.initPulseScale < 0 || params.initPulseScale > 1)
    error('params.initPulseScale must be [0,1]');
end
if(params.maxTimeStep == params.minTimeStep && params.timeGradPenalty ~=0)
    error('Time penalty must be 0 if maxTimeStep = minTimeStep');
end

%% Initializing propagator
%Prepare propagator
propagator = PropagatorObj2();
propagator.PulseLength = params.pulseLength;

%Normalization to get all pulses [-1,1]
if(params.maxTimeStep ~= params.minTimeStep)
    propagator.Hnat = HnatData.RotatingMatrix*(params.maxTimeStep + params.minTimeStep)/2;
    propagator.ControlMats{1} = HnatData.RotatingMatrix*(params.maxTimeStep - params.minTimeStep)/2;

    for ctMat=1:length(params.ControlMats)
        propagator.ControlMats{ctMat+1} = params.ControlMats{ctMat}*...
            params.minTimeStep*params.maxPower;
    end
else
    propagator.Hnat = HnatData.RotatingMatrix*params.maxTimeStep;
    for ctMat=1:length(params.ControlMats)
        propagator.ControlMats{ctMat} = params.ControlMats{ctMat}*...
            params.maxTimeStep*params.maxPower;
    end
end
propagator.rhoIn = params.rhoin;

%Amplitude discretization
propagator.AmpDisc = (-params.ampDisc:1:params.ampDisc)/params.ampDisc;

%Precalculation for matrix exponential calculation
%If Trotter number is a power of 2, calculations are faster
propagator.precalculation(params.trotter)

%% Starting Optimization
%Starting Goal
params.initStepSize = params.maxStepSize;
goalTotalTime = 0;
goalStepCount = 0;

%Create pulse
currentPulse = PulseObj2([length(propagator.ControlMats) params.pulseLength]);

scale = round([(1-params.initPulseScale)*params.ampDisc (1+params.initPulseScale)*params.ampDisc])+1;
if(params.randomPulseFlag)
    currentPulse.Index = currentPulse.makeRandomIndex(params.randPower,scale);
    if(params.maxTimeStep ~= params.minTimeStep) %Start from small length
        currentPulse.Index(1,:) = ceil(currentPulse.Index(1,:)/10);
    end
else
    currentPulse.Index = scale(2)*ones(length(propagator.ControlMats),params.pulseLength);
end

currentPulse.Pulse(:) = propagator.AmpDisc(currentPulse.Index(:));
bestPulse = currentPulse;

%Create fidelity
currentFidelity = Fidelity2();
currentFidelity.opGoal = params.opGoal;
currentFidelity.GradientOrder = params.gradientOrder;

%Overall best fidelity register
bestFidelity = Fidelity2();

%Direction optimization object
optDir = OptDirection2();

%Initialize some parameters for the main loop
stuckFlag = 0;
stepCount = 0;
currentTime = 0;
nextStepFlag = true;
annealingCount = 0;
oldStepSize = params.initStepSize;
timeShowCount = params.logCount;

%Initialize fidelities
propagator.ControlFields = currentPulse;
[propagator.opEnd,propagator.Unitary] = propagator.fullPropagation(params.optType);

[currentFidelity.Intensity currentFidelity.BaseIntensity] = ...
    currentFidelity.makeFidelity(propagator.opEnd,params,currentPulse.Pulse);
currentFidelity.OldIntensity = currentFidelity.Intensity;
bestFidelity.Intensity = currentFidelity.Intensity;
bestFidelity.OldIntensity = currentFidelity.Intensity;

%Calculate gradient for first step
optDir.newGradient = currentFidelity.makeGradient(params.optType, propagator);

while(nextStepFlag)
    currentFidelity.OldIntensity = currentFidelity.Intensity;
    stepCount = stepCount + 1;
    tic

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% MAIN PART %%%%%%%%%%%%%%%%%%%%%%%%%%   
    if(params.maxTimeStep ~= params.minTimeStep)
        %Time penalty
        optDir.newGradient(1,:) = optDir.newGradient(1,:) + ...
            params.timeGradPenalty/params.pulseLength;

        %Power penalty
        for ctField=2:size(optDir.newGradient,1)
            optDir.newGradient(ctField,:) = optDir.newGradient(ctField,:) + ...
                2*params.powerGradPenalty/params.pulseLength*currentPulse.Pulse(ctField,:);
        end
    else
        for ctField=1:size(optDir.newGradient,1)
            optDir.newGradient(ctField,:) = optDir.newGradient(ctField,:) + ...
                2*params.powerGradPenalty/params.pulseLength*currentPulse.Pulse(ctField,:);
        end
    end
    %CHECK AGAIN ABOVE NORMALIZATION
    
    %Find optimal direction
    directionImprov = optDir.Next(params.optMethod);

    %Find optimal step size
    optDir.deltaPulse = currentPulse.Pulse; %For BFGS
    [currentPulse, currentFidelity, currentStepSize] = GoldStepOpt2(params, propagator, currentPulse, currentFidelity, optDir);
    optDir.deltaPulse = currentPulse.Pulse - optDir.deltaPulse; 
    
    %Calculate gradient for next step
    propagator.ControlFields = currentPulse; %use real pulse (not index) for calcul of HTOT
    optDir.newGradient = currentFidelity.makeGradient(params.optType,propagator);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Memorize the best pulse
    fidelImprov = currentFidelity.Intensity - bestFidelity.Intensity;
    if(fidelImprov <= 0)
        bestPulse = currentPulse;
        bestFidelity = currentFidelity;
    end

    %Init step size reevaluation
    params.initStepSize = min(params.maxStepSize,currentStepSize + oldStepSize);
    oldStepSize = currentStepSize;
    
    %Check if algorithm is stuck, flag=1 for gradient, flag=2 for BFGS and
    %Conjugate due to trying gradient when step = 0 with those last methods
    if(currentStepSize == 0)
        stuckFlag = stuckFlag + 1;
    else
        stuckFlag = 0;
    end

    %LOG
    currentTime = currentTime + toc;
    if(timeShowCount >= params.logCount)
        timeShowCount=0;
        goalTotalTime = goalTotalTime + currentTime;
        
        pulseDuration = 1e6*sum((params.maxTimeStep + params.minTimeStep)/2+...
            bestPulse.Pulse(1,:)*(params.maxTimeStep - params.minTimeStep)/2);

        disp(sprintf('Complete: %0.5g %%|Calcul: %0.3g s|Pulse length: %0.5g us', ...
            (1-currentFidelity.BaseIntensity)*100,currentTime,pulseDuration));
        currentTime = 0;
    end
    timeShowCount = timeShowCount + 1;
    goalStepCount = goalStepCount + 1;

    %Condition for next loop, be careful with the order below of conditions
    nextStepFlag =  directionImprov > params.improvechk;
    if(strcmp(params.optMethod,'Gradient'))
        nextStepFlag = nextStepFlag && stuckFlag < 1;
    else
        nextStepFlag = nextStepFlag && stuckFlag < 2;
    end
    
    nextStepFlag = nextStepFlag && params.fidelityImp < abs(currentFidelity.Intensity - currentFidelity.OldIntensity);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Annealing/Randomization for global search
    if(nextStepFlag == false && annealingCount < params.annealingSteps)
        annealingCount = annealingCount + 1;
        disp(sprintf('\nAnnealing step: %d', annealingCount));
        nextStepFlag = true;
        
        params.maxStepSize = params.maxStepSize*0.8;
        
        %Randomization with decay when close to fidelity = 0   
        %time optimization
        currentPulse.Index = currentPulse.Index + round(rand(size(currentPulse.Index))*...
            (0.5+currentFidelity.BaseIntensity^2));%*ceil(params.ampDisc/100*currentFidelity.BaseIntensity);
        currentPulse.Index(:) = max(1,min(2*params.ampDisc+1,currentPulse.Index(:)));

        currentPulse.Pulse(:) = propagator.AmpDisc(currentPulse.Index(:));
        
        %Reinitialization of direction and fidelity objects
        optDir = OptDirection2();
        
%         params.initStepSize = (params.initStepSize + params.maxStepSize/annealingCount)/2;
        
        %Reinitialize fidelities
        propagator.ControlFields = currentPulse;
        [propagator.opEnd,propagator.Unitary] = propagator.fullPropagation(params.optType);
        
        [currentFidelity.Intensity currentFidelity.BaseIntensity] = ...
            currentFidelity.makeFidelity(propagator.opEnd,params,currentPulse.Pulse);
        currentFidelity.OldIntensity = currentFidelity.Intensity;
        bestFidelity.Intensity = currentFidelity.Intensity;
        bestFidelity.OldIntensity = currentFidelity.Intensity;
        
        %Recalculate initial gradient
        optDir.newGradient = currentFidelity.makeGradient(params.optType,propagator);  
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Fidelity achieved (Base intensity for stop criteria, penalty in
    %fidelity not used here
    nextStepFlag = nextStepFlag && bestFidelity.BaseIntensity > params.fidelity;
    
end %End current pulse optimization

disp(sprintf('----> Total computing duration: %0.5g | Step count: %0.5g',goalTotalTime,stepCount));

%% Measurements of spin states
clear Measurement
clf

figure(1)
propagator.ControlFields = bestPulse;

Measurement = PlotSpins2();
for ctSpin=1:length(params.watchSpin)
    Measurement.WatchSpin(params.watchSpin{ctSpin},HnatData);
end
Measurement.Propagate(params, propagator);
Measurement.Plot();

if(params.maxTimeStep ~= params.minTimeStep)
    figure(4)
    time = (params.maxTimeStep + params.minTimeStep)/2 + ...
        (params.maxTimeStep - params.minTimeStep)/2 * bestPulse.Pulse(1,:);

    hold on
    color=['r','g','b','c','m','y','k','w'];
    for ctField = 2:size(bestPulse.Pulse,1)
        amplitude = bestPulse.Pulse(ctField,:)*params.maxPower*params.minTimeStep./time;
        plot(cumsum(time,2)*1e6,amplitude/1e6,color(ctField-1));
    end
    hold off
    xlabel('Pulse duration (us)');
    ylabel('Pulse amplitude (MHz)');

    figure(5)
    plot(time*1e9)
    xlabel('Pulse length');
    ylabel('Time step (ns)');
else
    figure(4)
    plot(bestPulse.Pulse','*-')
end

%% Tests on the final result
%Test of Trotter Expansion
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

profile off