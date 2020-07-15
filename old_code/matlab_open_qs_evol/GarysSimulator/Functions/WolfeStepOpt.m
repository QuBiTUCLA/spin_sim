%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Function to get best step size optimized for Wolfe Method %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bestPulse, bestFidelity, bestStepSize, optDir] = WolfeStepOpt(params, propagator, currentPulse, currentFidelity, optDir)

%% Parameters
%[Goldstein minimal fidelity update condition]
%Size of the fidelity interval acceptable for an alpha
%Value must be 0 < rho < 1/2 with 1/2 equivalent to having interval = zero
%In truth, due to Goldstein upper limit replaced here by Wolfe condition,
%Interval will not be zero even with rho = 1/2, but condition still
%necessary
rho = 0.1;

%[Wolfe condition]
%Reduction in gradient between two step. Low sigma means high reduction in
%gradient, ie harder condition to respect but better convergence
%Value must be rho < sigma < 1 (strict inequalities)
sigma = 0.4;

%[Interval finder]
%Speed of increase of upper limit in first phase of line search
%Value must be expandSpeed > 1
expandSpeed = 9;

%[Step size finder]
%Speed of collapse of the interval to find a good step size
%Low for lower boundary and up for upper boundary
%Value must be 0 < collapseSpeedLow < collapseSpeedUp < 1/2
collapseSpeedLow = 0.1;
collapseSpeedUp = 0.5;

%Minimum increase in fidelity
epsFid = 5e-5;

%% Initialization
%Best solution
bestPulse = currentPulse;
bestFidelity = currentFidelity;
bestStepSize = 0;

%Last change on pulse direction
maxDir = max(abs(optDir.newDirection(:)));
if(maxDir == 0 || isnan(maxDir) || isinf(maxDir)) %Singular matrix warning problem to really solve
    pulseDirection = 0;
else
    pulseDirection = (params.maxPower/10)*(optDir.newDirection/maxDir);
end

%Gradient of fidelity reduced to single dimension problem of stepSize
grad0 = optDir.newGradient(:).'*pulseDirection(:);
oldGrad = grad0;
oldFidelity = currentFidelity.Intensity;

%Upper limit for interval: make Goldstein limit lower than the maximum
maxStep = (1 - currentFidelity.Intensity)/(rho*grad0);
interval = [0 maxStep];

%Initial guess for step size
oldStep = 0;

newStep = (2*max(currentFidelity.Intensity-...
    currentFidelity.OldIntensity,10*epsFid)/grad0 + maxStep)/2;

bestFidelity.OldIntensity = currentFidelity.Intensity;
    
%% Fidelity and Gradient calculation method
%We calculate both at the same time as main cost is unitary calculation
%and full propagation, both being done in gradient
function [newPulse, newFidelity, newGradient] = updatePulseFidelity(stepSize)
    %New Pulse
    newPulse = currentPulse;
    newPulse.Pulse = currentPulse.Pulse + stepSize*pulseDirection;

    %Limit power
    maxPulse = max(abs(newPulse.Pulse(:)));
    if(maxPulse > params.maxPower && maxPulse ~= 0)
         newPulse.Pulse(:) = newPulse.Pulse(:)*params.maxPower/maxPulse;
    elseif(maxPulse == 0)
        ERROR('Pulse is equal to zero');
    end
    
    %Make new gradient (and propagate)
    propagator.ControlFields = newPulse;
    [newGradient, opEnd] = currentFidelity.makeGradient(params.optType,propagator);
    
    %Check fidelity
    newFidelity = currentFidelity.makeFidelity(opEnd,params.optType);        
end

%% Start optimization
%Optimization possible only for positive gradient
if(grad0 <= 0)
    disp('Gradient <= 0, cannot optimize step');
    return;
end

%First phase: loop find an acceptable interval to find the best step
while(true)
    %Update pulse, fidelity and gradient
    [newPulse, newFidelity, newGradient] = updatePulseFidelity(newStep);
    
    %Calculate gradient for 'reduced' fidelity
    newGrad = newGradient(:).'*pulseDirection(:);
    
    %Fidelity cannot be larger than 1
    if(newFidelity >= 1)
        disp('Error: fidelity cannot be larger than 1');
        return;
    end
    
    %Condition to check if interval between two last step has a maximum
    if(newFidelity < currentFidelity.Intensity + newStep*grad0 ||...
            newFidelity <= oldFidelity)
        interval = [oldStep newStep];
        break; %Go to second phase
    end 
    
    %Strong Wolfe condition of gradient decrease
    %If ok, all conditions are respected and end of line search
    if(abs(newGrad) <= sigma*grad0)
        bestStepSize = newStep;
        bestFidelity.Intensity = newFidelity;
        bestPulse = newPulse;
        optDir.newGradient = newGradient;
        return;
    end
    
    %Condition of maximum included in interval, ie upper limit has negative
    %gradient
    if(newGrad <= 0)
        interval = [oldStep newStep];
        break; %Go to second phase
    end
    
    %Stop algorithm if newStep = maxStep, limit already reached and all
    %test failed
    if(newStep == maxStep)
        disp('Cannot find adequate interval');
        return;
    end
    
    %The interval is not good: increase upper limit of interval
    %Increase is increasingly fast with iteration
    if(maxStep <= 2*newStep - oldStep) %If close to max
        oldStep = newStep;
        newStep = maxStep;
    else
        %Otherwise use cubic interpolation
        old2Step = oldStep;
        oldStep = newStep;
        newStep = maxCspline([2*oldStep-old2Step min(maxStep,...
            oldStep+expandSpeed*(oldStep - old2Step))],[old2Step oldStep],...
            [oldFidelity newFidelity],[oldGrad newGrad]);
    end
    
    %Memorize fidelity and gradient
    oldFidelity = newFidelity;
    oldGrad = newGrad;    
end

%Fidelity and gradient at interval boundaries
fidelInt = [oldFidelity newFidelity];
gradInt = [oldGrad newGrad];

%Second phase: find step maximizer using previous calculated interval and
%reduce it until step respect Wolfe conditions 
while(true)
    %Update step size
    
    newStep = maxCspline([interval(1)+collapseSpeedLow*(interval(2)-interval(1)) ...
        interval(2)-collapseSpeedUp*(interval(2)-interval(1))],...
        interval, fidelInt, gradInt);
    
    %Make sure that progress is expected in fidelity
    if(abs((interval(1)-newStep)*gradInt(1)) <= epsFid)
        bestStepSize = interval(1);
        bestFidelity.Intensity = fidelInt(1);
        
        bestPulse.Pulse = currentPulse.Pulse + bestStepSize*pulseDirection;
        maxPulse = max(abs(bestPulse.Pulse(:)));
        if(maxPulse > params.maxPower)
             bestPulse.Pulse(:) = bestPulse.Pulse(:)*params.maxPower/maxPulse;
        end
        optDir.newGradient = newGradient;
        return;
    end
    
    %Update pulse, fidelity and gradient
    [newPulse, newFidelity, newGradient] = updatePulseFidelity(newStep);
    
    %Calculate gradient for 'reduced' fidelity
    newGrad = newGradient(:).'*pulseDirection(:);
    
    %Check if step is not too large
    if(newFidelity < currentFidelity.Intensity + rho*newStep*grad0 ||...
            newFidelity <= fidelInt(1))
        interval(2) = newStep;
        fidelInt(2) = newFidelity; 
        gradInt(2) = newGrad;
        continue;
    end
    
    %Strong Wolfe condition of gradient decrease
    %If ok, all conditions are respected and end of line search
    if(abs(newGrad) <= sigma*grad0)
        bestStepSize = newStep;
        bestFidelity.Intensity = newFidelity;
        bestPulse = newPulse;
        optDir.newGradient = newGradient;
        return;
    end
    
    %Step size is not good so reduce interval using following boundaries  
    if(newGrad <= 0)
        interval(2) = newStep;
        fidelInt(2) = newFidelity; 
        gradInt(2) = newGrad;
        continue;
    else
        interval(1) = newStep;
        fidelInt(1) = newFidelity; 
        gradInt(1) = newGrad;
        continue;
    end       
end

end

