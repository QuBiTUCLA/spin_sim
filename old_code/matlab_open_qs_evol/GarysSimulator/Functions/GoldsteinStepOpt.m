%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Function to get best step size optimized for Goldstein %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [bestTempPulse, bestTempFidelity, bestStepSize] = GoldsteinStepOpt(params, propagator, currentPulse, currentFidelity, optDir)

%Temporary pulse and fidelity for section method
tempPulse = currentPulse;
tempFidelity = currentFidelity;

bestTempPulse = currentPulse;
bestTempFidelity = currentFidelity;
bestStepSize = 0;

%Last change on direction
maxDir = max(abs(optDir.newDirection(:)));
if(maxDir == 0 || isnan(maxDir) || isinf(maxDir)) %Singular matrix warning problem to really solve
    pulseDirection = 0;
else
    pulseDirection = (params.maxPower/10)*(optDir.newDirection/maxDir);
end

%Function to update pulse and fidelity
    function updatePulseFidelity(stepSize)
        %New Pulse
        tempPulse.Pulse = currentPulse.Pulse + stepSize*pulseDirection;

        %Limit power
        maxPulse = max(abs(tempPulse.Pulse(:)));
        if(maxPulse > params.maxPower)
             tempPulse.Pulse(:) = tempPulse.Pulse(:)*params.maxPower/maxPulse;
        end
        
        %Propagate (exact and approximate)
        propagator.ControlFields = tempPulse;
        opEnd = propagator.fullPropagation(params.optType);
        
        %Check fidelity
        tempFidelity.Intensity = tempFidelity.makeFidelity(opEnd,params.optType);        
    
        %Update bestPulse
        if(tempFidelity.Intensity > bestTempFidelity.Intensity)
            bestTempPulse = tempPulse;
            bestTempFidelity = tempFidelity;
            bestStepSize = stepSize;
        end
    end

%Initial parameters
alpha1 = 0;
alpha2 = params.initStepSize;
alpha = alpha2;
rho = 0.2;

mainLoopFlag = true;
while(mainLoopFlag)
    updatePulseFidelity(alpha);
    
    %Lower limit
    if(tempFidelity.Intensity >= currentFidelity.Intensity +...
            rho*alpha*optDir.newGradient(:).'*pulseDirection(:))
        %Upper limit
        if(tempFidelity.Intensity <= currentFidelity.Intensity +...
            (1-rho)*alpha*optDir.newGradient(:).'*pulseDirection(:))
            break;
        else
            alpha1 = alpha;
            alpha = (alpha1 + alpha2)/2;
        end
    else
        alpha2 = alpha;
        alpha = (alpha1 + alpha2)/2; 
    end
    
    %Break if interval too small
    mainLoopFlag = mainLoopFlag & ...
        alpha2-alpha1 >= params.minStepInterval*params.initStepSize;
end

end

