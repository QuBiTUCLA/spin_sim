%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Function to get best step size optimized for Golden Method%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [bestTempPulse, bestTempFidelity, bestStepSize] = ...
    GoldStepOpt2(params, propagator, currentPulse, currentFidelity, optDir)

%Temporary pulse and fidelity for section method
tempPulse = currentPulse;
tempFidelity = currentFidelity;

bestTempPulse = currentPulse;
bestTempFidelity = currentFidelity;
bestStepSize = 0;
noImprovCount = 0;

initStepSize = params.initStepSize;

goldenNumber = (sqrt(5)-1)/2;

%Last change on direction
pulseDirection = zeros(size(optDir.newDirection));
for ctField = 1:size(optDir.newDirection,1)
    maxDir = max(abs(optDir.newDirection(ctField,:)));
    if(maxDir == 0 || isnan(maxDir) || isinf(maxDir)) 
        pulseDirection(ctField,:) = 0;
    else
        pulseDirection(ctField,:) = optDir.newDirection(ctField,:)/maxDir/25;
    end
end

%Function to update pulse and fidelity
    function updatePulseFidelity(stepSize)
        %New Pulse
        tempPulse.Index = currentPulse.Index + ...
            round(stepSize*pulseDirection*params.ampDisc); %Floor or round?
        tempPulse.Index(:) = max(1,min(2*params.ampDisc+1,tempPulse.Index(:)));
        tempPulse.Pulse(:) = propagator.AmpDisc(tempPulse.Index(:));        
        
        %Propagate (exact and approximate)
        propagator.ControlFields = tempPulse; %propagator use index
        [opEnd,U] = propagator.fullPropagation(params.optType);
        
        %Check fidelity
        [tempFidelity.Intensity tempFidelity.BaseIntensity] = tempFidelity.makeFidelity(opEnd,params,tempPulse.Pulse);

        %Update bestPulse
        if(tempFidelity.Intensity < min(fidelities(:,2)))
            bestTempPulse = tempPulse;
            bestTempFidelity = tempFidelity;
            bestStepSize = stepSize;
            
            propagator.Unitary = U;
            propagator.opEnd = opEnd;
            
            noImprovCount = 0;
        else
            noImprovCount = noImprovCount + 1;
        end
        
        fidelities(end+1,:) = [stepSize tempFidelity.Intensity];
    end

%Initialization: 0-InitStepsize boundary + first two interior boundaries
%Boundary = 0, ie currentPulse, currentFidelity
fidelities(1,:) = [0 currentFidelity.Intensity];
lowBound = 0;

%Boundary = InitStepsize
updatePulseFidelity(initStepSize); 
highBound = initStepSize;

noImprovCount = 0;

%Boundary = first lower inside boundary
updatePulseFidelity((1-goldenNumber)*initStepSize);
lowIndex = 3;

%Boundary = first upper inside boundary
updatePulseFidelity(goldenNumber*initStepSize);
highIndex = 4;

%Start main loop
stepSizeLoopFlag = 1;
while(stepSizeLoopFlag)

    %Update boundaries and recalculate fidelities
    if(fidelities(lowIndex,2) > fidelities(highIndex,2))
        highBound = fidelities(highIndex,1);
        highIndex = lowIndex;
        updatePulseFidelity(lowBound + (1-goldenNumber)*(highBound - lowBound));
        lowIndex = size(fidelities,1);
    else
        lowBound = fidelities(lowIndex,1);
        lowIndex = highIndex;
        updatePulseFidelity(lowBound + goldenNumber*(highBound - lowBound));
        highIndex = size(fidelities,1);
    end
    
    %Stop loop if bestfidelity do not change after a few loop, prevent long
    %search for little or no improvment (see update function for noImprovCount)
    stepSizeLoopFlag = noImprovCount < 4;
    
    %Make sure that interval is not too small
    stepSizeLoopFlag = stepSizeLoopFlag &&...
        ((highBound - lowBound) > params.minStepInterval*initStepSize);
end

% %Interpolation to find the best step size
% if(size(fidelities,1)>3 && bestStepSize ~= 0 && bestStepSize ~= initStepSize)
%     Xinterp = linspace(fidelities(1,1),fidelities(2,1),size(fidelities,1)*10);
%     [~,maxStepIndex] = min(interp1(fidelities(:,1),fidelities(:,2),Xinterp,'spline')); 
%     updatePulseFidelity(Xinterp(maxStepIndex));
% end

end

