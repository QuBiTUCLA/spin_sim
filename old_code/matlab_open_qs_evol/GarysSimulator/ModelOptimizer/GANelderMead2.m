function [mapping,mappingCost] = GANelderMead2(simplex,...
    mapping,mappingCost,modelCost,propagator,HnatModel,maxIter)
%Nelder Mead simplified for short local investigation in GA

%Define Nelder-Mead coefficients
refCoef = 1; %Reflection coefficient
shrCoef = 0.8; %Small shrink every generation
mappingSize = min(3,size(simplex.Points,1));

%Shrink step
meanPts = repmat(mean(simplex.Points,2),1,size(simplex.Points,2));
simplex.Points = meanPts - shrCoef*(meanPts-simplex.Points);
simplex.SimplexEvaluation(modelCost,propagator,HnatModel);
simplex.Sort();

mapping = [mapping simplex.Points(1:mappingSize,:)];
mappingCost = [mappingCost;simplex.Costs];

%Start NM reflections
iterNb = 0;
while(simplex.Costs(end)-simplex.Costs(1) > 1e-4 && iterNb<maxIter)
    iterNb = iterNb + 1;
    reflexNumber = 1;

    while (reflexNumber <= simplex.SpaceDim)
        %Evaluate reflected points
        centerPt = mean(simplex.Points(:,1:end-reflexNumber),2);
        centerPts = repmat(centerPt,1,reflexNumber);
        
        refPoints = centerPts + refCoef*(centerPts - ...
                     simplex.Points(:,end-reflexNumber+1:end));

        refCosts = zeros(reflexNumber,1);
        for ctPts = 1:reflexNumber
            point = refPoints(:,ctPts);
            refCosts(ctPts) = modelCost.MakeCost(propagator,HnatModel,point);
        end

        mapping = [mapping refPoints(1:mappingSize,:)];
        mappingCost = [mappingCost;refCosts]; %#ok<*AGROW>

        %Improvement ?
        minRefCost = min(refCosts);
        if(minRefCost < simplex.Costs(end-reflexNumber+1))

            %Update points
            simplex.Points(:,end-reflexNumber+1:end) = refPoints;
            simplex.Costs(end-reflexNumber+1:end) = refCosts;
            simplex.Sort();

            break; %get out of reflection loop
        else
            reflexNumber = reflexNumber + 1;
        end
    end
    if(reflexNumber > simplex.SpaceDim)
        break;
    end
end