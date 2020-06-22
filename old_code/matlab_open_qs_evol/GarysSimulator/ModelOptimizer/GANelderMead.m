function [mapping,mappingCost] = GANelderMead(simplex,...
    mapping,mappingCost,modelCost,propagator,HnatModel,maxIter)
%Nelder Mead simplified for short local investigation in GA

%Define Nelder-Mead coefficients
refCoef = 1; %Reflection coefficient
expCoef = 2; %Expansion coefficient
outConCoef = 1/2; %Outside contraction coefficient
inConCoef = -1/2; %Inside contraction coefficient

%Evaluate,sort and map
simplex.SimplexEvaluation(modelCost,propagator,HnatModel);
simplex.Sort();
mappingSize = min(3,size(simplex.Points,1));
mapping = [mapping simplex.Points(1:mappingSize,:)];
mappingCost = [mappingCost;simplex.Costs];

%Start Nelder Mead
iterNb = 0;
shrinkFlag = false;
while(simplex.Costs(end)-simplex.Costs(1) > 1e-4 && iterNb<maxIter)
    iterNb = iterNb + 1;
    %Reflect
    centerPt = 1/simplex.SpaceDim*sum(simplex.Points(:,1:end-1),2);

    refPt = centerPt + refCoef*(centerPt-simplex.Points(:,end));
    refCost = modelCost.MakeCost(propagator,HnatModel,refPt);

    mapping = [mapping refPt(1:mappingSize,:)];
    mappingCost = [mappingCost;refCost]; %#ok<*AGROW>

    if(refCost >= simplex.Costs(1) && refCost < simplex.Costs(end-1))
        simplex.Points(:,end) = refPt;
        simplex.Costs(end) = refCost;

    %Expand    
    elseif(refCost < simplex.Costs(1))
        expPt = centerPt + expCoef*(centerPt-simplex.Points(:,end));
        expCost = modelCost.MakeCost(propagator,HnatModel,expPt);

        mapping = [mapping expPt(1:mappingSize,:)];
        mappingCost = [mappingCost;expCost];

        if(expPt < refPt)
            simplex.Points(:,end) = expPt;
            simplex.Costs(end) = expCost;
        else
            simplex.Points(:,end) = refPt;
            simplex.Costs(end) = refCost;
        end

    %Contract    
    elseif(refCost >= simplex.Costs(end-1))
        %Outside
        if(refCost < simplex.Costs(end))
            outConPt = centerPt + outConCoef*(centerPt-simplex.Points(:,end));
            outConCost = modelCost.MakeCost(propagator,HnatModel,outConPt);

            mapping = [mapping outConPt(1:mappingSize,:)];
            mappingCost = [mappingCost;outConCost];

            if(outConCost <= refCost)
                simplex.Points(:,end) = outConPt;
                simplex.Costs(end) = outConCost;
            else
                shrinkFlag = true;
            end

        %Inside    
        else
            inConPt = centerPt + inConCoef*(centerPt-simplex.Points(:,end));
            inConCost = modelCost.MakeCost(propagator,HnatModel,inConPt);

            mapping = [mapping inConPt(1:mappingSize,:)];
            mappingCost = [mappingCost;inConCost];

            if(inConCost <= simplex.Costs(end))
                simplex.Points(:,end) = inConPt;
                simplex.Costs(end) = inConCost;
            else
                shrinkFlag = true;
            end
        end
    end

    %Shrink
    if(shrinkFlag)
        iterNb = maxIter;
    end
    simplex.Sort();
    
end