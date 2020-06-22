function [bestPoint,bestCost,mapping,mappingCost] = NelderMead(currentPoint,...
    mapping,mappingCost,modelCost,propagator,HnatModel,localEdgeLength,costError, maxIter)

%Define Nelder-Mead coefficients
refCoef = 1; %Reflection coefficient
expCoef = 2; %Expansion coefficient
outConCoef = 1/2; %Outside contraction coefficient
inConCoef = -1/2; %Inside contraction coefficient
shrCoef = 1/2; %Shrinkage coefficient

%Create a local simplex
simplex = LocalSimplex(size(currentPoint,1));
simplex.InitSimplex(currentPoint, localEdgeLength);
mappingSize = min(3,size(simplex.Points,1)); %CHANGE THIS

%Evaluate simplex costs
disp('Creating local simplex...');
for ctPts = 1:length(simplex.Costs)
    point = simplex.Points(:,ctPts);
    simplex.Costs(ctPts) = modelCost.MakeCost(propagator,HnatModel,point);
end
simplex.Sort();
oldMeanCost = mean(simplex.Costs);

%Save 3D map
mapping = [mapping simplex.Points(1:mappingSize,:)];
mappingCost = [mappingCost;simplex.Costs];

sizeDisp = (simplex.Costs(end)-simplex.Costs(1))/costError;
disp(sprintf('Simplex size: %0.5g', sizeDisp));

iterNb = 0;
shrinkFlag = false;
while(simplex.Costs(end)-simplex.Costs(1) > costError && iterNb<maxIter)
    iterNb = iterNb + 1;
    %Reflect
    centerPt = 1/simplex.SpaceDim*sum(simplex.Points(:,1:end-1),2);

    refPt = centerPt + refCoef*(centerPt-simplex.Points(:,end));
    refCost = modelCost.MakeCost(propagator,HnatModel,refPt);

    mapping = [mapping refPt(1:mappingSize,:)];
    mappingCost = [mappingCost;refCost];

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

    %Kelley restart (Test against stagnation)
    meanCost = mean(simplex.Costs);

    diffFun = simplex.Costs(2:end) - simplex.Costs(1);
    simpDir = simplex.Points(:,2:end)-repmat(simplex.Points(:,1),1,...
        size(simplex.Points,2)-1);
    grad = (simpDir.')\diffFun;

    %Improve condition not respected -> restart
    if(meanCost - oldMeanCost > - 1e-4*norm(grad)^2)
        shrinkFlag = false; %Kelley's step is already a shrink step

        dist = zeros(length(simplex.Costs)-1,1);
        for ct=2:length(simplex.Costs)
            dist(ct-1) = norm(simplex.Points(:,1)-simplex.Points(:,ct));
        end
        sigma_m = min(dist);
        beta = 0.5*sigma_m*sign(grad);

        %Replace points
        kelleyPts = repmat(simplex.Points(:,1),1,length(simplex.Costs)-1);
        for ct=1:size(kelleyPts,2)
            kelleyPts(ct,ct) = kelleyPts(ct,ct) + beta(ct);
        end

        kelleyCosts = zeros(size(kelleyPts,2),1);
        for ctPts=1:size(kelleyPts,2)
            point = kelleyPts(:,ctPts);
            kelleyCosts(ctPts) = modelCost.MakeCost(propagator,HnatModel,point); 
        end
        %Save 3D map
        mapping = [mapping kelleyPts(1:mappingSize,:)];
        mappingCost = [mappingCost;kelleyCosts]; %#ok<*AGROW>

        simplex.Points(:,2:end) = kelleyPts;
        simplex.Costs(2:end) = kelleyCosts;
    end

    %Shrink
    if(shrinkFlag)
        shrinkFlag = false;

        %Evaluate shrink points
        shrPts = repmat(simplex.Points(:,1),1,simplex.SpaceDim);
        shrPts = shrPts + shrCoef*(simplex.Points(:,2:end)-shrPts);
        shrCosts = zeros(size(shrPts,2),1);

        for ctPts = 1:size(shrPts,2)
            %Transform point to tensor
            point = shrPts(:,ctPts);
            shrCosts(ctPts) = modelCost.MakeCost(propagator,HnatModel,point);
        end
        %Save 3D map
        mapping = [mapping shrPts(1:mappingSize,:)];
        mappingCost = [mappingCost;shrCosts];

        simplex.Points(:,2:end) = shrPts;
        simplex.Costs(2:end) = shrCosts;
    end

    simplex.Sort();
    oldMeanCost = meanCost;

    if((simplex.Costs(end)-simplex.Costs(1))/costError < sizeDisp/5)
        sizeDisp = (simplex.Costs(end)-simplex.Costs(1))/costError;
        disp(sprintf('Simplex size: %0.5g', sizeDisp));
    end
    
    bestPoint = simplex.Points(:,1);
    bestCost = simplex.Costs(1);
end
end