function [results,parameters] = LauncherModelOpt(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Launcher: Hamiltonian model optimization %%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Gary Wolfowicz, May 2011 %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Simplex Coding Genetic Algorithm Method %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
% profile on

% clear all;
disp('Initializing...');

%Load Settings
% if(nargin == 1)
%     eval(varargin{1}); %Change in settings
% end
ModelOptSettings;
if(nargin == 1)
    eval(varargin{1}); %Change in settings
end

%Create Spins and Interactions, get
%'HnatModel','carbonIntTensors','averageCoupling'
SpinCreationOpt;

%Use the defined interactions to create a full, exact Hamiltonian
%Carbon interaction are zero at initialization
HnatModel.createFullHamiltonian();

%Create Rotating frame and define resonance in model
HnatModel.createRotatingFrame(params.Resonance);

%Create propagator for real sample if simulated by NV for testing
realSystem = PropagatorModelObj();
realSystem.rhoIn = params.rhoin;
realSystem.Hnat = HnatData.RotatingMatrix;
realSystem.TimeStep = params.timeStep;
realSystem.TimeLength = params.timeLength;
realSystem.ControlMats = params.ControlMats;

%Create propagator of model
propagator = PropagatorModelObj();
% rhoin = kron(kron(diag([0 1]),eye(3)/3),[1 0.5;0.5 1]/2);
% propagator.rhoIn = rhoin/trace(rhoin);
propagator.rhoIn = params.rhoin;
propagator.Hnat = HnatModel.RotatingMatrix;
propagator.TimeStep = params.timeStep;
propagator.TimeLength = params.timeLength;
propagator.ControlMats = params.ControlMats;

%Projection operator
realSystem.ProjOp = HnatData.expandOperator('NV',diag([0 1]));
propagator.ProjOp = realSystem.ProjOp;

%Simulate data
pulseData = cell(params.numberOfTrial,1);
sampleData = zeros(params.numberOfTrial,params.timeLength);
% pulse = PulseObj([length(propagator.ControlMats) params.timeLength]);

for ctTrial=1:params.numberOfTrial
%     pulseData{ctTrial} = pulse.makeRandomPulse(params.randPower/ctTrial,...
%         params.initPulseScale);
    
    pulseData{ctTrial} = params.initPulseScale*pulseFiltering(...
        2*rand(params.timeLength,length(propagator.ControlMats))-1,...
        params.timeStep,'low pass',params.LPFilter).';
    
    realSystem.ControlFields = pulseData{ctTrial};
    sampleData(ctTrial,:) = realSystem.fullPropagation2();
end

%Noise in sample (sampleData = [0,1])
sampleData = sampleData + params.noiseRatio*(rand(size(sampleData))-1/2);

%Precalculation of pulse exponentials
propagator.ControlFields = pulseData;
propagator.precalculation(params.trotter);

%Create parameter space to tensor transformation matrices
%Define number of parameter to optimize
space2Tensor = cell(length(carbonIntTensors),1); %#ok<*NODEF>
paramNb = 0; %Number of parameter
for ctTensor=1:length(carbonIntTensors)
    if(~strcmp(carbonIntTensors{ctTensor}.Object1.Name,'NV'))
        space2Tensor{ctTensor} = paramNb + [1 2 3; %Index
                                            2 4 5;
                                            3 5 6];
        paramNb = paramNb + 6;
    else %Secular approximation SzIx, SzIy, SzIz only
%         if(params.symmetryFlag)
            space2Tensor{ctTensor} = [0         0         0;
                                      0         0         0;
                                      paramNb+1 paramNb+1 paramNb+2];
            paramNb = paramNb + 2;
%         else
%             space2Tensor{ctTensor} = [0         0         0;
%                                       0         0         0;
%                                       paramNb+1 paramNb+2 paramNb+3];
%             paramNb = paramNb + 3;
%         end
    end
end

if(params.initPopulation < paramNb)
    error('Initial population should be above number of parameter (%d)',...
        paramNb);
end

%Create Model Cost object (similar to (in)fidelity)
modelCost = ModelCost(sampleData,carbonIntTensors,space2Tensor,averageCoupling);

%% Optimization
tic

%Initialization
disp('Creating initial population...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SCGA defined parameters
initialPopulationAccept = 0.7;
maxGenerations = params.maxGenerations;
spaceSize = params.valueRange;
NMiteration = 10; %increase every generation
selecMaxCoef = 1.2; %Selective pressure [1 2]
crossoverProba = 0.8; %[0 1] high crossover = many parents = longer
mutationProba = 0.1; %[0 1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SCGA initial population
population = Population(params.initPopulation,paramNb);

%Create initial points for chromosomes
initPoints = zeros(paramNb,population.Size);
initPoints(:,1) = (2*rand(paramNb,1)-1)*spaceSize;
nbInit = 1;
trialCount = 0;
normalDensity = 1;
while(nbInit < population.Size)
    trialCount = trialCount + 1;
    trialPoint = (2*rand(paramNb,1)-1)*spaceSize;
    dist = zeros(paramNb,nbInit);
    for ct=1:nbInit
        dist(:,ct) = abs(trialPoint(:)-initPoints(:,ct))/(2*spaceSize);  
    end
    %Minimum allowed distance between each simplex
    proba = min(max(dist,[],1))*normalDensity;
    if(proba > 1-initialPopulationAccept)
        nbInit = nbInit + 1;
        initPoints(:,nbInit) = trialPoint;
    end   
    
    %If population too large for spaceSize
    if(trialCount > 100*population.Size)
        normalDensity = normalDensity*2;
        trialCount = 0;
    end
end

%Create chromosomes
edgeLength = spaceSize/(params.initPopulation^(1/paramNb))/normalDensity;
for ctPop=1:population.Size;
    population.Chromosomes{ctPop}.InitSimplex(initPoints(:,ctPop), edgeLength);   
end

%Fast Nelder Mead to the population
for ctPop=1:population.Size
    [mapping,mappingCost] = GANelderMead(population.Chromosomes{ctPop},...
        [],[],modelCost,propagator,HnatModel,NMiteration);
end
population.Sort();

%% GA loop
disp('Starting genetic algorithm...');
disp(sprintf('Initial average cost: %0.5g', population.AverageCost()));
disp(sprintf('\nPopulation size: %d', population.Size));
selecMinCoef = 2 - selecMaxCoef;

genNb = 0;
while(genNb < maxGenerations && ...
        population.Chromosomes{1}.Costs(1) > params.goodCost && ...
        (population.AverageCost()-population.Chromosomes{1}.Costs(1))>1e-2)
        
    genNb = genNb + 1;
    disp(sprintf('Generation: %d', genNb));
    
    %Selection (Baker's linear ranking selection, "roulette wheel")
    selection = Population(population.Size,paramNb);
    
    if(population.Size > 2)
        selecWheel = 1/population.Size*(selecMaxCoef - (selecMaxCoef-selecMinCoef)*...
            ((1:population.Size)-1)/(population.Size-1));

        selectedPop = zeros(population.Size,1);
    
        %For good selection, random must be as good
        randVal = rand(population.Size,1);
        for ctPop=1:population.Size
            ctPop2 = 1;
            proba = selecWheel(1);
            while(randVal(ctPop) > proba)
                ctPop2 = ctPop2 + 1;
                proba = proba + selecWheel(ctPop2);
            end

            selectedPop(ctPop) = ctPop2;

            %Copy selected chromosomes (careful as handle)
            selection.Chromosomes{ctPop} = population.Chromosomes{ctPop2}.Copy();
        end
    else
        for ctPop=1:population.Size
            selection.Chromosomes{ctPop} = population.Chromosomes{ctPop}.Copy();
        end
    end
    
    %Crossover and mutation: new generation
    generationSize = 0;
    newGeneration = {};
    
    %Choose parents (at least 2)
    parentIdx = [];
    while(length(parentIdx) < 2)
        parentIdx = find(rand(selection.Size,1) < crossoverProba);
        parentIdx = unique(parentIdx);
    end 
    
    while(length(parentIdx) > 1)
        
        %Choose mating parents
        if(length(parentIdx) ~= 3)
            matingSize = min(length(parentIdx),randi([2 3],1));
            if(length(parentIdx) == 4)
                matingSize = 2;
            end
            
            %Find minimum distance between parents
            dist = zeros(length(parentIdx)-1,1);
            parent1 = mean(selection.Chromosomes{parentIdx(1)}.Points,2);
            for ctPop=2:length(parentIdx)
                parent2 = mean(selection.Chromosomes{parentIdx(ctPop)}.Points,2);
                dist(ctPop) = norm(parent1-parent2);
            end
            [~,minDistIdx] = sort(dist);
            parentsMating = parentIdx(minDistIdx(1:matingSize));
            
        else
            matingSize = 3;
            parentsMating = parentIdx;
        end
        parentIdx = setdiff(parentIdx,parentsMating);

        %Create children
        children = Population(matingSize,paramNb);
        
        %Find average parent
        meanPts = 0;
        for ctPop=1:children.Size
            meanPts = meanPts + selection.Chromosomes{...
                parentsMating(ctPop)}.Points;
        end
        meanPts = meanPts/children.Size;
        
        %Find maximum distance between parents
        dist = zeros(children.Size*(children.Size-1)/2,1);
        k=0;
        for ctPop=1:children.Size
            for ctPop2=ctPop+1:children.Size
                parent1 = mean(selection.Chromosomes{parentsMating(ctPop)}.Points,2);
                parent2 = mean(selection.Chromosomes{parentsMating(ctPop2)}.Points,2);
                k=k+1;
                dist(k) = norm(parent1-parent2);
            end
        end
        maxDist = max(dist);
        
        %Make new children positions
        for ctPop=1:children.Size
            children.Chromosomes{ctPop}.Points = meanPts + maxDist*...
                repmat(2*rand(size(meanPts,1),1)-1,1,size(meanPts,2));
        end
        
        %Mutation
        mutChildren = find(rand(children.Size,1) < mutationProba);
        mutChildrenNb = length(mutChildren);
        
        for ctPop=1:mutChildrenNb
            mutatedChild = children.Chromosomes{mutChildren(ctPop)};
            
            vtx2Ref = randi([1 size(mutatedChild.Points,2)],1);
            meanVtx = mean([mutatedChild.Points(:,1:vtx2Ref-1) ...
                mutatedChild.Points(:,vtx2Ref+1:end)],2);
            
            mutatedChild.Points(:,vtx2Ref) = meanVtx + (0.5 + rand(1))*...
                (meanVtx - mutatedChild.Points(:,vtx2Ref));
        end
        
        %Nelder-Mead short local optimization
        for ctPop=1:children.Size
            [mapping,mappingCost] = GANelderMead(children.Chromosomes{ctPop},...
                mapping,mappingCost,modelCost,propagator,HnatModel,NMiteration);
            
            generationSize = generationSize + 1;
            newGeneration{generationSize} = children.Chromosomes{ctPop}; %#ok<*SAGROW>
        end  
    end
    
    %Speed up convergence (to a limit)
    NMiteration = min(20,NMiteration + 1);
    
    %Replace 1/10 of worst old generation chromosomes
    %1/10 to preserve diversity
    allCosts = zeros(generationSize + population.Size,1);
    for ctGen=1:generationSize
        allCosts(ctGen) = newGeneration{ctGen}.Costs(1);
    end
    for ctPop=1:population.Size
        allCosts(generationSize + ctPop) = ...
            population.Chromosomes{ctPop}.Costs(1); %#ok<*SAGROW>
    end
    
    [~, index] = sort(allCosts);
    index = index(1:population.Size);
    
    genIdx = index(index <= generationSize);
    popIdx = setdiff(1:population.Size,...
        index(index > generationSize) - generationSize);
   
    for ctPop=1:min(length(popIdx),max(1,round(population.Size/10)))
        population.Chromosomes{popIdx(ctPop)} = ...
            newGeneration{genIdx(ctPop)}.Copy();
    end
    population.Sort();
    
    disp(sprintf('Average cost: %0.5g|Best cost: %0.5g', ...
        population.AverageCost(), population.Chromosomes{1}.Costs(1)));
    
    %Reduction of the population every 50/initpop generations, must have at
    %least 2 parents to mate after reduction
    if(rem(genNb,ceil(maxGenerations/params.initPopulation)) == 0 && ...
            population.Size > 2)
        newPopSize = population.Size - 1;
        population.Chromosomes = population.Chromosomes(1:newPopSize);
        population.Size = newPopSize;
        
        disp(sprintf('\nPopulation size: %d', population.Size));
    end       
end

%Final Nelder Mead
[bestPoint,bestCost,mapping,mappingCost] = NelderMead(...
    population.Chromosomes{1}.Points(:,1),mapping,mappingCost,...
    modelCost,propagator,HnatModel,1,params.costError, 400);

disp(sprintf('\nTotal time: %0.5g s', toc));
disp(sprintf('Function evaluations: %d', length(mappingCost)));
disp(sprintf('Final best cost: %0.5g', bestCost));

%% Checking
%Assign best result
for ctTensor=1:length(carbonIntTensors)
    carbonIntTensors{ctTensor}.Tensor = zeros(3);
    nzIdx = find(space2Tensor{ctTensor});
    carbonIntTensors{ctTensor}.Tensor(nzIdx) = averageCoupling(ctTensor)*...
        bestPoint(space2Tensor{ctTensor}(nzIdx)); %#ok<*AGROW>
end

HnatModel.createFullHamiltonian();
HnatModel.createRotatingFrame(params.Resonance);
propagator.Hnat = HnatModel.RotatingMatrix;

%Check trotter
normalCheck = zeros(params.numberOfTrial,params.timeLength);
trotterCheck = zeros(params.numberOfTrial,params.timeLength);
for ctTrial=1:params.numberOfTrial
    propagator.ControlFields = pulseData{ctTrial};
    normalCheck(ctTrial,:) = propagator.fullPropagation2();
    trotterCheck(ctTrial,:) = propagator.fullPropagation(ctTrial);
end
errorFun = normalCheck - trotterCheck;
disp(sprintf('Trotter error: %0.5g', norm(errorFun(:))));

%Radius checking (NV-Carbon 1)
for ct1=1:params.numberOfCarbon
    for ct2=1:params.numberOfCarbon
        %ZZ error
        ZZerrorFun(ct1,ct2) = abs((NVCarbonTensor{ct1}(3,3)-...
            carbonIntTensors{ct2}.Tensor(3,3))/NVCarbonTensor{ct1}(3,3)); %#ok<*USENS>
    
        %ZX/ZY error
        Ktheory = sqrt(NVCarbonTensor{ct1}(3,1)^2 + NVCarbonTensor{ct1}(3,2)^2);
        Kmodel = sqrt(carbonIntTensors{ct2}.Tensor(3,1)^2 + ...
            carbonIntTensors{ct2}.Tensor(3,2)^2);
        ZXYerrorFun(ct1,ct2) = abs((Ktheory-Kmodel)/Ktheory);
    end
end
if(params.numberOfCarbon == 1)
    disp(sprintf('\nCarbon tensor values checking'));
    disp(sprintf('ZZ error: %0.5g %%', ZZerrorFun(1,1)*100));
    disp(sprintf('ZX/ZY error: %0.5g %%', ZXYerrorFun(1,1)*100));
end

%% Plot
for ctTrial=1:params.numberOfTrial
    modelData(ctTrial,:) = propagator.fullPropagation(ctTrial);
end

figure(1)
time = (1:params.timeLength)*params.timeStep*1e6;
plot(time,sampleData,'r',time,modelData,'b');
title('Sample vs Model evolution');
xlabel('Time (us)');
ylabel('Density matrix projection to 0 state');

%Scatter plot of all points (3D map)
figure(2)
[mapping, mapIdx] = unique(mapping.','rows');
mapping = mapping.';
mappingCost = mappingCost(mapIdx);

showVal = log(mappingCost); %Cost normalization to see more clearly 
showVal = showVal - min(showVal) + 1e-5;

if(size(mapping,1) == 2)
scatter(params.NVC*mapping(1,:),params.NVC*mapping(2,:),showVal,showVal,'filled');
text(params.NVC*bestPoint(1),params.NVC*bestPoint(2),...
    '\leftarrow\fontsize{16}{\color{red}Best}');
xlabel('SzIx/Iy (Hz)');
ylabel('SzIz (Hz)');
elseif(size(mapping,1) == 3)
scatter3(params.NVC*mapping(1,:),params.NVC*mapping(2,:),...
    params.NVC*mapping(3,:),showVal,showVal,'filled');
text(params.NVC*bestPoint(1),params.NVC*bestPoint(2),...
    params.NVC*bestPoint(3),'\leftarrow\fontsize{16}{\color{red}Best}');
xlabel('SzIx (Hz)');
ylabel('SzIy (Hz)');
zlabel('SzIz (Hz)');
end

title('NV-Carbon 1 dipole tensor search map');
drawnow

% profile off

%% All results
parameters = {HnatData params};
results = {bestPoint bestCost mapping mappingCost ZZerrorFun ZXYerrorFun};

end







