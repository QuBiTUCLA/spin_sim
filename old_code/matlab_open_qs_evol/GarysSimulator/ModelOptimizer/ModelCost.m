classdef ModelCost
    properties
        SampleData %Data set from sample = nb trial x time length
        CarbonIntTensors %Cell of handle to interaction tensors
        Space2Tensor %Transformation from coordinate to tensor
        AverageCoupling %Weigth for each tensors
        Covariance
    end
    
    methods
        function obj = ModelCost(sampleData,carbonIntTensors,space2Tensor,averageCoupling)
            obj.SampleData = sampleData;
            obj.CarbonIntTensors = carbonIntTensors;
            obj.Space2Tensor = space2Tensor;
            obj.AverageCoupling = averageCoupling;
        end

        %Make cost, propagator could be in property as is a handle
        function outCost = MakeCost(obj,propagator,HnatModel,point)
            %Transform point to tensor
            for ctTensor=1:length(obj.CarbonIntTensors)
                obj.CarbonIntTensors{ctTensor}.Tensor = zeros(3);
                nzIdx = find(obj.Space2Tensor{ctTensor}); %nonzero index
                obj.CarbonIntTensors{ctTensor}.Tensor(nzIdx) = ...
                    point(obj.Space2Tensor{ctTensor}(nzIdx))*...
                    obj.AverageCoupling(ctTensor); 
            end
            
            %Update Hamiltonian
            HnatModel.createFullHamiltonian();
            HnatModel.createRotatingFrame(HnatModel.Resonance);
            propagator.Hnat = HnatModel.RotatingMatrix;

            %Evaluate cost
            modelData = zeros(size(obj.SampleData));
            for ctTrial=1:size(modelData,1)
                modelData(ctTrial,:) = propagator.fullPropagation(ctTrial);
            end
            errorFun = obj.SampleData - modelData;
            outCost = errorFun(:).'*errorFun(:)/numel(modelData)*1e3;
        end
    end
end


