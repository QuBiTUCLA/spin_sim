classdef Population < handle
    
    properties
        Chromosomes
        Size
    end
    
    methods
        function obj=Population(totalPopulation,dimension)
            obj.Size = totalPopulation;
            for ct=1:obj.Size
                obj.Chromosomes{ct} = GASimplex(dimension);
            end
        end
        
        %Evaluate cost of each vertex of each simplex
        function PopulationEvaluation(obj,modelCost,propagator,HnatModel)
            for ct=1:obj.Size
                obj.Chromosomes{ct}.SimplexEvaluation(modelCost,propagator,HnatModel);
            end
        end
        
        %Sort vertices in each chromosomes and then all chromosomes
        function Sort(obj)
            costs = zeros(obj.Size,1);
            for ct=1:obj.Size
                obj.Chromosomes{ct}.Sort();
                costs(ct) = obj.Chromosomes{ct}.Costs(1);
            end
            [~, index] = sort(costs);
            obj.Chromosomes = obj.Chromosomes(index);
        end
        
        function avCost = AverageCost(obj)
            avCost = 0;
            for ct=1:obj.Size
                avCost = avCost + mean(obj.Chromosomes{ct}.Costs);
            end
            avCost = avCost/obj.Size;
        end
    end   
end

