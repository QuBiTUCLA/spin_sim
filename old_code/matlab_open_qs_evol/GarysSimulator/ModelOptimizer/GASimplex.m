classdef GASimplex < handle
    %SIMPLEX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SpaceDim %Number of coordinate per point, parameter space dimension
        Points %Column = coordinate, row = different points (spacedim+1)
        Costs %Cost of different points      
    end
    
    methods
        function obj = GASimplex(dimension)
            obj.SpaceDim = dimension;
            obj.Points = zeros(obj.SpaceDim,obj.SpaceDim+1);
            obj.Costs = zeros(obj.SpaceDim+1,1);
        end
        
        %Make a random initial simplex
        function InitSimplex(obj, initPoint, edgeLength)
            obj.Points = repmat(initPoint,1,obj.SpaceDim+1);
            for ct=1:obj.SpaceDim
                obj.Points(ct,ct+1) = obj.Points(ct,ct) + edgeLength;
            end
        end
        
        %Sort points according to cost
        function Sort(obj)
            [obj.Costs, index] = sort(obj.Costs);
            obj.Points = obj.Points(:,index);
        end
        
        %Evaluate cost of all vertices
        function SimplexEvaluation(obj,modelCost,propagator,HnatModel)
            for ct=1:length(obj.Costs)
                obj.Costs(ct) = modelCost.MakeCost(propagator,HnatModel,...
                    obj.Points(:,ct));
            end
        end
        
        %Copy simplex (must be done this way as simplex = handle)
        function copy = Copy(obj)
            copy = GASimplex(obj.SpaceDim);
            copy.Points = obj.Points;
            copy.Costs = obj.Costs;
        end
        
    end
    
end

