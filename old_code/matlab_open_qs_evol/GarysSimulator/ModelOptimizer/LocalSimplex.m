classdef LocalSimplex < handle
    %SIMPLEX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SpaceDim %Number of coordinate per point, parameter space dimension
        Points %Column = coordinate, row = different points (spacedim+1)
        Costs %Cost of different points      
    end
    
    methods
        function obj = LocalSimplex(dimension)
            obj.SpaceDim = dimension;
        end
        
        %Make a random initial simplex
        function InitSimplex(obj, initPoint, edgeLength)
            obj.Points = repmat(initPoint,1,obj.SpaceDim+1);
            directions = orth(rand(obj.SpaceDim));
            obj.Points(:,2:end) = obj.Points(:,2:end)+...
                (rand(1)+1)/2*edgeLength*directions; %[L/2,L]
            obj.Costs = zeros(obj.SpaceDim+1,1);
        end
        
        %Sort points according to cost
        function Sort(obj)
            [obj.Costs, index] = sort(obj.Costs);
            obj.Points = obj.Points(:,index);
        end    
    end
    
end

