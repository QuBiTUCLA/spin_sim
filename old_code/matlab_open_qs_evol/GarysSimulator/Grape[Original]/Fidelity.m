classdef Fidelity %Must not be handle
    properties
        Intensity
        OldIntensity
        GradientOrder
        Size %Control and Time size
        opGoal;
    end
    
    methods
        function obj = Fidelity(sizeFieldTime)
            obj.Size = sizeFieldTime;
            obj.Intensity = 0;
            obj.OldIntensity = 0;
            obj.GradientOrder = 0;
        end 
        
        function [grad, opEnd] = makeGradient(obj,type,propagator)
            %order = approximation order in Hausdorff series
            switch(type)
                case 'Transfer'
                    [grad, opEnd] = obj.makeTransferGradient(propagator);
                case 'Unitary'
                    [grad, opEnd] = obj.makeUnitaryGradient(propagator);
                otherwise
                    ERROR('OPTIMIZATION TYPE UNKNOWN');
            end
        end
              
        %Make gradient for transfer type of optimization
        function [grad, opEnd] = makeTransferGradient(obj,propagator)
           U = propagator.unitary();
           maxField = obj.Size(1); %Called many time so faster to assign it
           order = obj.GradientOrder;
           pulse = propagator.ControlFields.Pulse;
           controlMats = propagator.ControlMats;
           timeStep = propagator.TimeStep;
           
           %u(1) -> U(1) + rhoIn -> rho1 so difference in the index of time
           %u(2) = u(1) + g(1) + rho1 -> rho2
           %u(N) + rhoN-1 -> rhoN

           %prepare backward and forward density matrix and operator
           rho = propagator.rhoIn;
           backOp = obj.opGoal';
           for ct=obj.Size(2):-1:1 
               backOp = U{ct}'*backOp*U{ct};
           end
           TrEnd = trace((backOp*propagator.rhoIn)');
           
           %make gradient for every control field and time
           grad = zeros(obj.Size);
           for ct=1:obj.Size(2)
               %Create Htotal for Hausdorff expansion
               HTOT = propagator.Hnat;
               for ctField=1:maxField
                   HTOT = HTOT + controlMats{ctField}.*pulse(ctField,ct);
               end
               
               for ctField=1:maxField
                   mat = Hausdorff(HTOT,controlMats{ctField},timeStep,order);
                   grad(ctField,ct) = trace(backOp*(mat*rho+rho*mat'));               
               end

               rho = U{ct}*rho*U{ct}';
               backOp = U{ct}*backOp*U{ct}';
           end
           
           grad = 2*real(grad*TrEnd);
           opEnd = rho;
        end 
        
        function [grad, opEnd] = makeUnitaryGradient(obj,propagator)
           U = propagator.unitary();
           
           %Called many time so faster to assign it (object call is slow)
           maxField = obj.Size(1); 
           order = obj.GradientOrder;
           controlMats = propagator.ControlMats;
           Hnat = propagator.Hnat;
           timeStep = propagator.TimeStep;
           pulse = propagator.ControlFields.Pulse;
           
           %prepare backU, frontU (P and U in GRAPE paper)
           %start backU=U2'...opGoal and frontU=U1
           %end backU=opGoal and frontU=UN...U1  
           frontU = eye(size(U{1}));
           
           backU = obj.opGoal';
           for ct=obj.Size(2):-1:1
               backU = backU*U{ct}; %Memorization of back U will make algorithm faster TO DO!!!!!
           end
           TrEnd = trace(backU');

           %make gradient for every control field and time
           grad = zeros(obj.Size);
           for ct=1:obj.Size(2)
               %Create Htotal for Hausdorff expansion
               HTOT = Hnat;
               for ctField=1:maxField
                   HTOT = HTOT + controlMats{ctField}.*pulse(ctField,ct);
               end
                   
               for ctField=1:maxField
                   %trace is called too many time and sum(diag()) is faster
                   grad(ctField,ct) = sum(diag(backU*...
                       Hausdorff(HTOT,controlMats{ctField},timeStep,order)*...
                       frontU));
               end

               frontU = U{ct}*frontU;
               backU = backU*U{ct}';
           end
           
           grad = 2*real(grad*TrEnd);
           opEnd = frontU;
        end 
        
        %Measure fidelity from propagated operator and goal operator
        function out = makeFidelity(obj,opEnd,type)
            switch(type)
                case 'Transfer'
                    out = abs(trace(obj.opGoal'*opEnd))^2;
                case 'Unitary'
                    out = abs(trace(obj.opGoal'*opEnd/size(obj.opGoal,1)))^2;
                    %size = trace of identity if end = goal
                otherwise
                    ERROR('OPTIMIZATION TYPE UNKNOWN');
            end       
        end  
    end 
end

