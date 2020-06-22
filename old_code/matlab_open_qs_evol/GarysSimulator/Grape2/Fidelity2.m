classdef Fidelity2 %Must not be handle

    properties
        Intensity
        OldIntensity
        GradientOrder
        BaseIntensity
        opGoal;
    end
    
    methods
        function obj = Fidelity2()
            obj.Intensity = 1;
            obj.OldIntensity = 1;
            obj.GradientOrder = 0;
        end 
        
        function grad = makeGradient(obj,type,propagator)
            %order = approximation order in Hausdorff series
            switch(type)
                case 'Transfer'
                    grad = obj.makeTransferGradient(propagator);
                case 'Unitary'
                    grad = obj.makeUnitaryGradient(propagator);
                otherwise
                    ERROR('OPTIMIZATION TYPE UNKNOWN');
            end
        end
              
        %Make gradient for transfer type of optimization
        function grad = makeTransferGradient(obj,propagator)
           order = obj.GradientOrder;
           pulse = propagator.ControlFields.Pulse;
           maxField = size(pulse,1);
           controlMats = propagator.ControlMats;
           
           %prepare backward and forward density matrix and operator
           rho = propagator.rhoIn;
           
           if(~iscell(propagator.Unitary)) 
               propagator.Unitary = propagator.unitary();
           end
           
           backOp = obj.opGoal';
           for ct=size(pulse,2):-1:1 
               backOp = propagator.Unitary(:,:,ct)'*backOp*propagator.Unitary(:,:,ct);
           end
           
           TrEnd = sum(diag((backOp*propagator.rhoIn)'));
           
           %make gradient for every control field and time
           grad = zeros(size(pulse));
           for ct=1:size(pulse,2)
               Uct = propagator.Unitary(:,:,ct);
               
               %Create Htotal for Hausdorff expansion
               HTOT = propagator.Hnat;
               for ctField=1:maxField
                   HTOT = HTOT + controlMats{ctField}.*pulse(ctField,ct); %sum unmix=sum mix property
               end
               
               for ctField=1:maxField
                   mat = Hausdorff2(HTOT,controlMats{ctField},order);
                   grad(ctField,ct) = trace(backOp*(mat*rho+rho*mat'));               
               end

               rho = Uct*rho*Uct';
               backOp = Uct*backOp*Uct';
           end
           
           grad = -2*real(grad*TrEnd);
        end 
        
        function grad = makeUnitaryGradient(obj,propagator)
           %Called many time so faster to assign it (object call is slow)
           order = obj.GradientOrder;
           controlMats = propagator.ControlMats;
           pulse = propagator.ControlFields.Pulse;
           maxField = size(pulse,1);
           
           %prepare backU, frontU (P and U in GRAPE paper)
           %start backU=U2'...opGoal and frontU=U1
           %end backU=opGoal and frontU=UN...U1
           if(~iscell(propagator.Unitary)) 
               propagator.Unitary = propagator.unitary();
               
               backU = obj.opGoal';
               for ct=size(pulse,2):-1:1
                   backU = backU*propagator.Unitary(:,:,ct);
               end
           else
               backU = obj.opGoal'*propagator.opEnd;
           end
           TrEnd = trace(backU');
           
           frontU = eye(size(backU));
           
           %make gradient for every control field and time
           grad = zeros(size(pulse));
           for ct=1:size(pulse,2)
               Uct = propagator.Unitary(:,:,ct);
               
               %Create Htotal for Hausdorff expansion
               HTOT = propagator.Hnat;
               for ctField=1:maxField
                   HTOT = HTOT + controlMats{ctField}.*pulse(ctField,ct); %sum unmix=sum mix property
               end

               for ctField=1:maxField
                   %trace is called too many time and sum(sum(A.'.*(B)) is faster
                   grad(ctField,ct) = sum(sum((backU.').*...
                       (Hausdorff2(HTOT,controlMats{ctField},order)*...
                       frontU)));
               end
               frontU = Uct*frontU;
               backU = backU*Uct';
           end
           
           grad = -2*real(grad*TrEnd);
        end 
        
        %Measure INfidelity from propagated operator and goal operator
        function [out, out2] = makeFidelity(obj, opEnd, params, pulse)
            %pulse must be unmixed pulse.Pulse
            
            %Transfer/unitary transfomation fidelity
            switch(params.optType)
                case 'Transfer'
                    out = 1-abs(trace(obj.opGoal'*opEnd))^2;
                case 'Unitary'
                    out = 1-abs(trace(obj.opGoal'*opEnd/size(obj.opGoal,1)))^2;
                    %size = trace of identity if end = goal
                otherwise
                    ERROR('OPTIMIZATION TYPE UNKNOWN');
            end
            out2 = out;
            
            if(params.maxTimeStep ~= params.minTimeStep)
                %Time penalty                
                if(params.timeGradPenalty ~= 0)
                    out = out +  params.timeGradPenalty*...
                        sum(pulse(1,:))/params.pulseLength;
                end
                
                %Power penalty
                if(params.powerGradPenalty ~= 0)
                    powerPenalty = 0;
                    for ctField=2:size(pulse,1)
                        powerPenalty = powerPenalty + ...
                            pulse(ctField,:)*pulse(ctField,:).';
                    end
                    out = out + params.powerGradPenalty*powerPenalty/...
                        (params.pulseLength*(size(pulse,1)-1));
                end
            else
                %Power penalty
                if(params.powerGradPenalty ~= 0)
                    out = out + pulse(:).'*pulse(:)/numel(pulse)*...
                        params.powerGradPenalty;
                end
            end
        end  
    end 
end

