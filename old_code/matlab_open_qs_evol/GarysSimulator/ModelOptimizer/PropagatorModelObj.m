classdef PropagatorModelObj < handle
    
    properties
        Hnat
        ControlFields
        ControlMats
        TimeStep
        TimeLength
        rhoIn
        
        ProjOp
        
        %precalculation
        PulseExp
        PCTrotter
    end
    
    methods
        
        %Constructor
        function obj = PropagatorModelObj()
        end
        
        function precalculation(obj,trotter)
            obj.PCTrotter = trotter;
            maxTrial = length(obj.ControlFields);
            Hdim = length(obj.Hnat);
            obj.PulseExp = cell(maxTrial);
            
            %calculate exponential of control matrix for all times
            for ctTrial = 1:maxTrial
                obj.PulseExp{ctTrial} = zeros(Hdim,Hdim,obj.TimeLength);
                for ct = 1:obj.TimeLength
                    control = obj.ControlFields{ctTrial}(1,ct)*obj.ControlMats{1};
                    for ctField = 2:length(obj.ControlMats);
                        control = control + ...
                            obj.ControlFields{ctTrial}(ctField,ct)*obj.ControlMats{ctField};
                    end
                    obj.PulseExp{ctTrial}(:,:,ct) = ...
                        expm(-2*1i*pi*obj.TimeStep*control/2^trotter);
                end    
            end
        end 
        
        %Special function for model opt, much faster than naive (prop2)
        %Is memory costly due to precalculation (cost of U could be avoided) 
        function projEnds = fullPropagation(obj,ctTrial)
            %Calcul of unitaries
            HnatExp = expm(-2*1i*pi*obj.TimeStep*obj.Hnat/2^obj.PCTrotter);
            U = mtimesx(HnatExp,obj.PulseExp{ctTrial});

            for i=1:obj.PCTrotter
                U = mtimesx(U,U);
            end
            
            %Propagation and projection
%             proj = conj(obj.ProjOp); %(proj').'
            proj = obj.ProjOp(:).';
            projEnds = zeros(1,obj.TimeLength);

            opEnd = obj.rhoIn;
            for ct=1:obj.TimeLength
                opEnd = U(:,:,ct)*opEnd*U(:,:,ct)';
%                 projEnds(1,ct) = sum(sum(proj.*opEnd));
                %If ProjOp is sym/real, otherwise use above commented formula
                %Also use commented proj
                projEnds(1,ct) = proj*opEnd(:);
            end
            projEnds = abs(projEnds);
        end
        
       %WITHOUT PRECALCULATION, for single pulse
        function U = unitary2(obj)
            Hdim = length(obj.Hnat);
            U = zeros(Hdim,Hdim,obj.TimeLength);
            pulse = obj.ControlFields;
            maxField = size(pulse,1);
            controlMat = obj.ControlMats;
            timeStep = obj.TimeStep;
            
            for ct = 1:obj.TimeLength
                HTOT = obj.Hnat;
                for ctField = 1:maxField
                    HTOT = HTOT + pulse(ctField,ct)*controlMat{ctField};
                end 
                U(:,:,ct) = expm(-2*1i*pi*timeStep*HTOT); 
            end
        end
        
        function projEnds = fullPropagation2(obj)
        %ControlFields must be set
            U = obj.unitary2();
            proj = obj.ProjOp';
            projEnds = zeros(1,obj.TimeLength);

            opEnd = obj.rhoIn;
            for ct=1:obj.TimeLength
                opEnd = U(:,:,ct)*opEnd*U(:,:,ct)';
                projEnds(1,ct) = abs(trace(proj*opEnd));
            end
        end
    end
end