classdef PropagatorSimObj < handle
    
    properties
        Hnat
        ControlFields
        ControlMats
        TimeStep
        TimeLength
        rhoIn
    end
    
    methods
        
        %Constructor
        function obj = PropagatorSimObj()
        end
        
       %Function to make unitary operators for all time
       %USE ROTFRAME
        function [U,HTOT] = unitary(obj)
            U = cell(obj.TimeLength,1);
            HTOT = cell(obj.TimeLength,1);
            pulse = obj.ControlFields.Pulse;
            rotFrame = obj.ControlFields.RotFrame;
            maxField = size(pulse,1);
            controlMat = obj.ControlMats;
            timeStep = obj.TimeStep;
            HnatMat = obj.Hnat;
            
            for ct = 1:obj.TimeLength
                HTOT{ct} = HnatMat{rotFrame(ct)};
                for ctField = 1:maxField
                    HTOT{ct} = HTOT{ct} + pulse(ctField,ct)*controlMat{ctField};
                end 
                U{ct} = expm(-2*1i*pi*timeStep*HTOT{ct}); 
            end
            
        end
        
        function opEnd = fullPropagation(obj,optType)
        %ControlFields must be set
            [U,~] = obj.unitary();

            switch(optType)
                case 'Transfer'
                    opEnd = obj.rhoIn;
                    for ct=1:length(U)
                        opEnd = U{ct}*opEnd*U{ct}';
                    end

                case 'Unitary'
                    opEnd = 1;
                    for ct=1:length(U)
                        opEnd = U{ct}*opEnd; %opEnd is the full propagator not just at time N
                    end

                otherwise
                    ERROR('OPTIMIZATION TYPE UNKNOWN');
            end      
        end
        
    end
end