classdef PropagatorObj2 < handle
    
    properties
        ControlFields
        ControlMats
        Hnat
        
        PulseLength
        
        rhoIn
        
        %For exponential calculations
        PCTrotter
        AmpDisc
        PCExpMats
        
        %Memorize current U propagator cell, avoid recalculating in
        %gradient, same for opEnd, propObj is handle so this is shared and
        %avoid many instances of U which is very memory costly
        Unitary 
        opEnd
    end
    
    methods
        
        %Constructor
        function obj = PropagatorObj2()
            obj.Unitary = 0;
        end
        
        function precalculation(obj,trotter) %Half of discretization could be memorized due to neg = transpose
            %Precalculation of matrix exponentials
            obj.PCTrotter = trotter; %trotter is the power of 2
            obj.PCExpMats = cell(length(obj.ControlMats));
            for ctField=1:length(obj.ControlMats)
                obj.PCExpMats{ctField} = zeros(size(obj.Hnat,1),size(obj.Hnat,2),length(obj.AmpDisc));
            end
            
            %Exp of internal Hamiltonian
            sym = length(obj.ControlMats) + 1;
            H0Exp =  expm(-2*1i*pi*obj.Hnat/sym/2^trotter); %/sym for 2nd order Trotter
            
            for ctAmp = 1:length(obj.AmpDisc)
                obj.PCExpMats{1}(:,:,ctAmp) = H0Exp*expm(-2*1i*pi*...
                    obj.AmpDisc(ctAmp)*obj.ControlMats{1}/2^trotter)*H0Exp;
                
                for ctField=2:length(obj.ControlMats)
                    obj.PCExpMats{ctField}(:,:,ctAmp) = expm(-2*1i*pi*...
                        obj.AmpDisc(ctAmp)*obj.ControlMats{ctField}/2^trotter)*H0Exp;
                end  
            end
            %NEGATIVE AMPLITUDE CAN BE CALCULATED FROM HERMITIAN, CHANGE
            %CODE
            
        end

        function U = unitary(obj)
            %New expm calculation using trotter product and some
            %precalculated matrices, used for gradient as U must be
            %memorized
            
            %Some calculation tricks (pre-allocation, class access seems slow)
            index = obj.ControlFields.Index;
            maxField = size(index,1);
            expMats = obj.PCExpMats;
            trotter = obj.PCTrotter;
            
            %Propagators
            U = expMats{1}(:,:,index(1,:));
            for ctField = 2:maxField                
                U = mtimesx(U,expMats{ctField}(:,:,index(ctField,:)));
            end
            
            for i=1:trotter
                U = mtimesx(U,U);
            end

%             %If not using mtimesx
%             U = expMats{1}(:,:,index(1,:));
%             trotter = 2^trotter;
%             for ct = 1:obj.PulseLength
%                 for ctField = 2:maxField               
%                     U(:,:,ct) = U(:,:,ct)*expMats{ctField}(:,:,index(ctField,ct));
%                 end         
%                 U(:,:,ct) = U(:,:,ct)^trotter;
%             end
        end
        
        function [opEnd,U] = fullPropagation(obj,optType)
        %ControlFields must be set
            U = obj.unitary();

            switch(optType)
                case 'Transfer'
                    opEnd = obj.rhoIn;
                    for ct=1:size(U,3)
                        opEnd = U(:,:,ct)*opEnd*U(:,:,ct)';
                    end

                case 'Unitary'
                    opEnd=ndfun('mprod',U);
                    
%                     opEnd = U(:,:,1);
%                     for ct=2:length(U)
%                         opEnd = U(:,:,ct)*opEnd;
%                     end
                      
                otherwise
                    ERROR('OPTIMIZATION TYPE UNKNOWN');
            end      
        end
        
         
       %Functions for testing
       %Function to make unitary operators for all time        
        function [opEnd,U] = fullPropagation2(obj,optType)
        %ControlFields must be set
            U = cell(obj.PulseLength,1);
            pulse = obj.ControlFields.Pulse;
            maxField = size(pulse,1);
            controlMat = obj.ControlMats;
            
            for ct = 1:obj.PulseLength
                HTOT = obj.Hnat;
                for ctField = 1:maxField
                    HTOT = HTOT + pulse(ctField,ct)*controlMat{ctField};
                end 
                U{ct} = expm(-2*1i*pi*HTOT); 
            end

            switch(optType)
                case 'Transfer'
                    opEnd = obj.rhoIn;
                    for ct=1:length(U)
                        opEnd = U{ct}*opEnd*U{ct}';
                    end

                case 'Unitary'
                    opEnd = eye(size(U{1}));
                    for ct=1:length(U)
                        opEnd = U{ct}*opEnd;
                    end

                otherwise
                    ERROR('OPTIMIZATION TYPE UNKNOWN');
            end      
        end
        
    end
end