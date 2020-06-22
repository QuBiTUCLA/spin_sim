classdef PropagatorObj < handle
    
    properties
        Hnat
        ControlFields
        ControlMats
        TimeStep
        TimeLength
        rhoIn
        
        %For exponential calculations
        PCTrotter
        PCExpDiagVec %Vector form of above for faster calculation
        PCExpMats %Base change matrices (with H0)
    end
    
    methods
        
        %Constructor
        function obj = PropagatorObj()
        end
        
        function precalculation(obj,trotter)
            %Precalculation of matrix exponentials
            obj.PCTrotter = trotter;
            obj.PCExpMats = cell(length(obj.ControlMats) + 1);
            obj.PCExpDiagVec = zeros(length(obj.ControlMats),length(obj.Hnat));
            PCExpDiag = cell(length(obj.ControlMats));
            
            %Exp of internal Hamiltonian
            sym = length(obj.ControlMats) + 1;
            H0Exp =  expm(-2*1i*pi*obj.TimeStep*obj.Hnat/sym/2^trotter); %/sym for 2nd order Trotter
            
            %Diagonalization
            [Pnew PCExpDiag{1}] = eig(expm(-2*1i*pi*obj.ControlMats{1}/2^trotter));
            
            obj.PCExpMats{1} = H0Exp*Pnew;
            Pold = Pnew;
            
            if(length(obj.ControlMats)>1)
                for ctField=2:length(obj.ControlMats)
                    [Pnew PCExpDiag{ctField}] = eig(expm(-2*1i*pi*obj.ControlMats{ctField}/2^trotter));
                    obj.PCExpMats{ctField} = Pold\H0Exp*Pnew;
                    Pold = Pnew;
                end  
            end
            
            obj.PCExpMats{length(obj.ControlMats) + 1} = Pnew\H0Exp;
            
            for ctField=1:length(PCExpDiag)
                %Diag into vector for fast access in unitary
                obj.PCExpDiagVec(ctField,:) =  PCExpDiag{ctField}(1:length(PCExpDiag{1})+1:end);
            end
           
%             %Cleaning up and sparsing matrices
%             for ctField=1:length(obj.PCExpMats)
%                 obj.PCExpMats{ctField}(abs(obj.PCExpMats{ctField}) < max(max(abs(obj.PCExpMats{ctField})))*1e-5) = 0;
%                 obj.PCExpMats{ctField} = sparse(obj.PCExpMats{ctField});
%             end
            
        end
        
        function U = unitary(obj)
            %New expm calculation using trotter product and some
            %precalculated matrices, used for gradient as U must be
            %memorized
            U = cell(obj.TimeLength,1);
            
            %Some calculation tricks (pre-allocation, class access seems slow)
            pulse = obj.TimeStep*obj.ControlFields.Pulse;
            maxField = size(pulse,1);
            Hlen = length(obj.Hnat);
            diagElem = 1:Hlen+1:Hlen^2;
            diagVec = obj.PCExpDiagVec;
%             diagMat = speye(Hlen);
            diagMat = eye(Hlen);
            expMats = obj.PCExpMats;
            trotter = 2^obj.PCTrotter;
            
            %Propagators
            for ct = 1:obj.TimeLength
                U{ct} = expMats{1};
                
                for ctField = 1:maxField                  
                    diagMat(diagElem) = diagVec(ctField,:).^pulse(ctField,ct);  %#ok<*SPRIX>
                    U{ct} = U{ct}*diagMat*expMats{ctField+1};
                end         
                U{ct} = U{ct}^trotter;
            end
        end
        
        function opEnd = fullPropagation(obj,optType)
        %Special full propagation tailored for unitary control and trotter
        %Faster and low memory
        %ControlFields must be set
            switch(optType)
                case 'Unitary'    
                    %Some calculation tricks (pre-allocation, class access seems slow)
                    pulse = obj.TimeStep*obj.ControlFields.Pulse;
                    maxField = size(pulse,1);
                    Hlen = length(obj.Hnat);
                    diagElem = 1:Hlen+1:Hlen^2;
                    diagVec = obj.PCExpDiagVec;
                    diagMat = speye(Hlen);
                    expMats = obj.PCExpMats;
                    trotter = 2^obj.PCTrotter;

                    opEnd = 1;
                    %Propagators
                    for ct = 1:obj.TimeLength
                        U = expMats{1};
                        for ctField = 1:maxField                  
                            diagMat(diagElem) = diagVec(ctField,:).^pulse(ctField,ct);
                            U = U*diagMat*expMats{ctField+1};
                        end         
                        opEnd = U^trotter*opEnd;
                    end

                case 'Transfer'
                    U = obj.unitary();
                    opEnd = obj.rhoIn;
                    for ct=1:length(U)
                        opEnd = U{ct}*opEnd*U{ct}';
                    end

                otherwise
                    ERROR('OPTIMIZATION TYPE UNKNOWN');
            end
        end
       
       %Functions for testing
       %Function to make unitary operators for all time
        function [U,HTOT] = unitary2(obj)
            %HTOT should be sparse (Hnat and mat already)
            U = cell(obj.TimeLength,1);
            HTOT = cell(obj.TimeLength,1);
            pulse = obj.ControlFields.Pulse;
            maxField = size(pulse,1);
            controlMat = obj.ControlMats;
            timeStep = obj.TimeStep;
            
            for ct = 1:obj.TimeLength
                HTOT{ct} = obj.Hnat;
                for ctField = 1:maxField
                    HTOT{ct} = HTOT{ct} + pulse(ctField,ct)*controlMat{ctField};
                end 
                U{ct} = expm(-2*1i*pi*timeStep*HTOT{ct}); 
            end
            
        end
        
        function opEnd = fullPropagation2(obj,optType)
        %ControlFields must be set
            [U,~] = obj.unitary2();

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