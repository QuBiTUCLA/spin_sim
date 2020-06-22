classdef Hamiltonian < handle
    
    properties
        Matrix 
        Spins = Spin.empty
        Interactions = Interaction.empty
        Field = Field.empty
        RotatingMatrix
        EnergyLevels
        Resonance
    end
    
    methods
        
        %Constructor
        function H = Hamiltonian()
        end
        
        %Methods for adding spins or fields
        function addSpin(H,S)
            H.Spins(end+1) = S;
        end
            
        function addSpins(H,varargin)
            for k=1:length(varargin)
                H.Spins(end+1) = varargin{k};
            end
        end
        
        function addField(H,F)
            H.Field(end+1) = F;
        end
        
        %Method for returning a spin or field
        function obj = findObject(H,ObjName)
            obj = [];
            
            %Find spin
            for k=1:length(H.Spins),
                if strcmp(H.Spins(k).Name,ObjName),
                    obj = H.Spins(k);
                    break
                end
            end
            
            %Find field
            for k=1:length(H.Field),
                if strcmp(H.Field(k).Name,ObjName),
                    obj = H.Field(k);
                    break
                end
            end
            
            %Find Interaction
            for k=1:length(H.Interactions),
                if strcmp(H.Interactions(k).Name,ObjName),
                    obj = H.Interactions(k);
                    break
                end
            end
            
            if isempty(obj),
                error('No Object found with name %s',ObjName);
            end
        end
        
        %Method for creating and adding iteractions
        function createInteraction(H,Object1,Object2,InteractionType,InteractionTensor)
            obj1 = H.findObject(Object1);
            %Sort out whether this is one or two object interaction
            if ~strcmp(class(Object2),'double'),
                obj2 = H.findObject(Object2);
                myInt = Interaction(obj1,obj2,InteractionType,InteractionTensor);
                myInt.Name = strcat(Object1,Object2,InteractionType);
            else
                myInt = Interaction(obj1);
                myInt.Name = strcat(Object1,InteractionType);
                myInt.setType(InteractionType);
                myInt.setTensor(InteractionTensor);
            end
            %Add it to the Hamiltonian object
            H.Interactions(end+1) = myInt;           
        end
        
        function d = HilbertDim(H) %total dimension spin1*spin2*...
            d = 1;
            for k=1:length(H.Spins),
                d = d*(2*H.Spins(k).SpinAngularMomentum + 1);
            end
        end
        
        function d = spinDimensions(H)
            d = zeros(1,length(H.Spins));
            for k=1:length(H.Spins),
                d(k) = H.Spins(k).dimH;
            end
        end
        
        function pos = findSpinPosition(H,Name)
            for k=1:length(H.Spins),
                if strcmp(H.Spins(k).Name,Name),
                    pos = k;
                    break
                else
                    pos = 0;
                end
            end
        end
        
        %Create the full Hamiltonian matrix in the full Hibert space
        function createFullHamiltonian(H)            
            H.Matrix = zeros(H.HilbertDim());
          
            % now we loop over all the interactions and generate matricies
            % Then we order them correctly according to the Spin Ordering
            for intct=1:length(H.Interactions)
                %Create the interaction matrix from tensor
                H.Interactions(intct).createMatrixForm();
                
                % check what kind of interaction it is
                switch H.Interactions(intct).Type
                    %Single spin interactions
                    case {'Zeeman','ZFS'}
                        iName = H.Interactions(intct).Object1.Name;
                        H.Matrix = H.Matrix + H.expandOperator(iName, H.Interactions(intct).Matrix);
                        
                    %Two spins interactions    
                    case {'Hyperfine','Dipolar'}
                        Name1 = H.Interactions(intct).Object1.Name;
                        Name2 = H.Interactions(intct).Object2.Name;
                        p1 = H.findSpinPosition(Name1);
                        p2 = H.findSpinPosition(Name2);
                        pXOR = setxor([p1,p2],1:length(H.Spins));
                        H.Matrix = H.Matrix + expandHilbertSpace(H.Interactions(intct).Matrix,[p1,p2],pXOR,H.spinDimensions());
                end
            end
        end

        function M = expandOperator(H,Name,Op)
            % get the spin position
            spinPos = findSpinPosition(H,Name);

            %Use expand hilbert space to put it together
            M = expandHilbertSpace(Op,spinPos,setxor(spinPos,1:length(H.Spins)),H.spinDimensions());
            %setxor gives the number in the two vectors that are not in
            %both at the same time: here it is the vector 1:length without spinpos
        end 
        
        function createRotatingFrame(H,varargin)
            NV = H.Spins(1); %Rotating frame for the first spin (normally NV)
            envSize = round(H.HilbertDim/NV.dimH);   
            
            %Make RWA on internal Hamiltonian
            HnatRotating = H.Matrix.*kron(eye(NV.dimH),ones(envSize));
            
            if(isempty(varargin))
                %Find eigenenergies
                [~,H.EnergyLevels] = eig(HnatRotating);
                
                partialState = partialtrace(H.EnergyLevels,NV.dimH);%get ms=1 mixed hyperfine
                H.Resonance = max(partialState(:))/envSize;%

                %Get to the closest real energy level
                [~,index] = min(abs(H.EnergyLevels(:)-H.Resonance));
                H.Resonance = H.EnergyLevels(index);
                
            elseif(length(varargin) == 1)
                H.Resonance = varargin{1};                
            else
                error('Too many input arguments');
            end
            
            %Make rotating frame (only for 1st spin defined, ms=1 and -1 must be resolved)
            RotatingOp = H.Resonance*H.expandOperator(NV.Name,NV.Sz);
            
            %Final Hamiltonian
            H.RotatingMatrix = HnatRotating - RotatingOp;
        end
        
    end
end