classdef Interaction < InteractionObject
    
    properties
        Type
        Object1
        Object2
        Tensor %Tensor of interactions values = "J"
        Matrix %Matrix = S*TENSOR*S
    end
    
    methods
        function myInt = Interaction(varargin) %varargin = optional arguments in case
            %1spin interaction (+field, ...) or 2 spin interaction
            if nargin == 0;
                myInt;
            end
            if nargin >=1,
                myInt.Object1 = varargin{1};
            end
            
            if nargin >= 2,
                myInt.Object2 = varargin{2};
            end
            if nargin >=3,
                myInt.Type = varargin{3};
            end
                
            if nargin >=3,
                myInt.Tensor = varargin{4};
            end
       
        end
        
        function myInt = setTensor(myInt,T)
            myInt.Tensor = T;
        end
        
        function myInt = setType(myInt,type)
            myInt.Type = type;
        end
        
        function createMatrixForm(myInt)
            switch myInt.Type
                case 'ZFS'
                    myInt.Matrix = CreateZFS(myInt);
                case 'Zeeman'
                    myInt.Matrix = CreateZeeman(myInt);
                case 'Hyperfine'
                    myInt.Matrix = CreateHyperfine(myInt);
                otherwise
                    error('No Function for creating matrix of interaction type: %s',myInt.Type);
            end
        end
        
        function M = CreateZFS(myInt)
            M = zeros(myInt.Object1.dimH);
            for ct1 = 1:3
                for ct2 = 1:3
                    M = M + myInt.Object1.Svec{ct1}'*myInt.Tensor(ct1,ct2)*myInt.Object1.Svec{ct2};
                end
            end
        end

        function M = CreateZeeman(myInt)
            M = zeros(myInt.Object1.dimH);
            for ct1 = 1:3
                for ct2 = 1:3
                    M = M + myInt.Object2.FieldVector(ct1)*myInt.Tensor(ct1,ct2)*myInt.Object1.Svec{ct2};
                end
            end
        end
        
        
        function M = CreateHyperfine(myInt)
            M = zeros(myInt.Object1.dimH*myInt.Object2.dimH);
            for ct1 = 1:3
                for ct2 = 1:3
                    M = M + myInt.Tensor(ct1,ct2)*kron(myInt.Object1.Svec{ct1},myInt.Object2.Svec{ct2});
                end
            end
                   
        end
    
    end
end