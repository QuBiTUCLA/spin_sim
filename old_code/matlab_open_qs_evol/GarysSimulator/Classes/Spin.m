classdef Spin < InteractionObject %child of InteractionObject class
    
    properties
        SpinAngularMomentum
        Sx
        Sy
        Sz
        Svec
    end
    
    methods
        %Constructor
        function obj = Spin(SpinAngularMomentum)
            %3/4 spin is a trick for spin 1 reduced to spin 1/2
            if(SpinAngularMomentum ~= 3/4) 
                obj.SpinAngularMomentum = SpinAngularMomentum;
                obj.Sx = CreateCartesianOperatorMatrix(SpinAngularMomentum,'X');
                obj.Sy = CreateCartesianOperatorMatrix(SpinAngularMomentum,'Y');
                obj.Sz = CreateCartesianOperatorMatrix(SpinAngularMomentum,'Z');
                obj.Svec = {obj.Sx obj.Sy obj.Sz};
            else
                obj.SpinAngularMomentum = 1/2; %Matrix dimension
                
                obj.Sx = CreateCartesianOperatorMatrix(1/2,'X');
                obj.Sy = CreateCartesianOperatorMatrix(1/2,'Y');
                obj.Sz = [1 0;0 0]; %To get ZFS, TO CHANGE, bad trick !!!
                
                obj.Svec = {obj.Sx obj.Sy obj.Sz};
            end
        end
        
        %Size of hilbert space
        function dimension = dimH(Spin)
            dimension = 2*Spin.SpinAngularMomentum +1;
        end
        
    end
end