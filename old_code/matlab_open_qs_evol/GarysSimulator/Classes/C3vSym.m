classdef C3vSym
    %Rotation tensors of the C3v symmetry
    %Position transformation
    
    properties
        Sym
    end
    
    methods
        function obj = C3vSym()
            %Rotation 0°
            obj.Sym{1} = [1 0 0;
                          0 1 0;
                          0 0 1];
            %Rotation 120°         
            obj.Sym{2} = [-1/2       sqrt(3)/2 0;
                          -sqrt(3)/2 -1/2      0;
                          0          0         1];   
            %Rotation 240°        
            obj.Sym{3} = [-1/2      -sqrt(3)/2 0;
                          sqrt(3)/2 -1/2       0;
                          0         0          1];   
            %Axial symmetry    
            obj.Sym{4} = [1  0 0;
                          0 -1 0;
                          0  0 1];
            %Axial symmetry + 120°           
            obj.Sym{5} = [-1/2      sqrt(3)/2 0;
                          sqrt(3)/2 1/2       0;
                          0         0         1];
            %Axial symmetry + 240°         
            obj.Sym{6} = [-1/2       -sqrt(3)/2 0;
                          -sqrt(3)/2 1/2        0;
                          0          0          1];      
        end
        
        function carbonPos = getCarbonPosition(obj, currentCarbon, carbonSym)
            %Give carbon position from the hyperfine table, the number are
            %the column in the table
            carbonPos = [currentCarbon(4) currentCarbon(5) currentCarbon(6)];
            carbonPos = obj.Sym{carbonSym}*carbonPos';
            carbonPos = carbonPos'*1e-10; %in Angstrom
        end
    end
    
end

