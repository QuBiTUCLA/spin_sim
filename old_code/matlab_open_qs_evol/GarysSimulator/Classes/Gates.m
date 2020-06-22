classdef Gates < handle
    %Common logic gates
    
    properties
    end
    
    methods (Static)
        function obj = Gates()
        end
        
        function U = Rotation(Hnat,spinName,globalAngle,Xangle,Yangle,Zangle)
            %Single qubit transformations (WRONG FOR Sz in 3/4)
            spin = Hnat.findObject(spinName);
            rotX = expm(-1i*Xangle*spin.Sx); 
            rotY = expm(-1i*Yangle*spin.Sy);
            rotZ = expm(-1i*Zangle*spin.Sz);
            %/spin.SpinAngularMomentum*2
            totalOperator = exp(1i*globalAngle)*rotZ*rotY*rotX;
            U = Hnat.expandOperator(spinName,totalOperator);
        end
        
        function U = Rotation2(Hnat,spinName,globalAngle,Xangle,Yangle,Zangle)
            %Rotation without taking -1 into account (spin 1 seen as spin 1/2)
            sp2 = Spin(1/2);
            spin = Hnat.findObject(spinName);
            
            rotX = eye(spin.dimH);
            rotX(1:2,1:2) = expm(-1i*Xangle*sp2.Sx);
            
            rotY = eye(spin.dimH);
            rotY(1:2,1:2) = expm(-1i*Yangle*sp2.Sy);
            
            rotZ = eye(spin.dimH);
            rotZ(1:2,1:2) = expm(-1i*Zangle*sp2.Sz);

            totalOperator = exp(1i*globalAngle)*rotZ*rotY*rotX;
            U = Hnat.expandOperator(spinName,totalOperator);
        end
        
        function U = CNOT(Hnat,controlName,spinName)
           %control = 0/-1, spin=spin, control = 1, spin = |1><0|+|0><1|
           controlSpin = Hnat.findObject(controlName);
           controlSpinPos = Hnat.findSpinPosition(controlName);
           spin = Hnat.findObject(spinName);
           spinPos = Hnat.findSpinPosition(spinName);
           
           %spin are always 1 0 (-1) in order
           Hcurrent = eye(controlSpin.dimH*spin.dimH); %Reduced system here to control/spin
           Hcurrent(1:2,1:2) = [0 1;1 0];
           
           %Expand reduced system to all Hilbert space
           pXOR = setxor([controlSpinPos,spinPos],1:length(Hnat.Spins));
           U = expandHilbertSpace(Hcurrent,[controlSpinPos,spinPos],pXOR,Hnat.spinDimensions());
        end
        
        function U = SWAP(Hnat,spin1Name,spin2Name)
            %Swap gate: 00->00, 01->10, 10->01, 11->11
            spin1 = Hnat.findObject(spin1Name);
            spin1Pos = Hnat.findSpinPosition(spin1Name);
            spin2 = Hnat.findObject(spin2Name);
            spin2Pos = Hnat.findSpinPosition(spin2Name);
            if(spin1.dimH~=2||spin2.dimH~=2)
                error('SWAP gate only for spin of dimension 2');
            end
            
            %Gate reduced to two bits system
            Hcurrent = eye(spin1.dimH*spin2.dimH);
            Hcurrent(2:3,2:3) = [0 1;1 0];
            
            %Expand reduced system to all Hilbert space
            pXOR = setxor([spin1Pos,spin2Pos],1:length(Hnat.Spins));
            U = expandHilbertSpace(Hcurrent,[spin1Pos,spin2Pos],pXOR,Hnat.spinDimensions());
        end
    end
    
end

