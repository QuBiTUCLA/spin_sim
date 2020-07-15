%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Setup all atoms and field vector %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loading data

%NV center
NV = Spin(NVSpinNumber); %First spin defined will be used as rotating frame
NV.Name = 'NV';    

Nitrogen = Spin(1);
Nitrogen.Name = 'Nitrogen';

%% Add spins (NV and random C)

Hnat.addSpin(NV)
Hnat.addSpin(Nitrogen);

j = 0; %number of carbon added index
symetries = C3vSym();
carbonListPosition = zeros(numberOfCarbon,3);
carbonListData = zeros(numberOfCarbon,size(carbonData,2));
while(j < numberOfCarbon) %Add randomly a carbon from the list until numberOfCarbon reached
    
    randCarbon = randi(size(carbonData,1)); %random integer between 1 and size
    randProba = randi(6);
    symNumber = 0;
    addFlag = true;
    
    currentCarbon = carbonData(randCarbon,:);
    
    for i=1:size(carbonListData,1) %check if already all symmetric sites are taken
        if(currentCarbon(1) == carbonListData(i,1))
            symNumber = symNumber + 1;
        end
    end
    if(symNumber >= currentCarbon(2))
        addFlag = false;
    end
    
    %Normalization of the random choice of carbon due to symmetry
    addFlag = addFlag && (randProba <= currentCarbon(2)); 
    
    %Skip number 1 in the list which is the N14 of NV and not a C
    addFlag = addFlag && (currentCarbon(1)~=1);
    
    %Only chooses carbon that are inside the carbonDistRadius value (in A)
    addFlag = addFlag && (currentCarbon(3) <= carbonDistRadius(2));
    addFlag = addFlag && (currentCarbon(3) >= carbonDistRadius(1));

    if(addFlag)
        j = j + 1;
        carbonList(j) = Spin(1/2);
        carbonList(j).Name = strcat('Carbon',num2str(j));
        Hnat.addSpins(carbonList(j));

        carbonListData(j,:) = currentCarbon(:);
        carbonListPosition(j,:) = symetries.getCarbonPosition(currentCarbon,randProba);
        disp(sprintf('Carbon position (A): %0.5g', carbonListData(j,3)));
    end
end

clear carbonData;
clear currentCarbon;







