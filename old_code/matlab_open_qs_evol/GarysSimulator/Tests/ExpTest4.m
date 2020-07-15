clear all;

%% Loading H0

NV = Spin(3/4);
Carbon = Spin(1/2);
Nitrogen = Spin(1);

%% Start Expm vs Trotter expansion Test
maxField = 6;
maxCarbon = 4;
dt = 4e-9;
Len = 100;
power = 1e7;
disc = 30;
scale = [0 2*disc]+1;
CPUtime = zeros(maxCarbon+1,maxField);
CPUtimeTrot = zeros(maxCarbon+1,maxField);
TrotterError = zeros(maxCarbon+1,maxField);

for ctCarbon = 0:maxCarbon
    ctCarbon
    LauncherHamiltonian(3/4,ctCarbon,[0 5],150*1e-4);
    load('HnatData.mat') %Hamiltonian from simulator, variable HnatData
    
    controls{1} = HnatData.expandOperator('NV',NV.Sx);
    controls{2} = HnatData.expandOperator('Nitrogen',Nitrogen.Sx)*1e-3;
    for ctC2 = 1:ctCarbon
        controls{ctC2+2} = HnatData.expandOperator(strcat('Carbon',num2str(ctC2)),Carbon.Sx)*1e-3;
    end
    
    for ctField = 1:length(controls)
        controls2{ctField} = controls{ctField}*dt*power;
    end
      
    for ctField = 1:length(controls)
        ctField
        
        propagator = PropagatorObj();
        propagator.TimeLength = Len;
        propagator.TimeStep = dt;
        propagator.Hnat = HnatData.RotatingMatrix;
        propagator.ControlMats = controls(1:ctField);
        
        propagator2 = PropagatorObj2();
        propagator2.PulseLength = Len;
        propagator2.AmpDisc = (-disc:1:disc)/disc;
        propagator2.Hnat = HnatData.RotatingMatrix*dt;
        propagator2.ControlMats = controls2(1:ctField);
        propagator2.precalculation(5);
        
        currentPulse2 = PulseObj2([ctField Len]);
        currentPulse2.Index(:) = currentPulse2.makeRandomIndex(0.1,scale);
        currentPulse2.Pulse(:) = propagator2.AmpDisc(currentPulse2.Index(:));
        propagator2.ControlFields = currentPulse2;
        
        currentPulse = PulseObj([ctField Len]);
        currentPulse.Pulse = currentPulse2.Pulse*power;
        propagator.ControlFields = currentPulse;

        tic
        Uexp = propagator.unitary2();
        CPUtime(ctCarbon+1,ctField) = toc;
        clear propagator
        
        tic
        Utrot = propagator2.unitary();
        CPUtimeTrot(ctCarbon+1,ctField) = toc;
        clear propagator2
        
        opEnd = 1;
        for ct=1:length(Uexp)
            opEnd = Uexp{ct}*opEnd;
        end
        opEndTrot = 1;
        for ct=1:size(Utrot,3)
            opEndTrot = Utrot(:,:,ct)*opEndTrot;
        end
        err = (1-abs(trace(opEnd'*opEndTrot/length(opEndTrot))))^2;
        TrotterError(ctCarbon+1,ctField) = err;
    end
end
imagesc((CPUtime./CPUtimeTrot).')
