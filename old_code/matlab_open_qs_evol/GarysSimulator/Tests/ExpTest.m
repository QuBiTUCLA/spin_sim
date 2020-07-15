clear all;

%% Loading propagator

%% Start Expm vs Trotter expansion Test
maxSpin = 5;
maxtrial = 10;
CPUtimeTrot = zeros(maxSpin+1,1);
CPUtimeExp = zeros(maxSpin+1,1);
TrotterError = zeros(maxSpin+1,1);

for spinNum = 0:maxSpin
    spinNum
    
    LauncherHamiltonian(3/4,spinNum,[0 8],150*1e-4);
    load('HnatData.mat') %Hamiltonian from simulator, variable HnatData
    
    clear propagator
    propagator = PropagatorObj();
    propagator.TimeStep = 1e-9;
    propagator.TimeLength = 100;
    propagator.Hnat = HnatData.RotatingMatrix;

    NV = Spin(3/4);    
    propagator.ControlMats{1} = HnatData.expandOperator('NV',NV.Sx);
    propagator.precalculation(2^5);
    
    for trial = 1:maxtrial
        trial

        currentPulse = PulseObj([1 propagator.TimeLength]);
        currentPulse.Pulse = currentPulse.makeRandomPulse(0.1,1e7);
        propagator.ControlFields = currentPulse;

        tic
        Utrot=propagator.unitary();
        CPUtimeTrot(spinNum+1)=CPUtimeTrot(spinNum+1)+toc;

        tic
        Uexp=propagator.unitary2();
        CPUtimeExp(spinNum+1)=CPUtimeExp(spinNum+1)+toc; 
        
        opEnd = 1;
        for ct=1:length(Uexp)
            opEnd = Uexp{ct}*opEnd;
        end
        opEndTrot = 1;
        for ct=1:length(Utrot)
            opEndTrot = Utrot{ct}*opEndTrot;
        end
        err = (1-abs(trace(opEnd'*opEndTrot/length(opEndTrot))))^2;
        TrotterError(spinNum+1)= TrotterError(spinNum+1) + err;
    end
    
end
CPUtimeTrot = CPUtimeTrot/maxtrial;
CPUtimeExp = CPUtimeExp/maxtrial;
TrotterError = TrotterError/maxtrial;

plot(0:maxSpin,CPUtimeTrot,'r',0:maxSpin,CPUtimeExp,'b')
