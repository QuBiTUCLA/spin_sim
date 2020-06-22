clear all;

%% Loading propagator

%% Start Expm vs Trotter expansion Test
maxField = 6;
maxtrial = 30;
CPUtimeTrot = zeros(maxField,1);
CPUtimeExp = zeros(maxField,1);
TrotterError = zeros(maxField,1);

LauncherHamiltonian(3/4,1,[0 8],150*1e-4);
load('HnatData.mat') %Hamiltonian from simulator, variable HnatData
NV = Spin(3/4);
Carbon = Spin(1/2);
Nitrogen = Spin(1);
controls{1} = HnatData.expandOperator('NV',NV.Sx);
controls{2} = HnatData.expandOperator('NV',NV.Sy);
controls{3} = HnatData.expandOperator('Nitrogen',Nitrogen.Sx);
controls{4} = HnatData.expandOperator('Nitrogen',Nitrogen.Sy);
controls{5} = HnatData.expandOperator('Carbon1',Carbon.Sx);
controls{6} = HnatData.expandOperator('Carbon1',Carbon.Sy);

for fieldNum = 1:maxField
    fieldNum
    
    clear propagator
    propagator = PropagatorObj();
    propagator.TimeStep = 1e-9;
    propagator.TimeLength = 500;
    propagator.Hnat = HnatData.RotatingMatrix;
    
    propagator.ControlMats = controls(1:fieldNum);
    propagator.precalculation(2^5);
    
    for trial = 1:maxtrial

        currentPulse = PulseObj([fieldNum propagator.TimeLength]);
        currentPulse.Pulse = currentPulse.makeRandomPulse(0.1,1e7);
        propagator.ControlFields = currentPulse;

        tic
        Utrot=propagator.unitary();
        CPUtimeTrot(fieldNum)=CPUtimeTrot(fieldNum)+toc;

        tic
        Uexp=propagator.unitary2();
        CPUtimeExp(fieldNum)=CPUtimeExp(fieldNum)+toc; 
        
        opEnd = 1;
        for ct=1:length(Uexp)
            opEnd = Uexp{ct}*opEnd;
        end
        opEndTrot = 1;
        for ct=1:length(Utrot)
            opEndTrot = Utrot{ct}*opEndTrot;
        end
        err = (1-abs(trace(opEnd'*opEndTrot/length(opEndTrot))))^2;
        TrotterError(fieldNum)= TrotterError(fieldNum) + err;
    end
    
end
CPUtimeTrot = CPUtimeTrot/maxtrial;
CPUtimeExp = CPUtimeExp/maxtrial;
TrotterError = TrotterError/maxtrial;

plot(1:maxField,CPUtimeTrot,'r',1:maxField,CPUtimeExp,'b')
