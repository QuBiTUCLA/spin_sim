clear all;

%% Loading propagator

%% Start Expm vs Trotter expansion Test
maxField = 1;
maxtrial = 5;

CPUtimeTrot = zeros(10,1);
CPUtimeExp = zeros(10,1);

LauncherHamiltonian(3/4,1,[0 8],150*1e-4);
load('HnatData.mat') %Hamiltonian from simulator, variable HnatData
NV = Spin(3/4);
k=0;
for ctTime=1000:1000:10000
    ctTime
    k=k+1;
    for trial = 1:maxtrial
        trial
        clear propagator
        propagator = PropagatorModelObj();
        propagator.TimeStep = 1e-9;
        propagator.TimeLength = ctTime;
        propagator.Hnat = HnatData.RotatingMatrix;
        propagator.rhoIn = HnatData.expandOperator('NV',diag([0 1])); 
        propagator.ControlMats{1} = HnatData.expandOperator('NV',NV.Sx);
        propagator.ProjOp = HnatData.expandOperator('NV',diag([0 1]));

        pulseData{1} = 1e7*pulseFiltering(2*rand(ctTime,1)-1,1e-9,'low pass',1e9).';
        propagator.ControlFields = pulseData;

        propagator.precalculation(5);

        tic
        propagator.fullPropagation(1);
        CPUtimeTrot(k)=CPUtimeTrot(k)+toc;

        propagator.ControlFields = pulseData{1};

        tic
        propagator.fullPropagation2();
        CPUtimeExp(k)=CPUtimeExp(k)+toc; 
    end
end
CPUtimeTrot = CPUtimeTrot/maxtrial;
CPUtimeExp = CPUtimeExp/maxtrial;

plot(1000:1000:10000,CPUtimeTrot,'r+-',1000:1000:10000,CPUtimeExp,'bo-');
