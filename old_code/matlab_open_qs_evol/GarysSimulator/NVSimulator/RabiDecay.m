%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Experiment Rabi decay

clear all;
load('HnatData.mat');

NV = Spin(1/2);

%Prepare propagator
propagator = PropagatorObj();
propagator.Hnat = HnatData.RotatingMatrix;
propagator.ControlMats{1} = NV.Sx;
propagator.ControlMats{2} = NV.Sz;
propagator.rhoIn = diag([0 1]);

%Params
propagator.TimeStep = 1e-9;
propagator.TimeLength = 5000;
bathCorel = 1e7;
NVbath = 3e6;
trials = 5000;

output = propagator.rhoIn(:).';

propagator.precalculation(5);

currentPulse = PulseObj([length(propagator.ControlMats) propagator.TimeLength]);
currentPulse.Pulse(1,:) = 1e7;

results = zeros(propagator.TimeLength+1,1);
results(1)=trials;
for meanCt = 1:trials
    meanCt

    randGauss = normrnd(0,1,propagator.TimeLength,1);
    bathField = zeros(propagator.TimeLength,1);
    bathField(1) = NVbath*normrnd(0,1);
    for ct=2:propagator.TimeLength
        bathField(ct) = bathField(ct-1)*exp(-propagator.TimeStep*bathCorel) +...
            NVbath*sqrt((1-exp(-2*propagator.TimeStep*bathCorel)))*randGauss(ct);
    end
    
    currentPulse.Pulse(2,:) = bathField;
    propagator.ControlFields = currentPulse;
    
    U = propagator.unitary();
    rho = propagator.rhoIn;
    for ct=1:length(U)
        rho = U{ct}*rho*U{ct}';
        results(ct+1) = results(ct+1) + output*rho(:);
    end
end
results = abs(results/trials);

x = (0:propagator.TimeLength)*propagator.TimeStep;
T2 = NVbath^2*bathCorel/1e7^2;
hold on
plot(x,0.5)
plot(x,(1+exp(-x*T2))/2)
plot(x,results,'r');