%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Test of the Magnus expansion %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
clear all;

%Change directory for own computer
addpath(genpath('C:\Users\GaryW\Desktop\Simulator\'));

%Load Data
load('HnatData.mat') %Hamiltonian from simulator, variable HnatData

%Load Settings
SpectrumSettings; %Prepare all settings, get 'params'

%% Starting

%Prepare propagator
iMin = 10;
iStep = 100;
iMax = 510;
fidelity1 = zeros(iMax,1);
fidelity2 = zeros(iMax,1);
fidelity3 = zeros(iMax,1);
for timeLength = iMin:iStep:iMax
    timeLength
    
    propagator = PropagatorSpectrum();
    params.timeLength = timeLength;
    params.timeStep = 1e-9;
    params.freqStep = 1/(params.timeLength*params.timeStep);
    params.freqLength = round(1/(5*params.timeStep*params.freqStep));
    
    propagator.TimeStep = params.timeStep;
    propagator.FreqStep = params.freqStep;
    propagator.FreqMin = params.freqMin;

    HnatData.toggleInteractionFrame(1); %Rotating frame
    propagator.Hnat = HnatData.InteractionHamiltonian;

    propagator.ControlMats = params.ControlMats;
    propagator.rhoIn = params.rhoin;
    propagator.TimeLength = params.timeLength;

    %Precalculate common used cosin values
    propagator.Precalculation(params.freqLength); 

    pulse = PulseGrape([length(propagator.ControlMats) params.freqLength]); 
    pulse.Pulse = pulse.makeRandomPulse(1,1e6);
    propagator.ControlFields = pulse;
    opEnd1 = propagator.fullPropagation(params.optType);

    freqAxis = 2*pi*(0:params.freqLength-1).'*params.freqStep;
    

    opEnd2 = eye(size(propagator.Hnat));
    ctTimeStep = timeLength;
    for ctTime=1:ctTimeStep
        T0 = params.timeLength/ctTimeStep*(ctTime-1)*params.timeStep;
        Tend = params.timeLength/ctTimeStep*ctTime*params.timeStep;
        omega = propagator.Hnat*(Tend-T0) + pulse.Pulse(1,1)*propagator.ControlMats{1}*(Tend-T0);
%         omega = omega + pulse.Pulse(2,1)*propagator.ControlMats{2}*(Tend-T0);
        for i = 2:length(freqAxis)
            omega = omega + pulse.Pulse(1,i)*...
                ((sin(freqAxis(i)*Tend)-sin(freqAxis(i)*T0))/freqAxis(i)*propagator.ControlMats{1}-...
                (-sin(freqAxis(i)*T0)*freqAxis(i)*T0-2*cos(freqAxis(i)*T0)+...
                sin(freqAxis(i)*T0)*Tend*freqAxis(i)-T0*freqAxis(i)*sin(freqAxis(i)*Tend)+...
                2*cos(freqAxis(i)*Tend)+sin(freqAxis(i)*Tend)*Tend*freqAxis(i))/freqAxis(i)^2*...
                commut(propagator.Hnat,propagator.ControlMats{1})/2);
%             omega = omega + pulse.Pulse(2,i)*...
%                 (cos(freqAxis(i)*T0)-cos(freqAxis(i)*Tend))/freqAxis(i)*propagator.ControlMats{2};   
        end
    
        opEnd2 = expm(-1i*2*pi*omega)*opEnd2;
    end
    
    fidelity1(timeLength) = abs(trace(opEnd1'*opEnd2/length(omega)))^2;
    opEnd3=opEnd2;
    
    Tend = params.timeLength*params.timeStep;
    omega = propagator.Hnat*Tend + pulse.Pulse(1,1)*propagator.ControlMats{1}*Tend;
%     omega = omega + pulse.Pulse(2,1)*propagator.ControlMats{2}*Tend;
    for i = 2:length(freqAxis)
    omega = omega + pulse.Pulse(1,i)*(sin(freqAxis(i)*Tend)/freqAxis(i)*propagator.ControlMats{1}+...
        ((2-2*cos(freqAxis(i)*Tend))/freqAxis(i)^2 - sin(freqAxis(i)*Tend)/freqAxis(i)*Tend)*...
        commut(propagator.Hnat,propagator.ControlMats{1})/2);
    
%     omega = omega + pulse.Pulse(2,i)*(1-cos(freqAxis(i)*Tend))/freqAxis(i)*propagator.ControlMats{2};
    end
    opEnd2 = expm(-1i*2*pi*omega);

%     opEnd2 = eye(size(propagator.Hnat));
%     ctTimeStep = 20;
%     for ctTime=1:ctTimeStep
%         T0 = params.timeLength/ctTimeStep*(ctTime-1)*params.timeStep;
%         Tend = params.timeLength/ctTimeStep*ctTime*params.timeStep;
%         omega = propagator.Hnat*(Tend-T0) + pulse.Pulse(1,1)*propagator.ControlMats{1}*(Tend-T0);
%         omega = omega + pulse.Pulse(2,1)*propagator.ControlMats{2}*(Tend-T0);
%         for i = 2:length(freqAxis)
%             omega = omega + pulse.Pulse(1,i)*...
%                 (sin(freqAxis(i)*Tend)-sin(freqAxis(i)*T0))/freqAxis(i)*propagator.ControlMats{1};
%             omega = omega + pulse.Pulse(2,i)*...
%                 (cos(freqAxis(i)*T0)-cos(freqAxis(i)*Tend))/freqAxis(i)*propagator.ControlMats{2};   
%         end
%     
%         opEnd2 = expm(-1i*2*pi*omega)*opEnd2;
%     end


    
    fidelity2(timeLength) = abs(trace(opEnd1'*opEnd2/length(omega)))^2;
    
    fidelity3(timeLength) = abs(trace(opEnd3'*opEnd2/length(omega)))^2;
end
% A=smooth(fidelity1(iMin:iStep:iMax),30);
% plot(iMin:iStep:iMax,A);
% figure(1)
% plot(iMin:iStep:iMax,fidelity1(iMin:iStep:iMax));
% figure(2)
% plot(iMin:iStep:iMax,fidelity2(iMin:iStep:iMax));
figure(3)
plot(iMin:iStep:iMax,fidelity3(iMin:iStep:iMax));