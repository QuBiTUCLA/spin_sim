close all;
clear
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Params that need to be changed

%simulation loaded below uses a pi flop for a certain rabi frequency. Enter the pi flop
%time below
pi_time = 25e-9;

% repetitions per point
repsperpt = 1e5; 
% averages
avgs = 50;

%State to be compared: Psi = amp*B + (1-amp)*D
%amp = 0.5; % equal superposition of B and D
%amp = 1; % init at B - worse case picture
amp = 0; %init at D

%using max number of photons in one shot = 3
%using threshold for pi-detection 0 photon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('cla_2d_arr.mat');

display('loaded data ...');

%get number of detection lengths and waiting times
%a(1) is nb of detection lengths
%a(2) is nb of waiting times
a = size(res);

results_simulationPIAMP = zeros(a(1),a(2));
results_simulationSTDPIAMP = zeros(a(1),a(2));
results_simulationUSUAL = zeros(a(1),a(2));
results_simulationSTDUSUAL = zeros(a(1),a(2));
results_simulationPIAMPSTATS = zeros(a(1),a(2));
% results_simulation(# det times, # waiting times)


%quantum_efficiency = 5e-4; % this is a very rough assumption by assuming 200ns detection time and 20 kcps photon count rate for the bright state and comparing that to the sum of NbarB (=7) given by the simulation.
quantum_efficiency = 1; 

carrier_decay = 77e6 * quantum_efficiency;
cross_decay = 1.5e6 * quantum_efficiency;

max_no_of_photons = 20;


%for kappa = 1:10:a(2) %loop over waiting times
for kappa = a(2)    
    kappa/a(2)
    
tdetection_arr = [];

for kk=1:1:a(1)
tdetection_arr = [tdetection_arr res{kk,kappa}.detection_duration];
end


for auxi = 1:1:a(1)
    
    %res{auxi,kappa}.waiting_time %just to see where we are
    
    tdetection = tdetection_arr(auxi);
    
    %%%%%%%%%%%%%%%%%%%%%% Bright state
    [min_difference, B_array_position_end_det1] = min(abs(res{auxi,kappa}.init0.T - tdetection)); %usual detection
    [min_difference, B_array_position_begin_det2] = min(abs(res{auxi,kappa}.init0.T - (tdetection + res{auxi,kappa}.waiting_time + pi_time))); %2nd det for pi-det
    
    %%%%%%%%%%%%%%%%%%%%%% Dark state
    [min_difference, D_array_position_end_det1] = min(abs(res{auxi,kappa}.init1.T - tdetection));
    [min_difference, D_array_position_begin_det2] = min(abs(res{auxi,kappa}.init1.T - (tdetection + res{auxi,kappa}.waiting_time + pi_time)));
 
   
% Now I want the avg nb of photons for each tstep
%these below are the avg value of the distributions, assumed to be
%Poissonian bc delta t small
NbarB = zeros(1,B_array_position_end_det1);
NbarD = zeros(1,D_array_position_end_det1);
%Bright
for q = 1:1:B_array_position_end_det1
  tstep = res{auxi,kappa}.init0.T(q+1) - res{auxi,kappa}.init0.T(q);
  NbarB(q) = (res{auxi,kappa}.init0.pop_zero_e(q) + res{auxi,kappa}.init0.pop_mone_e(q) + res{auxi,kappa}.init0.pop_one_e(q))*(carrier_decay)*tstep + cross_decay*(res{auxi,kappa}.init0.pop_mone_e(q) + res{auxi,kappa}.init0.pop_one_e(q))*tstep + cross_decay*(2*res{auxi,kappa}.init0.pop_zero_e(q))*tstep;
end
%Dark
for q = 1:1:D_array_position_end_det1
  tstep = res{auxi,kappa}.init1.T(q+1) - res{auxi,kappa}.init1.T(q);
  NbarD(q) = (res{auxi,kappa}.init1.pop_zero_e(q) + res{auxi,kappa}.init1.pop_mone_e(q) + res{auxi,kappa}.init1.pop_one_e(q))*(carrier_decay)*tstep + cross_decay*(res{auxi,kappa}.init1.pop_mone_e(q) + res{auxi,kappa}.init1.pop_one_e(q))*tstep + cross_decay*(2*res{auxi,kappa}.init1.pop_zero_e(q))*tstep;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cheating to test code
%NbarB = zeros(1, length(NbarB));
%NbarB(1) = 1;
%NbarD = zeros(1, length(NbarD));
%NbarD(1) = 0.0001;
%%%%%%%%%%%%%%%%%%%%%%%%%%


poi_distrB = []; % list of poisson distr
poi_distrD = []; 
AvgB = 0;
AvgD = 0;
%poi_distr(nb of photons k, time interval)
karray = 0:1:max_no_of_photons;
for k = karray % k is the nb of photons in the histo; assuming never to go beyond 5
    %Cla
      poi_distrB(k+1) = sum(((NbarB).^k).*exp(-NbarB)./factorial(k))/length(NbarB);
      poi_distrD(k+1) = sum(((NbarD).^k).*exp(-NbarD)./factorial(k))/length(NbarD);

    %Boe
	%  poi_distrB(k+1) = (sum(NbarB).^k).*exp(-sum(NbarB))./factorial(k);
    %  poi_distrD(k+1) = (sum(NbarD).^k).*exp(-sum(NbarD))./factorial(k);

      AvgB = AvgB + k*poi_distrB(k+1);
      AvgD = AvgD + k*poi_distrD(k+1);
end





%%%%%%%%%%%%%%%%%%%%

% Now I want the avg nb of photons for each tstep FOR PI DETECTION
%these below are the avg value of the distributions, assumed to be
%Poissonian bc delta t small
NbarBpi = zeros(1,length(res{auxi,kappa}.init0.T)-B_array_position_begin_det2);
NbarDpi = zeros(1,length(res{auxi,kappa}.init1.T)-D_array_position_begin_det2);

%Bright
for q = B_array_position_begin_det2:1:length(res{auxi,kappa}.init0.T)
  tstep = res{auxi,kappa}.init0.T(q) - res{auxi,kappa}.init0.T(q-1);
  NbarBpi(q) = (res{auxi,kappa}.init0.pop_zero_e(q) + res{auxi,kappa}.init0.pop_mone_e(q) + res{auxi,kappa}.init0.pop_one_e(q))*(carrier_decay)*tstep + cross_decay*(res{auxi,kappa}.init0.pop_mone_e(q) + res{auxi,kappa}.init0.pop_one_e(q))*tstep + cross_decay*(2*res{auxi,kappa}.init0.pop_zero_e(q))*tstep;
end
%Dark
for q = D_array_position_begin_det2:1:length(res{auxi,kappa}.init1.T)
  tstep = res{auxi,kappa}.init1.T(q) - res{auxi,kappa}.init1.T(q-1);
  NbarDpi(q) = (res{auxi,kappa}.init1.pop_zero_e(q) + res{auxi,kappa}.init1.pop_mone_e(q) + res{auxi,kappa}.init1.pop_one_e(q))*(carrier_decay)*tstep + cross_decay*(res{auxi,kappa}.init1.pop_mone_e(q) + res{auxi,kappa}.init1.pop_one_e(q))*tstep + cross_decay*(2*res{auxi,kappa}.init1.pop_zero_e(q))*tstep;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cheating to test code
%NbarBpi = zeros(1, length(NbarBpi));
%NbarBpi(1) = 0.0001;
%NbarDpi = zeros(1, length(NbarDpi));
%NbarDpi(1) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%






poi_distrBpi = []; % list of poisson distr
poi_distrDpi = []; 
AvgBpi = 0;
AvgDpi = 0;
%poi_distr(nb of photons k, time interval)
karray = 0:1:max_no_of_photons;
for k = karray % k is the nb of photons in the histo; assuming never to go beyond 5
    %Cla
      poi_distrBpi(k+1) = sum(((NbarBpi).^k).*exp(-NbarBpi)./factorial(k))/length(NbarBpi);
      poi_distrDpi(k+1) = sum(((NbarDpi).^k).*exp(-NbarDpi)./factorial(k))/length(NbarDpi);
    %Boe
     % poi_distrBpi(k+1) = (sum(NbarBpi).^k).*exp(-sum(NbarBpi))./factorial(k); % this is the distribution for the init0
     % poi_distrDpi(k+1) = (sum(NbarDpi).^k).*exp(-sum(NbarDpi))./factorial(k); % this is the distribution for the init1

      AvgBpi = AvgBpi + k*poi_distrBpi(k+1);
      AvgDpi = AvgDpi + k*poi_distrDpi(k+1);
end









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BETTER COMPARISON OUR METHOD + PI DET METHOD
% via Monte Carlo

%CDF for 1st detection, taking into account chosen initial state,
%parametrized by defined amplitude amp in B state
% First_cdf(nb of detected photons+1), so indices go from 1 to 4
First_cdf = amp*cumsum(poi_distrB) + (1-amp)*cumsum(poi_distrD);

%************************** OUR METHOD ************************************
Exp_array = zeros(avgs,repsperpt);
%Throw the dice
Dice = rand(avgs,repsperpt);
for p=1:1:avgs
   for l=1:1:repsperpt
      %if Dice(p,l) <= First_cdf(1)  
          %Do nothing bc array already has zeros
      no_of_photons = min(find(Dice(p,l) <= First_cdf))-1;
      if ~isempty(no_of_photons)
        Exp_array(p, l) = no_of_photons;
      else
        Exp_array(p, l) = max_no_of_photons;
      end          
   end
end

MEDIA = zeros(1,avgs); %mean nb of photons for each avg
for p=1:1:avgs
   MEDIA(p) = mean(Exp_array(p,:));
end

% figure(3);
% hist(MEDIA)
% ylabel('occurrences')
% xlabel('mean number of photons per average')

%Get values of averaged nb of photons for B,D from refs (ie, poisson distrs)
%scaling MEDIA vector taking into account MEDIAB and MEDIAD
MEDIA_NORM = (MEDIA - AvgD)./(AvgB - AvgD);
figure(4000);
hist(MEDIA_NORM, [-1:0.1:2])
ylabel('occurrences')
xlabel(['mean simulated amp =' num2str(mean(MEDIA_NORM))  ' for an input amp = ' num2str(amp) ', stdev amp = ' num2str(std(MEDIA_NORM))])
title(['Current detection: Reps per point =' num2str(repsperpt) ', Avgs =' num2str(avgs) ', Detection time =' num2str(tdetection)])

%****************************** PI DETECTION ******************************

%Using for the first detection event the above Exp_array
%threshold is 0 photon
% in First_state_array, 0 codes D state, 1 codes B state
First_state_array = zeros(avgs,repsperpt);
for p=1:1:avgs
for l=1:1:repsperpt
    if Exp_array(p,l) ~= 0
        First_state_array(p,l) = 1;       
    end
end
end

%Make a pi flip and throw the dice again

%CDF for 2st detection, taking into account chosen initial state,
%parametrized by defined amplitude amp in B state
% First_cdf(nb of detected photons+1), so indices go from 1 to 4
%already all the switches for pi pulse have been made

%original
%Second_cdf = (1-amp)*cumsum(poi_distrBpi) + (amp)*cumsum(poi_distrDpi);
%test, since we already do the flip in the state
Second_cdf = (amp)*cumsum(poi_distrBpi) + (1-amp)*cumsum(poi_distrDpi);


Exp_array_pi = zeros(avgs,repsperpt);
%Throw the dice
Dice = rand(avgs,repsperpt);
for p=1:1:avgs
   for l=1:1:repsperpt
      %if Dice(p,l) <= First_cdf(1)  
          %Do nothing bc array already has zeros
      no_of_photons = min(find(Dice(p,l) <= Second_cdf))-1;
      if ~isempty(no_of_photons)
        Exp_array_pi(p, l) = no_of_photons;
      else
        Exp_array_pi(p, l) = max_no_of_photons;
      end          
   end
end

%Using for the second detection event the above Exp_array_pi
%threshold is 0 photon
% in Second_state_array, 0 codes D state, 1 codes B state
Second_state_array = zeros(avgs,repsperpt);
for p=1:1:avgs
for l=1:1:repsperpt
    if Exp_array_pi(p,l) ~= 0
        Second_state_array(p,l) = 1;       
    end
end
end

%Compare state arrays
% where the Result array is one is an event to use; 
% where the Result array is zero is an event to discard.
Result = xor(First_state_array,Second_state_array);
PIAMP = [];
PIAMP_stat = []; % corresponds to the amount of statistics left
for p=1:1:avgs
	indices = [];
	indices = find(Result(p,:) > 0);
	clean_array = [];
	for l=1:1:length(indices)
		clean_array = [clean_array First_state_array(p, indices(l))];
	end
	PIAMP(p) = mean(clean_array); % this takes the mean of the repetitions per avg
	PIAMP_stat(p) = length(clean_array)/length(First_state_array(p, :)); %how much is kept
end

figure(5000);
hist(PIAMP, [-1:0.1:2])
ylabel('occurrences')
xlabel(['mean simulated amp =' num2str(mean(PIAMP))  ' for an input amp = ' num2str(amp) ', stdev amp = ' num2str(std(PIAMP))])
title(['Pi detection: Reps per point =' num2str(repsperpt) ', Avgs =' num2str(avgs) ', Detection time =' num2str(tdetection) ', waitt = ' num2str(res{auxi,kappa}.waiting_time)])


% addition by McBoe
% 
% figure(6);
% subplot(2,1,1);
% bar([0:3], log10(repsperpt * [poi_distrB.' poi_distrD.']))
% subplot(2,1,2);
% bar([0:3], log10(repsperpt * [poi_distrBpi.' poi_distrDpi.']))
% 
% figure(7);
% full_evol_B = zeros(7, 2*length(tpops{1}));
% full_evol_D = zeros(7, 2*length(tpops{1}));
% for k = 1:7
%     full_evol_B(k, :) = [tpops{k} tpopspi{k}];
%     full_evol_D(k, :) = [tpops2{k} tpopspi2{k}];
% end;
% subplot(2,1,1);
% plot(full_evol_B.');
% subplot(2,1,2);
% plot(full_evol_D.');

% save stuff
results_simulationPIAMP(auxi,kappa) = mean(PIAMP);
results_simulationSTDPIAMP(auxi,kappa) = std(PIAMP);
results_simulationPIAMPSTATS(auxi,kappa) = mean(PIAMP_stat);
results_simulationUSUAL(auxi,kappa) = mean(MEDIA_NORM);
results_simulationSTDUSUAL(auxi,kappa) = std(MEDIA_NORM);

end %end loop over detection times

figure(kappa);
plot(tdetection_arr, results_simulationUSUAL(:,kappa), ...
    tdetection_arr, results_simulationPIAMP(:,kappa), ...
    tdetection_arr, results_simulationSTDUSUAL(:,kappa), ...
    tdetection_arr, results_simulationSTDPIAMP(:,kappa))
legend('Mean Curr', 'Mean Pi', 'Std Curr', 'Std Pi');
xlabel('Detection Time');
title(['Waiting time = ' num2str(res{auxi,kappa}.waiting_time)])

end %end loop over waiting times



