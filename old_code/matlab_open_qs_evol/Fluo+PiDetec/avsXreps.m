close all;
clear
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Params that need to be changed

%total number of points, N (must be >10000)
N = 2000000;

%State to be compared: Psi = amp*B + (1-amp)*D
%amp = 0.5; % equal superposition of B and D
amp = 0; % init at B - worse case picture

%using max number of photons in one shot = 3
%using threshold for pi-detection 0 photon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('cla_2d_arr.mat');

%get number of detection lengths and waiting times
%a(1) is nb of detection lengths
%a(2) is nb of waiting times
a = size(res);

results_simulationUSUAL = [];
results_simulationSTDUSUAL = [];
% results_simulation(# det times, # waiting times)

kappa = 1; %loop over waiting times; waiting time = 0 (not used, was for pi-detection)
auxi = 10; %corresponds to 100ns acquisition time; =20 corresponds to 200ns
tdetection = res{auxi,kappa}.detection_duration;
    
    %%%%%%%%%%%%%%%%%%%%%% Bright state
    [min_difference, B_array_position_end_det1] = min(abs(res{auxi,kappa}.init0.T - tdetection)); %usual detection
   
    %%%%%%%%%%%%%%%%%%%%%% Dark state
    [min_difference, D_array_position_end_det1] = min(abs(res{auxi,kappa}.init1.T - tdetection));
    
% Now I want the avg nb of photons for each tstep
%these below are the avg value of the distributions, assumed to be
%Poissonian bc delta t small
NbarB = zeros(1,B_array_position_end_det1);
NbarD = zeros(1,D_array_position_end_det1);
%Bright
for q = 1:1:B_array_position_end_det1
  tstep = res{auxi,kappa}.init0.T(q+1) - res{auxi,kappa}.init0.T(q);
  NbarB(q) = (res{auxi,kappa}.init0.pop_zero_e(q) + res{auxi,kappa}.init0.pop_mone_e(q) + res{auxi,kappa}.init0.pop_one_e(q))*(77e6)*tstep + 1.5e6*(res{auxi,kappa}.init0.pop_mone_e(q) + res{auxi,kappa}.init0.pop_one_e(q))*tstep + 1.5e6*(2*res{auxi,kappa}.init0.pop_zero_e(q))*tstep;
end
%Dark
for q = 1:1:D_array_position_end_det1
  tstep = res{auxi,kappa}.init1.T(q+1) - res{auxi,kappa}.init1.T(q);
  NbarD(q) = (res{auxi,kappa}.init1.pop_zero_e(q) + res{auxi,kappa}.init1.pop_mone_e(q) + res{auxi,kappa}.init1.pop_one_e(q))*(77e6)*tstep + 1.5e6*(res{auxi,kappa}.init1.pop_mone_e(q) + res{auxi,kappa}.init1.pop_one_e(q))*tstep + 1.5e6*(2*res{auxi,kappa}.init1.pop_zero_e(q))*tstep;
end

poi_distrB = []; % list of poisson distr
poi_distrD = []; 
AvgB = 0;
AvgD = 0;
%poi_distr(nb of photons k, time interval)
karray = 0:1:3;
for k = karray % k is the nb of photons in the histo; assuming never to go beyond 5
      poi_distrB(k+1) = sum(((NbarB).^k).*exp(-NbarB)./factorial(k))/length(NbarB);
      poi_distrD(k+1) = sum(((NbarD).^k).*exp(-NbarD)./factorial(k))/length(NbarD);
      AvgB = AvgB + k*poi_distrB(k+1);
      AvgD = AvgD + k*poi_distrD(k+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% via Monte Carlo

%CDF for 1st detection, taking into account chosen initial state,
%parametrized by defined amplitude amp in B state
% First_cdf(nb of detected photons+1), so indices go from 1 to 4
First_cdf = amp*cumsum(poi_distrB) + (1-amp)*cumsum(poi_distrD);
kk = 0;
ratioRtoA = [];
for repp = 10000:5000:N
    
    repp/N
    
    kk = kk + 1; 
    ave = floor(N/repp);
    ratioRtoA = [ratioRtoA repp/ave];

%************************** OUR METHOD ************************************
Exp_array = zeros(ave,repp);
%Throw the dice
Dice = rand(ave,repp);
for p=1:1:ave
   for l=1:1:repp
      %if Dice(p,l) <= First_cdf(1)  
          %Do nothing bc array already has zeros
      no_of_photons = min(find(Dice(p,l) <= First_cdf))-1;
      if ~isempty(no_of_photons)
        Exp_array(p, l) = no_of_photons;
      else
        Exp_array(p, l) = 3;
      end          
   end
end

MEDIA = zeros(1,ave); %mean nb of photons for each avg
for p=1:1:ave
   MEDIA(p) = mean(Exp_array(p,:));
end

%Get values of averaged nb of photons for B,D from refs (ie, poisson distrs)
%scaling MEDIA vector taking into account MEDIAB and MEDIAD
MEDIA_NORM = (MEDIA - AvgD)./(AvgB - AvgD);
% figure(4000);
% hist(MEDIA_NORM, [-1:0.1:2])
% ylabel('occurrences')
% xlabel(['mean simulated amp =' num2str(mean(MEDIA_NORM))  ' for an input amp = ' num2str(amp) ', stdev amp = ' num2str(std(MEDIA_NORM))])
% title(['Current detection: Reps per point =' num2str(repsperpt) ', Avgs =' num2str(avgs) ', Detection time =' num2str(tdetection)])

% save stuff
results_simulationUSUAL(kk) = mean(MEDIA_NORM);
results_simulationSTDUSUAL(kk) = std(MEDIA_NORM);

end

figure(1);
plot(ratioRtoA, results_simulationUSUAL, ...
    ratioRtoA, results_simulationSTDUSUAL)
legend('Mean Curr', 'Std Curr')
xlabel('Ratio Reps to Avgs');