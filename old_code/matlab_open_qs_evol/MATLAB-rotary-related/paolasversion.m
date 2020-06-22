%===========
close all
clear all

%% data sets to be loaded
name{1} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-11-201924Uno50cycles.mat';
name{2} ='1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-17-202727Uno75cycles.mat';
name{3} ='1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-17-163452Uno87cyclesincomplete.mat';%alternative
name{4} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-14-174935Uno117cycles.mat';
name{5} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-15-034143Uno-164cycles.mat';

%% Load and process the data
for aux = 1:5

load(name{aux});

%% PARAMETERS
ncycles(aux)=Scan.Variable_values{8}.value;
tcycle(aux)=Scan.Variable_values{6}.value*2;
IntT(aux)= ncycles(aux)*tcycle(aux);
ave(aux)=Scan.Averages;
rep(aux)=Scan.Repetitions;

%% SIGNAL
% load signal
% the errors need to be /srqt(Repetitions) bc at the time of the exps we
% were not yet using the std dev of the mean yet
Sig = Scan.ExperimentData{1}{1};
experrorSIG = Scan.ExperimentDataError{1}{1}/sqrt(rep(aux));
% load microwave frequency w [Hz]
w = Scan.vary_begin:((Scan.vary_end-Scan.vary_begin)/(Scan.vary_points-1)):Scan.vary_end;

%% MOVING AVG
lagpts = 16; %needs to be an even number, >=2

Sigma = tsmovavg(Sig, 's', lagpts);
Sigma=Sigma(lagpts:1:end);
% errors get moved too for consistency
errSIG = tsmovavg(experrorSIG, 's', lagpts);
errSIG = errSIG(lagpts:1:end);
w = w(floor(lagpts/2):1:end-floor(lagpts/2));


%% dS/dw
diffe = tsmovavg(diff(Sigma),'s',lagpts)/(w(2)-w(1));
diffe=diffe(lagpts-1:end);
w0 = w(floor(lagpts/2):1:end-floor(lagpts/2));
errSIG = errSIG(floor(lagpts/2):1:end-floor(lagpts/2)); %OR Stadevi(1:1:end-1);

%% Delta(dw)
% I could even put mean(errSIG) and it does not make any difference.
ARR=(errSIG)./abs(diffe);
% Minimum dw: approximate range by trial/error
[deltaw(aux),array_poschosenw] = min(ARR(40+10*aux:70+10*aux));

%% PLOTS AND REGIONS
plot(w0(40+10*aux:70+10*aux),diffe(40+10*aux:70+10*aux)*(w(2)-w(1)),'.')
hold on
plot(w0,diffe*(w(2)-w(1)))
plot(w(40+10*aux+7:77+10*aux),Sigma(40+10*aux+7:77+10*aux)-mean(Sigma),'k.')
plot(w0(40+10*aux:70+10*aux),errSIG(40+10*aux:70+10*aux)-mean(errSIG),'.m')
plot(w,Sigma-mean(Sigma),'k')
plot(w0,errSIG-mean(errSIG),'m')
plot(w0(array_poschosenw+40+10*aux),errSIG(array_poschosenw+40+10*aux)-mean(errSIG),'or')
plot(w0(array_poschosenw+40+10*aux),errSIG(array_poschosenw+40+10*aux)-mean(errSIG),'*r')
pause
hold off
end

%% Ideal Sensitivity
tempi = 1e-6:0.1e-6:10e-6;
%deltawid=exp((tempi/7e-6).^2/2)./tempi/sqrt(ave(1)*rep(1))/0.0086/2/pi;
deltawid=exp((tempi/5e-6).^2)./tempi/sqrt(ave(1)*rep(1))/0.0086/2/pi;

%% PLOT
figure(10)
hold on
plot(IntT,deltaw,'-*k')
plot(tempi,deltawid,'r')



%==========