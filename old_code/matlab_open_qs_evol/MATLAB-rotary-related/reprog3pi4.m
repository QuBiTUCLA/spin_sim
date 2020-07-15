%%%%%THIS IS THE FILE USED FOR THE SENSITIVITY PLOT, 3PI/4
%%%%%% NO OTHER

%3pi/4
clear all

doma = 1; %do moving avg for derivative too
donoma = 1;

IntT = 42e-9*[35 49 75];% 90];% 100];

%adaptative lagpts
frequencyscales = [35 49 75 90 100]/35;
varystep = [10000, 10000,10000, 5000, 10000];
initnumber =35; %35 is the correct number to start from if we start pi-rot at 50cycles with 24
lagpontos = round((initnumber-1)*10000./frequencyscales./varystep  +1);

name{1} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-23-212840Uno3pi435cycles.mat';
name{2} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-18-113453Uno3pi449cycles.mat';
name{3} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-19-145847Uno3pip475cycles.mat';
name{4} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-18-072131Uno3pi490cycles.mat';
name{5} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-18-155453Uno3pi4100cyclesnotverygood.mat';

% gyromagnetic ratio
gyro = 1.760859708e11; % s-1 T-1

realangle = 3*pi/4*1.04768;

for aux = 1:1:3

%%% LOAD
clear Stadevi diffe signal rabi lowref oneref zeroref w array_pos arrp Stadevibestfitsignal deriva pbest
load(name{aux});

% load signal (called 'rabi'), 0 ref, 1 ref, errors;
% the errors need to be /srqt(Repetitions) bc at the time of the exps we
% were not yet using the std dev of the mean yet
rabi = Scan.ExperimentData{1}{1};
lowref = Scan.ExperimentData{1}{2};
oneref = Scan.ExperimentData{1}{3};
zeroref = Scan.ExperimentData{1}{4};
experrorONE = Scan.ExperimentDataError{1}{3}/sqrt(100000); %doing moving avg for repetitions as well; 
experrorZERO = Scan.ExperimentDataError{1}{4}/sqrt(100000);  %doing moving avg for repetitions as well
experrorSIG = Scan.ExperimentDataError{1}{1}/sqrt(100000);  %doing moving avg for repetitions as well

% initial, lower and higher bounds for oscillation frequency fitting (param called p(2))
a(aux) = 2*pi*2*IntT(aux)*sin(realangle/2)/realangle;
b(aux) = 2*pi*1.9*IntT(aux)*sin(realangle/2)/realangle; %1.9
c(aux) = 2*pi*2.1*IntT(aux)*sin(realangle/2)/realangle; %2.1



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Paola's correction
% Reps(aux) =  Scan.Repetitions;
% Avgs(aux) = Scan.Averages;
% for k=1:4
%     for na=1:Scan.Averages
%        sig2{na}{k}=(1/sqrt(100000))*Scan.ExperimentDataErrorEachAvg{1}{na}{k}.^2*Reps(aux)*(Reps(aux)-1)+Scan.ExperimentDataEachAvg{1}{na}{k}.^2;    
%     end
% end
% 
% for k=1:4
%     stdd{k}=0;
%     for na=1:Scan.Averages
%         stdd{k}=stdd{k}+sig2{na}{k};
%     end
%     stdd{k}=sqrt(stdd{k}-Scan.ExperimentData{1}{k}.^2)/sqrt(Reps(aux)*Avgs(aux)*(Reps(aux)*Avgs(aux)-1));
% end
% 
% experrorSIG = stdd{1};
% experrorONE = stdd{3};
% experrorZERO = stdd{4};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load varying variable, microwave w [Hz]
w = Scan.vary_begin:((Scan.vary_end -Scan.vary_begin)/(Scan.vary_points-1)):Scan.vary_end;

% normalized signal
rabinorm = (rabi - oneref)./(zeroref - oneref); %rabi normalized

%%% MOVING AVG
lagpts = lagpontos(aux); %needs to be an even number, >=2

one1ma = tsmovavg(oneref, 's', lagpts);
one1ma = one1ma(lagpts:1:end);
zero1ma = tsmovavg(zeroref, 's', lagpts);
zero1ma = zero1ma(lagpts:1:end);
rabima = tsmovavg(rabi, 's', lagpts);
rabima = rabima(lagpts:1:end);




% errors get moved too for consistency
errONE = tsmovavg(experrorONE, 's', lagpts);
errONE = errONE(lagpts:1:end);
errZERO = tsmovavg(experrorZERO, 's', lagpts);
errZERO = errZERO(lagpts:1:end);
errSIG = tsmovavg(experrorSIG, 's', lagpts);
errSIG = errSIG(lagpts:1:end);

%figure(43)
%plot(experrorSIG,'b')
%hold on
%plot(errSIG,'r')
%hold off

%%% NORMALIZATION OF MOVING AVG SIGNAL AND ERRORS
rabinorm2 = (rabima - one1ma)./(zero1ma - one1ma); 
deltaA = errZERO; 
deltaB = errONE;
deltaC = errSIG;
% deltaDbarNEW below is the stdev of the normalized signal
deltaDbarNEW = sqrt(deltaA.^2.*((one1ma - rabima)./(zero1ma - one1ma).^2).^2 + deltaB.^2.*((1./(one1ma - zero1ma) + (rabima - one1ma)./(zero1ma - one1ma).^2).^2)+ deltaC.^2.*(1./(zero1ma - one1ma)).^2);

%%% RENAME, MOVE VECTOR OF FREQUENCIES W
rabinorm = rabinorm2;
if rem(lagpts,2)
    if lagpts == 1
    %w = w(1:1:end-1);  
    else
    w = w(floor(lagpts/2):1:end-floor(lagpts/2)-1);
    end
else
w = w(floor(lagpts/2):1:end-floor(lagpts/2));   
end
figure(1234)
plot(w,rabinorm,'b')
hold on

%%% FIT
myfun = @(p, w) p(1) * cos(p(2)*(p(3) - w) + p(8)*pi) + p(4) + p(1)*cos(p(2)*(p(6) + (p(3) - w)) + p(8)*pi) + p(1)*cos(p(2)*(p(6) - (p(3) - w)) + p(8)*pi) + 0*p(5) + 0*p(7);
pinit = [0.2,     a(aux),      3.157e9,   mean(rabinorm),   0, 2.17e6,   0,     0];
LB =    [0.01,    b(aux),      3.1565e9,             0.01,   0, 1.9e6,   0      0];
UB =    [0.4,     c(aux),      3.1575e9,              0.7,   0, 2.3e6,   0,     0];
[pbest,delta_p]=easyfit(w, rabinorm, pinit, myfun, LB, UB);


% fitted resonance frequency
res(aux) = pbest(3);
% fitted oscillation frequency
oscfreq(aux) = pbest(2);

HF = 2*pi*pbest(6);
pbest(6)
realangle = 3*pi/4*1.04768;
n0 = mean(zero1ma)*1000*100e-9;
n1 = mean(one1ma)*1000*100e-9;
time =  IntT(aux);
%cvector2(aux) = sqrt(1/(1 + 2*(zzz + umumum)/(zzz - umumum)^2));
%cvectorNEW(aux) = sqrt(1/(0.5*(3 + 4/(umumum-zzz) + cos(realangle) + 8*zzz/(zzz-umumum)^2/sin(realangle/2)^2)));
cvectorNEW(aux) = sqrt(1/(1.5 + (-11*n0 + 5*n1)/(2*(n0-n1)^2) +cos(realangle)/2*(1 - (n1+n0)/(n0-n1)^2) + 8*n0/((n0-n1)^2*sin(realangle/2)^2)))*(1/(3/(abs(1+2*cos(2*HF*time*sin(realangle/2)/realangle)))));




% if aux == 3
%   res75 = pbest(3); %- Scan.Variable_values{2}.value;
%   w75 = w;% - Scan.Variable_values{2}.value;
%   sig75 = rabinorm;
%   fit75 = myfun(pbest,w);
%   save('NEW75cycles3pi4.mat','w75','sig75','fit75','res75')
% end
% 
% if aux == 4
%   res90 = pbest(3);%- Scan.Variable_values{2}.value;
%   w90 = w;% - Scan.Variable_values{2}.value;
%   sig90 = rabinorm;
%   fit90 = myfun(pbest,w);
%   save('NEW90cycles3pi4.mat','w90','sig90','fit90','res90') 
% end
% 
% if aux == 5
%   res100 = pbest(3);%- Scan.Variable_values{2}.value;
%   w100 = w;% - Scan.Variable_values{2}.value;
%   sig100 = rabinorm;
%   fit100 = myfun(pbest,w);
%   save('NEW100cycles3pi4.mat','w100','sig100','fit100','res100') 
% end


wnumber(aux) = length(w);

noise(aux) = mean(deltaDbarNEW);
Arms(aux) = (max(myfun(pbest,w)) - min(myfun(pbest,w)))/2/sqrt(2);

%%% Calculate 
Stadevi = deltaDbarNEW; %STD signal 
Stadevi = Stadevi(2:1:end); %OR Stadevi(1:1:end-1);

diffe = abs(diff(rabinorm)./diff(w)); %ORIGINAL
%diffe = abs(diff(myfun(pbest,w))./diff(w)); %TEST
%diffefun = @(p,w) p(1)*2*p(2)*pi*sin(2*p(2)*pi*(p(3)-w)) + p(5)*2*p(2)*pi*sin(2*p(2)*pi*(p(3) + p(6) - w)) - p(7)*2*pi*p(2)*sin(2*p(2)*pi*(-p(3) + p(6) + w));
%diffe = diffefun(pbest,w);
%diffe = diffe(2:1:end);

%%%%%%%%%% MOV AVG FOR DIFFE TOO
if doma
lagpts2 = 2;
diffema = abs(tsmovavg(diff(rabinorm),'s',lagpts2)./diff(w));
diffema=diffema(lagpts2-1:end);
w2 = w(floor(lagpts2/2):1:end-floor(lagpts2/2));
Stadevima = tsmovavg(Stadevi,'s',lagpts2);
Stadevima = Stadevi(lagpts2-1:1:end);
rabinorma = tsmovavg(rabinorm,'s',lagpts2);
deltaDbarNEWma = tsmovavg(deltaDbarNEW,'s',lagpts2);
end

%%% Find minimum within 'nn' oscilation periods
nn = 1;
ene = nn;
figure(3333)
if donoma
[min_diff,arrayposresw] = min(abs(pbest(3)-w)); %position of resonance
[min_diff,array_posmaxw] = min(abs((pbest(3)+nn/2/pbest(2))-w)); 
[min_diff,array_posminw] = min(abs((pbest(3)-nn/2/pbest(2))-w));  
ARR = (Stadevi./abs(diffe));
plot(abs(diffe),'b')
end

if doma
[min_diffma,arrayposreswma] = min(abs(pbest(3)-w2)); %position of resonance
[min_diffma,array_posmaxwma] = min(abs((pbest(3)+nn/2/pbest(2))-w2)); 
[min_diffma,array_posminwma] = min(abs((pbest(3)-nn/2/pbest(2))-w2));  
ARRma = (Stadevima./abs(diffema));
hold on
plot(abs(diffema),'b')
end


%%%%%%%%%%%%%%%%%%%%%%%%%PLOT
figure(1111)
plot(w(1:1:end-1),ARR,'b')
hold on
%plot(w2,ARRma,'r')
hline1 = line([pbest(3) ;pbest(3)],[0 ;5e7]);
set(hline1,'Color','k');
hline2 = line([pbest(3)+ene/2/pbest(2) ;pbest(3)+ene/2/pbest(2)],[0 ;5e7]);
set(hline2,'Color','k');
set(hline2,'LineStyle','--');
hline3 = line([pbest(3)-ene/2/pbest(2) ;pbest(3)-ene/2/pbest(2)],[0 ;5e7]);
set(hline3,'Color','k');
set(hline3,'LineStyle','--');
%%%%%%%%%%%%%%%%%%%%%%%%
if donoma
if array_posmaxw > length(ARR)
    array_posmaxw = length(ARR);
end
end

if doma
if array_posmaxwma > length(ARRma)
    array_posmaxwma = length(ARRma);
end
end
% the uncertainty in frequency is myminNUM [Hz]
if donoma
[myminNUM(aux),array_poschosenw] = min(ARR(array_posminw:1:array_posmaxw));
end
if doma
[myminNUMma(aux),array_poschosenwma] = min(ARRma(array_posminwma:1:array_posmaxwma));
end

%%% Calculate stdev of the uncertainty in frequency at the chosen w that
%%% minimizes the uncertainty in frequency within 'nn' oscilation periods

%analytical, without MA
%stanNUMNEW = lagpts*sqrt( (  ( (1 - 2*rabinorm(1:1:end-1)).^2./(2*sqrt(rabinorm(1:1:end-1).*(1-rabinorm(1:1:end-1)))).^2  )  .* (1./(diffe) ) ).^2   .*    (deltaDbarNEW(1:1:end-1)).^2  );
stanNUMNEW = 1* (  ( (1 - 2*rabinorm(1:1:end-1))./(2*sqrt(rabinorm(1:1:end-1).*(1-rabinorm(1:1:end-1))))  )  .* (1./(diffe) ) )   .*    (deltaDbarNEW(1:1:end-1));
vectstNUMNEW(aux) = stanNUMNEW(array_poschosenw+array_posminw-1);

%analytical, with MA
%stanNUMNEWma = lagpts*lagpts2*sqrt( (  ( (1 - 2*rabinorma(1:1:end-1))./(2*sqrt(rabinorma(1:1:end-1).*(1-rabinorma(1:1:end-1))))  )  .* (1./(diffema) ) ).^2   .*    (deltaDbarNEWma(1:1:end-1)).^2  );
%vectstNUMNEWma(aux) = stanNUMNEWma(array_poschosenwma+array_posminwma-1);


figure(2459)
plot(w,rabinorm,'b')
hold on
hline1 = line([pbest(3) ;pbest(3)],[0 ;1]);
set(hline1,'Color','k');
hline2 = line([pbest(3)+ene/2/pbest(2) ;pbest(3)+ene/2/pbest(2)],[0 ;1]);
set(hline2,'Color','k');
set(hline2,'LineStyle','--');
hline3 = line([pbest(3)-ene/2/pbest(2) ;pbest(3)-ene/2/pbest(2)],[0 ;1]);
set(hline3,'Color','k');
set(hline3,'LineStyle','--');
hold off

end

figure(1492)

% multiply by sqrt(Total interrogation time * N), divide by gyro, convert
% to nT/sqrt(Hz) by multiplying by (2pi)^(3/2)
if donoma
etas = myminNUM*2*pi/gyro.*sqrt(IntT*40*100000)*1e6*sqrt(2*pi); 
etaserror = vectstNUMNEW/gyro*1e6*2*pi*sqrt(2*pi).*sqrt(IntT*40*100000);
end

if doma
etasma = myminNUMma*2*pi/gyro.*sqrt(IntT*40*100000)*1e6*sqrt(2*pi); 
%etaserrorma = vectstNUMNEWma/gyro*1e9*2*pi*sqrt(2*pi).*sqrt(IntT*40*100000);
end

SNfactor = Arms./noise*pi/sqrt(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if donoma
%plot(IntT/1e-6,etas,'b',IntT/1e-6,etas./SNfactor,'k',IntT/1e-6,etas/(lagpts/2),'y',IntT/1e-6,etas./SNfactor/(lagpts/2),'r')

%errorbar(IntT/1e-6,etas,etaserror,'k*') %ORIGINAL
errorbar(IntT/1e-6,etasma,etaserror,'r*')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on

if doma
%plot(IntT/1e-6,etasma,'r',IntT/1e-6,etasma./SNfactor,'g',IntT/1e-6,etasma/lagpts/lagpts2,'m')
end
%errorbar(IntT/1e-6,etas,etaserror,'b')
%errorbar(IntT/1e-6,etasma,etaserrorma,'r')
hold on

res %resonance positions
res-pbest(3)
noise
SNfactor = Arms./noise*pi/sqrt(6);

realangle = 3*pi/4*1.04768;

%Theoretical sensitivity
Cpi = mean(cvectorNEW);   
Cpiworse =  min(cvectorNEW);
Cpibest = max(cvectorNEW);

tempi = 0:0.1e-6:10e-6;
T2Ram =  3.0582e-6; %1.9096e-6; 

% in nT/sqrt(Hz)
%T2 input is Ramsey's T2
sensp = @(time,C,angle,T2) sqrt(2*pi)/C/gyro./sqrt(time)*angle/2/sin(angle/2)^2*1e6.*exp((time/(T2*(angle)/2/sin(angle/2))).^2); 
plot(tempi/1e-6,sensp(tempi,Cpi,realangle,T2Ram),'b',tempi/1e-6,sensp(tempi,Cpiworse,realangle,T2Ram),'--b',tempi/1e-6,sensp(tempi,Cpibest,realangle,T2Ram),'--b')
hold on
%plot all points individually
pointC1 = sensp(IntT(1),cvectorNEW(1),realangle,T2Ram);
pointC2 = sensp(IntT(2),cvectorNEW(2),realangle,T2Ram);
pointC3 = sensp(IntT(3),cvectorNEW(3),realangle,T2Ram);
%pointC4 = sensp(IntT(4),cvectorNEW(4),realangle,T2Ram);
%%%%%
plot(IntT(1)/1e-6,pointC1,'r*')
plot(IntT(2)/1e-6,pointC2,'r*')
plot(IntT(3)/1e-6,pointC3,'r*')
%plot(IntT(4)/1e-6,pointC4,'r*')

axis([0.1 7 1 50])
ylabel('Sensitivity (\muT/\surd Hz)','FontSize',20)
xlabel('Interrogation time (\mus)','FontSize',20)

