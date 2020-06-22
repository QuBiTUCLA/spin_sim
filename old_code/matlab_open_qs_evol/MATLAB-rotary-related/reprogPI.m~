
%%%%%THIS IS THE FILE USED FOR THE SENSITIVITY PLOT, PI


doma = 0; %do moving avg for derivative too
donoma = 1;

% total sequence times for each data set
IntT = 56e-9*[50 75 87 117 164];

%adaptative lagpts
frequencyscales = [50 75 87 117 164]/50;
varystep = [10000, 10000, 5000, 2500, 1250];
initnumber =24; %looks more or less good at those lag pts for 1st scan at 50 cycles
lagpontos = round((initnumber-1)*10000./frequencyscales./varystep  +1); %this the same percentage of mvg avg per cycle
%lagpontos = [10,12,30,30,86] %by hand, good for etasma
%lagpontos= [26,24,50,50,138]; %also by hand
lagpontos = [1,1,1,1,1];

% data sets to be loaded
name{1} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-11-201924Uno50cycles.mat';
name{2} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-12-103314Uno75cycles.mat';
name{3} ='1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-17-163452Uno87cyclesincomplete.mat';%alternative
name{4} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-14-174935Uno117cycles.mat';
name{5} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-15-034143Uno-164cycles.mat';

% initial, lower and higher bounds for oscillation frequency fitting (param called p(2))
a = 1e-6*[2, 3, 3.5, 4.5, 6.5];
b = 1e-6*[1, 2, 3,   4,   6];
c = 1e-6*[3, 4, 4,   5,   7];

% gyromagnetic ratio
gyro = 1.760859708e11; % s-1 T-1

for aux = 1:1:5

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
mean(experrorZERO)
klklk

zz = mean(zeroref)*1000*100e-9;
umum = mean(oneref)*1000*100e-9;
cvector(aux) = sqrt(1/(1 + 2*(zz + umum)/(zz - umum)^2));

% load varying variable, microwave w [Hz]
w = Scan.vary_begin:((Scan.vary_end -Scan.vary_begin)/(Scan.vary_points-1)):Scan.vary_end;

% normalized signal
rabinorm = (rabi - oneref)./(zeroref - oneref); %rabi normalized

%error before ma
deltaa = experrorZERO;
deltab = experrorONE;
deltac = experrorSIG;
deltaDbarBEFma = sqrt(deltaa.^2.*((oneref - rabi)./(zeroref - oneref).^2).^2 + deltab.^2.*((1./(oneref - zeroref) + (rabi - oneref)./(zeroref - oneref).^2).^2)+ deltac.^2.*(1./(zeroref - oneref)).^2);

%%% MOVING AVG
lagpts = lagpontos(aux);

one1ma = tsmovavg(oneref, 's', lagpts);
one1ma = one1ma(lagpts:1:end);
zero1ma = tsmovavg(zeroref, 's', lagpts);
zero1ma = zero1ma(lagpts:1:end);
rabima = tsmovavg(rabi, 's', lagpts);
rabima = rabima(lagpts:1:end);

realangle = pi*1.04768;
zzz = mean(zero1ma)*1000*100e-9;
umumum = mean(one1ma)*1000*100e-9;
cvector2(aux) = sqrt(1/(1 + 2*(zzz + umumum)/(zzz - umumum)^2));
cvectorNEW(aux) = sqrt(1/(0.5*(3 + 4/(umumum-zzz) + cos(realangle) + 8*zzz/(zzz-umumum)^2/sin(realangle/2)^2)));

% errors get moved too for consistency
errONE = tsmovavg(experrorONE, 's', lagpts);
errONE = errONE(lagpts:1:end);
errZERO = tsmovavg(experrorZERO, 's', lagpts);
errZERO = errZERO(lagpts:1:end);
errSIG = tsmovavg(experrorSIG, 's', lagpts);
errSIG = errSIG(lagpts:1:end);

%%% NORMALIZATION OF MOVING AVG SIGNAL AND ERRORS
rabinorm2 = (rabima - one1ma)./(zero1ma - one1ma); 
deltaA = errZERO; 
deltaB = errONE;
deltaC = errSIG;
% deltaDbarNEW below is the stdev of the normalized signal
deltaDbarNEW = sqrt(deltaA.^2.*((one1ma - rabima)./(zero1ma - one1ma).^2).^2 + deltaB.^2.*((1./(one1ma - zero1ma) + (rabima - one1ma)./(zero1ma - one1ma).^2).^2)+ deltaC.^2.*(1./(zero1ma - one1ma)).^2);

%figure(43)
%plot(experrorSIG,'b')
%hold on
%plot(errSIG,'r')
%hold off

%%% RENAME, MOVE VECTOR OF FREQUENCIES W
rabinorm = rabinorm2;
if rem(lagpts,2)
    if lagpts == 1
    %w = w(1:1:end-1);  
    else
    deltaDbarBEFma = deltaDbarBEFma(floor(lagpts/2):1:end-floor(lagpts/2)-1);
    w = w(floor(lagpts/2):1:end-floor(lagpts/2)-1);
    end
else
w = w(floor(lagpts/2):1:end-floor(lagpts/2));  
deltaDbarBEFma = deltaDbarBEFma(floor(lagpts/2):1:end-floor(lagpts/2));
end
figure(1234)
plot(w,rabinorm,'b')
hold on

figure(55)
plot(deltaDbarBEFma, 'b')
hold on
plot(deltaDbarNEW, 'r')
hold off

%%% FIT
myfun = @(p, w) p(1) * cos(2*pi*p(2)*(p(3) - w)) + p(4) + p(5)*cos(2*pi*p(2)*(p(6) + (p(3) - w))) + p(7)*cos(2*pi*p(2)*(p(6) - (p(3) - w)));
pinit = [0.2, a(aux),       3.157e9,  mean(rabinorm),  0.2, 2.1e6, 0.2];
LB = [0.1,    b(aux),      3.1569e9,             0.1,  0.1, 1.9e6, 0.1];
UB = [0.4,    c(aux),      3.1575e9,             0.7,  0.4, 2.4e6, 0.4];
[pbest,delta_p]=easyfit(w, rabinorm, pinit, myfun, LB, UB);
% fitted resonance frequency
res(aux) = pbest(3);
% fitted oscillation frequency
oscfreq(aux) = pbest(2);

% if aux == 1
%  res50 = pbest(3);%- Scan.Variable_values{2}.value;
%   w50 = w;% - Scan.Variable_values{2}.value;
%   sig50 = rabinorm;
%   fit50 = myfun(pbest,w);
%   save('NEW50cyclespi.mat','w50','sig50','fit50','res50')
% end
% 
% if aux == 3
%   res87 = pbest(3);%- Scan.Variable_values{2}.value;
%   w87 = w;% - Scan.Variable_values{2}.value;
%   sig87 = rabinorm;
%   fit87 = myfun(pbest,w);
%   save('NEW87cyclespi.mat','w87','sig87','fit87','res87') 
% end
% 
% if aux == 5
%   res164 = pbest(3);%- Scan.Variable_values{2}.value;
%   w164n = 3.1568*1e9:10:3.1576*1e9; %w;% - Scan.Variable_values{2}.value;
%   w164 = w;
%   sig164 = rabinorm;
%   fit164 = myfun(pbest,w164n);
%   save('NEW164cyclespi.mat','w164','sig164','fit164','res164','w164n') 
% end

%if aux==5
 %  figure(56)
  % wnov = 3.1568*1e9:10:3.1576*1e9;
   %plot(wnov, myfun(pbest,wnov),'b')
%end

wnumber(aux) = length(w);

noise(aux) = mean(deltaDbarNEW);
Arms(aux) = (max(myfun(pbest,w)) - min(myfun(pbest,w)))/2/sqrt(2);

%%% Calculate 
Stadevi = deltaDbarNEW; %STD signal 
Stadevi = Stadevi(2:1:end); %OR Stadevi(1:1:end-1);

%3 different ways of calculating the derivative
%diffe = abs(diff(myfun(pbest,w))./diff(w)); %TEST

%diffefun = @(p,w) p(1)*2*p(2)*pi*sin(2*p(2)*pi*(p(3)-w)) + p(5)*2*p(2)*pi*sin(2*p(2)*pi*(p(3) + p(6) - w)) - p(7)*2*pi*p(2)*sin(2*p(2)*pi*(-p(3) + p(6) + w));
%diffe = diffefun(pbest,w);
%diffe = diffe(2:1:end);

diffe = abs(diff(rabinorm)./diff(w)); %ORIGINAL

%%%%%%%%%% MOV AVG FOR DIFFE TOO
if doma
lagpts2 = 2;
diffema = abs(tsmovavg(diff(rabinorm),'s',lagpts2)./diff(w)/(lagpts2));
diffema=diffema(lagpts2-1:end);
%diffema=diffema(lagpts2:end-1);
%w2 = w(floor(lagpts2/2):1:end-floor(lagpts2/2));
Stadevima = tsmovavg(Stadevi,'s',lagpts2);
Stadevima = Stadevi(lagpts2-1:1:end);
rabinorma = tsmovavg(rabinorm,'s',lagpts2);
deltaDbarNEWma = tsmovavg(deltaDbarNEW,'s',lagpts2);

if rem(lagpts2,2)
    if lagpts2 == 1
    w2 = w(1:1:end-1);
    else
    w2 = w(floor(lagpts2/2):1:end-floor(lagpts2/2)-1);
    end
else
w2 = w(floor(lagpts2/2):1:end-floor(lagpts2/2));  
end
end

%%% Find minimum within 'nn' oscilation periods
nn = 1;
 if aux == 5
    nn = 2; 
 end
ene = nn;
figure(3333)
if donoma
[min_diff,arrayposresw] = min(abs(pbest(3)-w)); %position of resonance
[min_diff,array_posmaxw] = min(abs((pbest(3)+nn/2/pbest(2))-w)); 
[min_diff,array_posminw] = min(abs((pbest(3)-nn/2/pbest(2))-w)); 
ARR = (Stadevi./abs(diffe));
ARRBEFma = deltaDbarBEFma(2:1:end)./abs(diffe);
hold on
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
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%PLOT
figure(1111)
plot(w(1:1:end-1),ARR,'b')
hold on
%plot(w2,ARRma,'r')
hline1 = line([pbest(3) ;pbest(3)],[0 ;max(ARR)]);
set(hline1,'Color','k');
hline2 = line([pbest(3)+ene/2/pbest(2) ;pbest(3)+ene/2/pbest(2)],[0 ;max(ARR)]);
set(hline2,'Color','k');
set(hline2,'LineStyle','--');
hline3 = line([pbest(3)-ene/2/pbest(2) ;pbest(3)-ene/2/pbest(2)],[0 ;max(ARR)]);
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
    %ORIGINAL
[myminNUM(aux),array_poschosenw] = min(ARR(array_posminw:1:array_posmaxw));
[myminNUMBEF(aux),array_poschosenwBEFma] = min(ARRBEFma(array_posminw:1:array_posmaxw));

plot(w(array_poschosenw+array_posminw-1),ARR(array_poschosenw+array_posminw-1),'r*')
hold off

%TEST
%[blablu,array_poschosenw] = min(ARR(array_posminw:1:array_posmaxw));
%myminNUM(aux) = ARR(arrayposresw);
%[myminNUMBEF(aux),array_poschosenwBEFma] = min(ARRBEFma(array_posminw:1:array_posmaxw));
end
if doma
[myminNUMma(aux),array_poschosenwma] = min(ARRma(array_posminwma:1:array_posmaxwma));
end

%%% Calculate stdev of the uncertainty in frequency at the chosen w that
%%% minimizes the uncertainty in frequency within 'nn' oscilation periods

%analytical, without MA
stanNUMNEW = 1* (  ( (1 - 2*rabinorm(1:1:end-1))./(2*sqrt(rabinorm(1:1:end-1).*(1-rabinorm(1:1:end-1))))  )  .* (1./(diffe) ) )   .*    (deltaDbarNEW(1:1:end-1));
vectstNUMNEW(aux) = stanNUMNEW(array_poschosenw+array_posminw-1);

%analytical, with MA
%stanNUMNEWma = lagpts*lagpts2*sqrt( (  ( (1 - 2*rabinorma(1:1:end-1))./(2*sqrt(rabinorma(1:1:end-1).*(1-rabinorma(1:1:end-1))))  )  .* (1./(diffema) ) ).^2   .*    (deltaDbarNEWma(1:1:end-1)).^2  );
%vectstNUMNEWma(aux) = stanNUMNEWma(array_poschosenwma+array_posminwma-1);

%analytical DISCRETE, w/o MA
aa = diff(w);
delw = aa(1);
stanNUMM = [];
for m = 1:1:length(rabinorm)-1
stanNUMM(m) = delw*sqrt((rabinorm(m) + rabinorm(m+1) - 2*rabinorm(m)*rabinorm(m+1))^2/(4*rabinorm(m)*(1-rabinorm(m))*(rabinorm(m+1) - rabinorm(m))^4)*deltaDbarNEW(m)^2 + (rabinorm(m)*(1-rabinorm(m)))/(rabinorm(m+1) - rabinorm(m))^4*deltaDbarNEW(m+1)^2);
end
vectstNUMM(aux) = stanNUMM(array_poschosenw+array_posminw-1);



end
 
figure(1492)

% multiply by sqrt(Total interrogation time * N), divide by gyro, convert
% to nT/sqrt(Hz) by multiplying by (2pi)^(3/2)
if donoma
etas = myminNUM*2*pi/gyro.*sqrt(IntT*40*100000)*1e6*sqrt(2*pi); 
etasBEFma = myminNUMBEF*2*pi/gyro.*sqrt(IntT*40*100000)*1e6*sqrt(2*pi); 
etaserror = vectstNUMNEW/gyro*1e6*2*pi*sqrt(2*pi).*sqrt(IntT*40*100000);
end

if doma
etasma = myminNUMma*2*pi/gyro.*sqrt(IntT*40*100000)*1e6*sqrt(2*pi); 
%etaserrorma = vectstNUMNEWma/gyro*1e9*2*pi*sqrt(2*pi).*sqrt(IntT*40*100000);
end

SNfactor = Arms./noise*pi/sqrt(6)


klklkl

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if donoma
%plot(IntT/1e-6,etas,'b',IntT/1e-6,etasBEFma,'r');
%plot(IntT/1e-6,etas,'b',IntT/1e-6,etas./SNfactor,'k',IntT/1e-6,etas./(lagpontos),'y',IntT/1e-6,etas./SNfactor./(lagpontos),'r')

%errorbar(IntT/1e-6,etas,etaserror,'b*') %ORIGINAL
errorbar(IntT/1e-6,etasma,etaserror,'b*')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on

if doma
%plot(IntT/1e-6,etasma,'r',IntT/1e-6,etasma./SNfactor,'g',IntT/1e-6,etasma/(lagpts)/(lagpts2),'m')
end
%errorbar(IntT/1e-6,etas,etaserror,'b')
%errorbar(IntT/1e-6,etasma,etaserrorma,'r')
hold on

res %resonance positions

%Theoretical sensitivity
Cpi = mean(cvectorNEW);   
Cpiworse =  min(cvectorNEW);
Cpibest = max(cvectorNEW);

tempi = 1e-6:0.1e-6:10e-6;
T2Ram =  3.4632e-6;  
realangle = pi*1.04768;
factor =1;%sqrt(6)/pi/sqrt(mean(wnumber)); %sqrt(6)/pi/1/sqrt(mean(wnumber)/mean(lagpontos));%sqrt(3)/pi; %to see if any correction needed

% in nT/sqrt(Hz)
%T2 input is Ramsey's T2
sensp = @(time,C,angle,T2) factor * sqrt(2*pi)/C/gyro./sqrt(time)*angle/2/sin(angle/2)^2*1e6.*exp((time/(T2*(angle)/2/sin(angle/2))).^2); 
%plot(tempi/1e-6,sensp(tempi,1*Cpi,1.0868*pi,5.479e-6),'b',tempi/1e-6,sensp(tempi,1*Cpiworse,1.0868*pi,5.479e-6),'--b',tempi/1e-6,sensp(tempi,1*Cpibest,1.0868*pi,5.479e-6),'--b')
plot(tempi/1e-6,sensp(tempi,1*Cpi,realangle,T2Ram),'b',tempi/1e-6,sensp(tempi,1*Cpiworse,realangle,T2Ram),'--b',tempi/1e-6,sensp(tempi,1*Cpibest,realangle,T2Ram),'--b')


axis([1 10 0 8])
ylabel('Sensitivity (\muT/\surd Hz)','FontSize',20)
xlabel('Interrogation time (\mus)','FontSize',20)

% figure(8889)
% plot(IntT/1e-6,oscfreq/oscfreq(1),'k*--')
% hold on
% plot(IntT/1e-6,frequencyscales,'rd--')

mean(wnumber)