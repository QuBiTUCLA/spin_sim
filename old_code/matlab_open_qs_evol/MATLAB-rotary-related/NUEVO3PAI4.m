clear all

%15 cycles
name{1} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-06-10-180222ManuelzaoRot3pi415cycles.mat';

%Real angle
realangle = 3*pi/4*1.00959; 
%Rabi is 18.93MHz, giving true 3pi/4 pulse = 19.81ns. Here we're using
%3pi/4 pulse = 20ns (Rabi 18.75MHz).

% gyromagnetic ratio
gyro = 1.760859708e11; % s-1 T-1

%all plots taken centered around presumed resonance 3.024GHz

for aux = 1:1:1
    
%%% LOAD
clear ARR pbest delta_p Stadevi diffe signal rabi lowref oneref zeroref w array_pos arrp Stadevibestfitsignal deriva pbest experrorONE experrorZERO experrorSIG deltaDbarNEW deltaA deltaB deltaC n0 n1 rabinorm w
load(name{aux});

%IntT = complete cycle time * [nbcycle] = 52e-9*[14 ...];
IntT(aux) = 2*Scan.Variable_values{7}.value*Scan.Variable_values{9}.value;
Avgs(aux) = Scan.Averages;
Reps(aux) =  Scan.Repetitions;

% initial, lower and higher bounds for oscillation frequency fitting (param called p(2))
a(aux) = 2*pi*2*IntT(aux)*sin(realangle/2)/realangle;
b(aux) = 2*pi*1.9*IntT(aux)*sin(realangle/2)/realangle; %1.9
c(aux) = 2*pi*2.1*IntT(aux)*sin(realangle/2)/realangle; %2.1

% load signal (called 'rabi'), 0 ref, 1 ref, errors;
% the errors need to be /srqt(Repetitions) bc at the time of the exps we
% were not yet using the std dev of the mean yet
rabi = Scan.ExperimentData{1}{1};
oneref = Scan.ExperimentData{1}{2};
zeroref = Scan.ExperimentData{1}{3};

experrorONE = Scan.ExperimentDataError{1}{2}; 
experrorZERO = Scan.ExperimentDataError{1}{3};
experrorSIG = Scan.ExperimentDataError{1}{1}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Paola's correction

for k=1:3
    for na=1:Scan.Averages
       sig2{na}{k}=Scan.ExperimentDataErrorEachAvg{1}{na}{k}.^2*Reps(aux)*(Reps(aux)-1)+Scan.ExperimentDataEachAvg{1}{na}{k}.^2;    
    end
end

for k=1:3
    stdd{k}=0;
    for na=1:Scan.Averages
        stdd{k}=stdd{k}+sig2{na}{k};
    end
    stdd{k}=sqrt(stdd{k}-Scan.ExperimentData{1}{k}.^2)/sqrt(Reps(aux)*Avgs(aux)*(Reps(aux)*Avgs(aux)-1));
end

experrorSIG = stdd{1};
experrorONE = stdd{2};
experrorZERO = stdd{3};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalized signal
rabinorm = (rabi - oneref)./(zeroref - oneref); %rabi normalized

deltaA = experrorZERO; 
deltaB = experrorONE;
deltaC = experrorSIG;
% deltaDbarNEW below is the stdev of the normalized signal
deltaDbarNEW = sqrt(deltaA.^2.*((oneref - rabi)./(zeroref - oneref).^2).^2 + deltaB.^2.*((1./(oneref - zeroref) + (rabi - oneref)./(zeroref - oneref).^2).^2)+ deltaC.^2.*(1./(zeroref - oneref)).^2);

% load varying variable, microwave w [Hz]
w = Scan.vary_begin:((Scan.vary_end -Scan.vary_begin)/(Scan.vary_points-1)):Scan.vary_end;




figure(33)
plot(w,rabinorm,'r')
hold on
%%% FIT
%myfun = @(p, w) p(1) * cos(p(2)*(p(3) - w) + p(8)*pi) + p(4) +p(5)*cos(p(2)*(p(6) + (p(3) - w)) + p(8)*pi) + p(7)*cos(p(2)*(p(6) - (p(3)- w)) + p(8)*pi); %original
%pinit = [0.2,     a(aux),      3.024e9,   mean(rabinorm),   0.2, 2.1e6,  0.2,     0];
%LB =    [0.01,    b(aux),      3.022e9,             0.01,  0.01, 1.9e6, 0.01   -0.5];
%UB =    [0.4,     c(aux),      3.026e9,              0.7,   0.4, 2.4e6,  0.4,   0.5];
%test below forcing amps to be the same
myfun = @(p, w) p(1) * cos(p(2)*(p(3) - w) + p(8)*pi) + p(4) + p(1)*cos(p(2)*(p(6) + (p(3) - w)) + p(8)*pi) + p(1)*cos(p(2)*(p(6) - (p(3) - w)) + p(8)*pi) + 0*p(5) + 0*p(7);
pinit = [0.2,     a(aux),      3.024e9,   mean(rabinorm),   0, 2.1e6,   0,     0];
LB =    [0.01,    b(aux),      3.022e9,             0.01,   0, 1.9e6,   0      0];
UB =    [0.4,     c(aux),      3.026e9,              0.7,   0, 2.4e6,   0,     0];
[pbest,delta_p]=easyfit(w, rabinorm, pinit, myfun, LB, UB);
hold off

%CHECK if fitting good
% delta_p(1)/pbest(1) %PB %1st cos amp
% delta_p(2)/pbest(2)
% delta_p(3)/pbest(3) %PB %resonance freq
% delta_p(4)/pbest(4)
% delta_p(5)/pbest(5) %PB %2nd cos amp
% delta_p(6)/pbest(6) %PB %HF
% delta_p(7)/pbest(7) %PB %3rd cos amp
% delta_p(8)/pbest(8) %PB %phase

% fitted resonance frequency
res(aux) = pbest(3);
% fitted oscillation frequency
oscfreq(aux) = pbest(2);
% CHECK if fitted osci freq corresponds to expected; below nb should be
% close to 1
%(pbest(2)/2/pi)/(2*IntT*sin(realangle/2)/realangle)


n0 = mean(zeroref)*1000*200e-9; %now SPD acquiring for 200ns
n1 = mean(oneref)*1000*200e-9; %now SPD acquiring for 200ns
%cvectorNEW(aux) = sqrt(1/(0.5*(3 + 4/(n1-n0) + cos(realangle) + 8*n0/(n0-n1)^2/sin(realangle/2)^2)));
time = IntT(aux);
HF = 2*pi*pbest(6);
%cvectorNEW(aux) = (1/18).*((n0+(-1).*n1).^(-2).*(1+2.*cos(HF.*realangle.^(-1).*time.*sin((1/2).*realangle))).^(-4).*((1/36).*(3.*n0.*(1+n0)+n1+n1.^2+(n0+(-1).*n1).*(1+n0+n1).*cos(realangle)).*(1+2.*cos(HF.*realangle.^(-1).*time.*sin((1/2).*realangle))).^2+(-1/1296).*(3.*n0+n1+(n0+(-1).*n1).*cos(realangle)).^2.*(1+2.*cos(HF.*realangle.^(-1).*time.*sin((1/2).*realangle))).^4).*csc((1/2).*realangle).^2).^(-1/2);

cvectorNEW(aux) = sqrt(1/(1.5 + (-11*n0 + 5*n1)/(2*(n0-n1)^2) +cos(realangle)/2*(1 - (n1+n0)/(n0-n1)^2) + 8*n0/((n0-n1)^2*sin(realangle/2)^2)))*(1/(3/(abs(1+2*cos(2*HF*time*sin(realangle/2)/realangle)))));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STD SIGNAL (NUMERATOR)

% NUMERICAL
Stadevi = deltaDbarNEW; %STD signal 
Stadevi = Stadevi(2:1:end); %OR Stadevi(1:1:end-1);

% NUMERICAL shot-noise-limited
%Stadevi = sqrt(rabinorm.*(1-rabinorm));
%Stadevi = Stadevi(2:1:end); %OR Stadevi(1:1:end-1);

% ANALYTICAL (shot-noise limited; cannot use delta_p for full treatment bc delta_p too large)
%Stadevi = sqrt(myfun(pbest,w).*(1-myfun(pbest,w)));
%Stadevi = Stadevi(2:1:end); %OR Stadevi(1:1:end-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DERIV (DENO)

% NUMERICAL ALL THE WAY
%if aux == 1 || aux == 2
%diffe = abs(diff(rabinorm)./diff(w)); %ORIGINAL

% ANALYTICAL FOR dSignal, NUMERICAL for dw 
%else
diffe = abs(diff(myfun(pbest,w))./diff(w)); 

% ANALYTICAL
%else
%diffefun = @(p,w) p(1)*p(2)*sin(p(2)*(p(3)-w) + pi*p(8)) + p(5)*p(2)*sin(p(2)*(p(3) + p(6) - w) + pi*p(8)) - p(7)*p(2)*sin(p(2)*(-p(3) + p(6) + w) + pi*p(8));
%diffe = diffefun(pbest,w);
%diffe = diffe(2:1:end);
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END Calculate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Find minimum within 'nn' oscilation periods
nn = 1;
figure(3333)
[min_diff,arrayposresw] = min(abs(pbest(3)-w)); %position of resonance
[min_diff,array_posmaxw] = min(abs((pbest(3)+nn/2/(pbest(2)/2/pi))-w)); 
[min_diff,array_posminw] = min(abs((pbest(3)-nn/2/(pbest(2)/2/pi))-w));  
ARR = (Stadevi./abs(diffe));
plot(abs(diffe),'b')

%%%%%%%%%%%%%%%%%%%%%%%%%PLOT
figure(1111)
plot(w(1:1:end-1),ARR,'b')
hold on
hline1 = line([pbest(3) ;pbest(3)],[0 ;max(ARR)]);
set(hline1,'Color','k');
hline2 = line([pbest(3)+nn/2/(pbest(2)/2/pi) ;pbest(3)+nn/2/(pbest(2)/2/pi)],[0 ;max(ARR)]);
set(hline2,'Color','k');
set(hline2,'LineStyle','--');
hline3 = line([pbest(3)-nn/2/(pbest(2)/2/pi) ;pbest(3)-nn/2/(pbest(2)/2/pi)],[0 ;max(ARR)]);
set(hline3,'Color','k');
set(hline3,'LineStyle','--');
%%%%%%%%%%%%%%%%%%%%%%%%

if array_posmaxw > length(ARR)
    array_posmaxw = length(ARR);
end

% the uncertainty in frequency is myminNUM [Hz]
[myminNUM(aux),array_poschosenw] = min(ARR(array_posminw:1:array_posmaxw));
plot(w(array_poschosenw+array_posminw-1),ARR(array_poschosenw+array_posminw-1),'r*')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Errorbars

%%% Calculate stdev of the uncertainty in frequency at the chosen w that
%%% minimizes the uncertainty in frequency within 'nn' oscilation periods

stanNUMNEW = (  ( (1 - 2*rabinorm(1:1:end-1))./(2*sqrt(rabinorm(1:1:end-1).*(1-rabinorm(1:1:end-1))))  )  .* (1./(diffe) ) )   .*    (deltaDbarNEW(1:1:end-1)); %ORIGINAL
% fuc teste
fuc = myfun(pbest,w);
%stanNUMNEW = (  ( (1 - 2*fuc(1:1:end-1))./(2*sqrt(fuc(1:1:end-1).*(1-fuc(1:1:end-1))))  )  .* (1./(diffe) ) )   .*    (deltaDbarNEW(1:1:end-1));
%stanNUMNEW = (  ( (1 - 2*fuc(2:1:end))./(2*sqrt(fuc(2:1:end).*(1-fuc(2:1:end))))  )  .* (1./(diffe) ) )   .*    (deltaDbarNEW(2:1:end));

vectstNUMNEW(aux) = stanNUMNEW(array_poschosenw+array_posminw-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END Errorbars

figure(2459)
plot(w,rabinorm,'b')
hold on
hline1 = line([pbest(3) ;pbest(3)],[0 ;1]);
set(hline1,'Color','k');
hline2 = line([pbest(3)+nn/2/(pbest(2)/2/pi) ;pbest(3)+nn/2/(pbest(2)/2/pi)],[0 ;1]);
set(hline2,'Color','k');
set(hline2,'LineStyle','--');
hline3 = line([pbest(3)-nn/2/(pbest(2)/2/pi) ;pbest(3)-nn/2/(pbest(2)/2/pi)],[0 ;1]);
set(hline3,'Color','k');
set(hline3,'LineStyle','--');
hold off

% multiply by sqrt(Total interrogation time * N), divide by gyro, convert
% to nT/sqrt(Hz) by multiplying by (2pi)^(3/2)
etas(aux) = myminNUM(aux)*2*pi/gyro.*sqrt(IntT(aux)*Avgs(aux)*Reps(aux))*1e6*sqrt(2*pi);
etaserror(aux) = vectstNUMNEW(aux)/gyro*1e6*2*pi*sqrt(2*pi).*sqrt(IntT(aux)*Avgs(aux)*Reps(aux));

end

figure(1493)
errorbar(IntT/1e-6,etas,etaserror,'k*') %ORIGINAL
hold on

%Theoretical sensitivity
Cpi = mean(cvectorNEW);   
Cpiworse =  min(cvectorNEW);
Cpibest = max(cvectorNEW);

tempi = 0:0.1e-6:10e-6;
T2Ram =  2.19e-6; %3.0582e-6; %1.9096e-6; 

% in nT/sqrt(Hz)
%T2 input is Ramsey's T2
sensp = @(time,C,angle,T2) sqrt(2*pi)/C/gyro./sqrt(time)*angle/2/sin(angle/2)^2*1e6.*exp((time/(T2*(angle)/2/sin(angle/2))).^2); 
plot(tempi/1e-6,sensp(tempi,Cpi,realangle,T2Ram),'b',tempi/1e-6,sensp(tempi,Cpiworse,realangle,T2Ram),'--b',tempi/1e-6,sensp(tempi,Cpibest,realangle,T2Ram),'--b')
hold on
%plot all points individually
pointC1 = sensp(IntT(1),cvectorNEW(1),realangle,T2Ram);
%pointC2 = sensp(IntT(2),cvectorNEW(2),realangle,T2Ram);
%pointC3 = sensp(IntT(3),cvectorNEW(3),realangle,T2Ram);
%pointC4 = sensp(IntT(4),cvectorNEW(4),realangle,T2Ram);
%%%%%
plot(IntT(1)/1e-6,pointC1,'r*')
%plot(IntT(2)/1e-6,pointC2,'r*')
%plot(IntT(3)/1e-6,pointC3,'r*')
%plot(IntT(4)/1e-6,pointC4,'r*')

axis([0.1 7 1 50])
ylabel('Sensitivity (\muT/\surd Hz)','FontSize',20)
xlabel('Interrogation time (\mus)','FontSize',20)

%  figure(8889)
%  plot(IntT/1e-6,oscfreq/oscfreq(1),'k*--')
%   
%  figure(8809)
%  plot(IntT/1e-6,res,'k*--')
%  

%res
%oscfreq