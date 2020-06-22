clear all

%14 cycles.
%%%% using the one below
name{1} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-06-07-161630Manuelzao14cycles10avg.mat';

%name{1} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-06-07-162636Manuelzao14cycles12avgs300kreps.mat';
%below pb with avg #13, have to remove it from analysis
%name{1} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-06-06-165828Manuelzao14cycles20avgs.mat';
%name{1} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-06-07-184846Manuelzao14cycles39avgs.mat';

%42 cyckes
%name{2} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-06-07-213819Manuelzao42cycles6avgs.mat';
%using the one below
name{2} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-06-08-074640Manuelzao42cycles10avgs.mat';

%name{2} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-06-08-080035Manuelzao42cycles12avgs.mat';
%name{2} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-06-08-085622Manuelzao42cycles20avgs.mat';
%name{2} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-06-08-110825Manuelzao42cycles39avgs.mat';

%67 cycles
name{3} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-06-09-113905Manuelzao67cycles10avgs.mat';

%99 cycles
name{4} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-06-09-152313Manuelzao99cycles.mat';

%Real angle
realangle = pi*0.98436; 
%Rabi is 18.93MHz, giving true pi pulse = 26.4131ns. Here we're using pi
%pulse = 26ns (Rabi 19.23MHz).

% gyromagnetic ratio
gyro = 1.760859708e11; % s-1 T-1

%all plots taken centered around presumed resonance 3.024GHz

for aux = 1:4 %1:1:4
    
%%% LOAD
clear ARR pbest delta_p Stadevi diffe signal rabi lowref oneref zeroref w array_pos arrp Stadevibestfitsignal deriva pbest experrorONE experrorZERO experrorSIG deltaDbarNEW deltaA deltaB deltaC n0 n1 rabinorm w
load(name{aux});

%IntT = complete cycle time * [nbcycle] = 52e-9*[14 ...];
IntT(aux) = 2*Scan.Variable_values{7}.value*Scan.Variable_values{9}.value;

%if aux == 1
%Avgs(aux) = Scan.Averages - 1;
%else
Avgs(aux) = Scan.Averages;
%end
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
        if aux == 1
          if na ~= 13
              sig2{na}{k}=Scan.ExperimentDataErrorEachAvg{1}{na}{k}.^2*Reps(aux)*(Reps(aux)-1)+Scan.ExperimentDataEachAvg{1}{na}{k}.^2; 
          end
        else
            sig2{na}{k}=Scan.ExperimentDataErrorEachAvg{1}{na}{k}.^2*Reps(aux)*(Reps(aux)-1)+Scan.ExperimentDataEachAvg{1}{na}{k}.^2;
        end
        
        %sig2{na}{k}=Scan.ExperimentDataErrorEachAvg{1}{na}{k}.^2*Reps(aux)*(Reps(aux)-1)+Scan.ExperimentDataEachAvg{1}{na}{k}.^2;
        %end
    end
end

for k=1:3
    stdd{k}=0;
    for na=1:Scan.Averages
        if aux == 1
            if na ~= 13
                stdd{k}=stdd{k}+sig2{na}{k};
            end
        else
            stdd{k}=stdd{k}+sig2{na}{k};
        end
        
        %stdd{k}=stdd{k}+sig2{na}{k};
        %end
    end
    stdd{k}=sqrt(stdd{k}-Scan.ExperimentData{1}{k}.^2)/sqrt(Reps(aux)*Avgs(aux)*(Reps(aux)*Avgs(aux)-1));
end

experrorSIG = stdd{1};
experrorONE = stdd{2};
experrorZERO = stdd{3};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load varying variable, microwave w [Hz]
w = Scan.vary_begin:((Scan.vary_end -Scan.vary_begin)/(Scan.vary_points-1)):Scan.vary_end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if aux == 3
%     lagpts = 2;
% oneref = tsmovavg(oneref, 's', lagpts);
% oneref = oneref(lagpts:1:end);
% zeroref = tsmovavg(zeroref, 's', lagpts);
% zeroref = zeroref(lagpts:1:end);
% rabi = tsmovavg(rabi, 's', lagpts);
% rabi = rabi(lagpts:1:end);
% % errors get moved too for consistency
% experrorONE = tsmovavg(experrorONE, 's', lagpts);
% experrorONE = experrorONE(lagpts:1:end);
% experrorZERO = tsmovavg(experrorZERO, 's', lagpts);
% experrorZERO = experrorZERO(lagpts:1:end);
% experrorSIG = tsmovavg(experrorSIG, 's', lagpts);
% experrorSIG = experrorSIG(lagpts:1:end);
% 
% if rem(lagpts,2)
%     if lagpts == 1
%     w = w(1:1:end-1);  
%     else
%     w = w(floor(lagpts/2):1:end-floor(lagpts/2)-1);
%     end
% else
% w = w(floor(lagpts/2):1:end-floor(lagpts/2));  
% end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalized signal
rabinorm = (rabi - oneref)./(zeroref - oneref); %rabi normalized

deltaA = experrorZERO; 
deltaB = experrorONE;
deltaC = experrorSIG;
% deltaDbarNEW below is the stdev of the normalized signal
deltaDbarNEW = sqrt(deltaA.^2.*((oneref - rabi)./(zeroref - oneref).^2).^2 + deltaB.^2.*((1./(oneref - zeroref) + (rabi - oneref)./(zeroref - oneref).^2).^2)+ deltaC.^2.*(1./(zeroref - oneref)).^2);



%%%%%%%%%%%%%%%%%%%%%%%
%MV AVG
% if aux == 4 || aux == 3
% lagpts = 2;
% rabinorm = tsmovavg(rabinorm, 's', lagpts);
% rabinorm = rabinorm(lagpts:1:end);
% w = w(floor(lagpts/2):1:end-floor(lagpts/2));  
% end
%%%%%%%%%%%%%%%%%%%%%%%

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
pinit = [0.2,     a(aux),      3.024e9,   mean(rabinorm),   0, 2.17e6,   0,     0];
LB =    [0.01,    b(aux),      3.022e9,             0.01,   0, 1.9e6,   0      0];
UB =    [0.4,     c(aux),      3.026e9,              0.7,   0, 2.3e6,   0,     0];
[pbest,delta_p]=easyfit(w, rabinorm, pinit, myfun, LB, UB);
hold off


n0 = mean(zeroref)*1000*200e-9; %now SPD acquiring for 200ns
n1 = mean(oneref)*1000*200e-9; %now SPD acquiring for 200ns
HF = 2*pi*pbest(6);
time = IntT(aux);
cvectorOLD(aux) = sqrt(1/(0.5*(3 + 4/(n1-n0) + cos(realangle) + 8*n0/(n0-n1)^2/sin(realangle/2)^2)));
%cvectorNEW(aux) = (1/18).*((n0+(-1).*n1).^(-2).*(1+2.*cos(HF.*realangle.^(-1).*time.*sin((1/2).*realangle))).^(-4).*((1/36).*(3.*n0.*(1+n0)+n1+n1.^2+(n0+(-1).*n1).*(1+n0+n1).*cos(realangle)).*(1+2.*cos(HF.*realangle.^(-1).*time.*sin((1/2).*realangle))).^2+(-1/1296).*(3.*n0+n1+(n0+(-1).*n1).*cos(realangle)).^2.*(1+2.*cos(HF.*realangle.^(-1).*time.*sin((1/2).*realangle))).^4).*csc((1/2).*realangle).^2).^(-1/2);
cvectorNEW(aux) = sqrt(1/(1.5 + (-11*n0 + 5*n1)/(2*(n0-n1)^2) +cos(realangle)/2*(1 - (n1+n0)/(n0-n1)^2) + 8*n0/((n0-n1)^2*sin(realangle/2)^2)))*(1/(3/(abs(1+2*cos(2*HF*time*sin(realangle/2)/realangle)))));
conly(aux) = sqrt(1/(1.5 + (-11*n0 + 5*n1)/(2*(n0-n1)^2) +cos(realangle)/2*(1 - (n1+n0)/(n0-n1)^2) + 8*n0/((n0-n1)^2*sin(realangle/2)^2)));

C_A(aux) = (1/(3/(abs(1+2*cos(2*HF*time*sin(realangle/2)/realangle)))));

ctzero(aux) = n0;
ctone(aux) = n1;
AHF(aux) = HF;

%CHECK if fitting good
% delta_p(1)/pbest(1) %PB %1st cos amp
% delta_p(2)/pbest(2)
% delta_p(3)/pbest(3) %PB %resonance freq
% delta_p(4)/pbest(4)
% delta_p(5)/pbest(5) %PB %2nd cos amp
% delta_p(6)/pbest(6) %PB %HF
% delta_p(7)/pbest(7) %PB %3rd cos amp
% delta_p(8)/pbest(8) %PB %phase


%fitted hf
agaf(aux) = 2*pi*pbest(6);
% fitted resonance frequency
res(aux) = pbest(3);
% fitted oscillation frequency
oscfreq(aux) = pbest(2)/2/pi;
% CHECK if fitted osci freq corresponds to expected; below nb should be
% close to 1
%(pbest(2)/2/pi)/(2*IntT*sin(realangle/2)/realangle)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STD SIGNAL (NUMERATOR)

% NUMERICAL
Stadevi = deltaDbarNEW; %STD signal 
%Stadevi = abs(rabinorm - myfun(pbest,w));
Stadevi = Stadevi(2:1:end); %OR Stadevi(1:1:end-1);

% if aux == 4 || aux == 3
%     Stadevi = Stadevi(floor(lagpts/2):1:end-floor(lagpts/2));  
%     deltaDbarNEW = deltaDbarNEW(floor(lagpts/2):1:end-floor(lagpts/2));  
% end

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

% %%%%%%%%
% frequencyscales(aux) = Scan.Variable_values{9}.value/14;
% varystep(aux) =  Scan.vary_step_size;
% initnumber =10; 
% lagpontos(aux) = round((initnumber-1)*100000./frequencyscales(aux)./varystep(aux)  +1); %this the same percentage of mvg avg per cycle
% %%%%%%%%%%
% lagpts = lagpontos(aux);
% diffe = abs(tsmovavg(diff(rabinorm),'s',lagpts)./diff(w));
% diffe=diffe(lagpts-1:end);
% 
% %move everything with lag pts
% rabinorm = tsmovavg(rabinorm,'s', lagpts);
% rabinorm = rabinorm(lagpts-1:1:end);
% Stadevi = tsmovavg(Stadevi,'s', lagpts);
% Stadevi = Stadevi(lagpts-1:1:end);
% deltaDbarNEW = tsmovavg(deltaDbarNEW,'s', lagpts);
% deltaDbarNEW = deltaDbarNEW(lagpts-1:1:end);
% 
% if rem(lagpts,2)
%     if lagpts == 1
%     w = w(1:1:end-1);  
%     else
%     w = w(floor(lagpts/2):1:end-floor(lagpts/2)-1);
%     end
% else
% w = w(floor(lagpts/2):1:end-floor(lagpts/2));  
% end

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
%plot(w,ARR,'b')
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

%if aux == 1 || aux == 2
stanNUMNEW = (  ( (1 - 2*rabinorm(1:1:end-1))./(2*sqrt(rabinorm(1:1:end-1).*(1-rabinorm(1:1:end-1))))  )  .* (1./(diffe) ) )   .*    (deltaDbarNEW(1:1:end-1)); %ORIGINAL
%else
% fuc teste
%fuc = myfun(pbest,w);
%stanNUMNEW = (  ( (1 - 2*fuc(1:1:end-1))./(2*sqrt(fuc(1:1:end-1).*(1-fuc(1:1:end-1))))  )  .* (1./(diffe) ) )   .*    (deltaDbarNEW(1:1:end-1));
%stanNUMNEW = (  ( (1 - 2*fuc(2:1:end))./(2*sqrt(fuc(2:1:end).*(1-fuc(2:1:end))))  )  .* (1./(diffe) ) )   .*    (deltaDbarNEW(2:1:end));
%end
vectstNUMNEW(aux) = stanNUMNEW(array_poschosenw+array_posminw-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END Errorbars

figure(2459)
plot(w,rabinorm,'b')
hold on
errorbar(w,rabinorm,deltaDbarNEW,'k*')
%plot(w,rabinorm(1:1:end-1),'b')
hold on
%errorbar(w(array_poschosenw+array_posminw-1),rabinorm(array_poschosenw+array_posminw-1),stanNUMNEW(array_poschosenw+array_posminw-1))
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


%% save results
% % 
%  if aux == 1
% %  res14 = pbest(3);%- Scan.Variable_values{2}.value;
% %   w14 = w;% - Scan.Variable_values{2}.value;
% %   sig14 = rabinorm;
% %   fit14 = myfun(pbest,w);
% %   save('NUEVO14cyclespi.mat','w14','sig14','fit14','res14')
% err14 = mean(deltaDbarNEW);
% save('npi14.mat','err14')
%  end
% % 
%  if aux == 2
% %  res42 = pbest(3);%- Scan.Variable_values{2}.value;
% %   w42 = w;% - Scan.Variable_values{2}.value;
% %   sig42 = rabinorm;
% %   fit42 = myfun(pbest,w);
% %   save('NUEVO42cyclespi.mat','w42','sig42','fit42','res42')
% err42 = mean(deltaDbarNEW);
% save('npi42.mat','err42')
%  end
% % 
%  if aux == 3
% %  res67 = pbest(3);%- Scan.Variable_values{2}.value;
% %   w67 = w;% - Scan.Variable_values{2}.value;
% %   sig67 = rabinorm;
% %   fit67 = myfun(pbest,w);
% %   save('NUEVO67cyclespi.mat','w67','sig67','fit67','res67')
% err67 = mean(deltaDbarNEW);
% save('npi67.mat','err67')
%  end
% % 
%  if aux == 4
% %  res99 = pbest(3);%- Scan.Variable_values{2}.value;
% %   w99 = w;% - Scan.Variable_values{2}.value;
% %   sig99 = rabinorm;
% %   fit99 = myfun(pbest,w);
% %   save('NUEVO99cyclespi.mat','w99','sig99','fit99','res99')
% err99 = mean(deltaDbarNEW);
% save('npi99.mat','err99')
%  end



end

figure(1492)
errorbar(IntT/1e-6,etas,etaserror,'k*') %ORIGINAL
hold on

%Theoretical sensitivity
Cpi = mean(cvectorNEW) %5.9e-3
%std(cvectorNEW) %1.4e-3
%std(cvectorNEW)
Cpiworse =  min(cvectorNEW);
Cpibest = max(cvectorNEW);
%COLD = mean(cvectorOLD)
%std(cvectorNEW)
%Cpibest - Cpi
%Cpi - Cpiworse

tempi = 0:0.01e-6:6e-6;

%in use was below
%tempi = 0:52e-9:6e-6;
T2Ram = 2.19e-6; % as per fitting in plotramseymanuelzao
%old values, for Zorba Ramsey: %3.0582e-6;  %1.9096e-6; 

% in nT/sqrt(Hz)
%T2 input is Ramsey's T2
sensp = @(time,C,angle,T2) sqrt(2*pi)/C/gyro./sqrt(time)*angle/2/sin(angle/2)^2*1e6.*exp((time/(T2*(angle)/2/sin(angle/2))).^2); 
plot(tempi/1e-6,sensp(tempi,Cpi,realangle,T2Ram),'b',tempi/1e-6,sensp(tempi,Cpiworse,realangle,T2Ram),'--b',tempi/1e-6,sensp(tempi,Cpibest,realangle,T2Ram),'--b')
hold on
HF = mean(agaf);
HFmin = min(agaf);
HFmax = max(agaf);
%HF/2/pi
%std(agaf/2/pi)
Cpin = mean(conly); %mean(cvectorNEW);
%std(cvectorNEW)
Cpiworsen =  min(conly); %min(cvectorNEW);
Cpibestn = max(conly); %max(cvectorNEW);
sensp2 = @(time,C,angle,T2,HFi) sqrt(2*pi)/C./(1./(3./(abs(1+2*cos(2*HFi*time*sin(realangle/2)/realangle)))))/gyro./sqrt(time)*angle/2/sin(angle/2)^2*1e6.*exp((time/(T2*(angle)/2/sin(angle/2))).^2); 
%plot(tempi/1e-6,sensp2(tempi,Cpin,realangle,T2Ram,HF),'r',tempi/1e-6,sensp2(tempi,Cpiworsen,realangle,T2Ram,HF),'--r',tempi/1e-6,sensp2(tempi,Cpibestn,realangle,T2Ram,HF),'--r')
plot(tempi/1e-6,sensp2(tempi,Cpin,realangle,T2Ram,HF),'r')

%plot(tempi/1e-6,sensp(tempi,Cpi+std(cvectorNEW),realangle,T2Ram),'--m',tempi/1e-6,sensp(tempi,Cpi-std(cvectorNEW),realangle,T2Ram),'--m')
hold on
%plot all points individually
 pointC1 = sensp(IntT(1),cvectorNEW(1),realangle,T2Ram);
 pointC2 = sensp(IntT(2),cvectorNEW(2),realangle,T2Ram);
 pointC3 = sensp(IntT(3),cvectorNEW(3),realangle,T2Ram);
 pointC4 = sensp(IntT(4),cvectorNEW(4),realangle,T2Ram);
% %%%%%
 %plot(IntT(1)/1e-6,pointC1,'r*')
 %plot(IntT(2)/1e-6,pointC2,'r*')
 %plot(IntT(3)/1e-6,pointC3,'r*')
 %plot(IntT(4)/1e-6,pointC4,'r*')
%plot(tempi/1e-6,sensp(tempi,COLD,realangle,T2Ram),'g')

axis([0.1 5.5 1 30])
ylabel('Sensitivity (\muT/\surd Hz)','FontSize',20)
xlabel('Interrogation time (\mus)','FontSize',20)

 %figure(8889)
 %plot(IntT/1e-6,oscfreq/oscfreq(1),'k*--')
 %plot(IntT/1e-6,oscfreq,'k*--')
 %hold on
 %plot(IntT/1e-6,2*IntT*sin(realangle/2)/realangle,'g*--')
%   
 % figure(8809)
 % plot(IntT/1e-6,res-3.024e9,'k*--')
%  

%
C_A
mean(C_A)
std(C_A)

%res
mean(res)
std(res)
%oscfreq

%mean(ctzero)
%std(ctzero)

%mean(ctone)
%std(ctone)

%mean(AHF)/2/pi
%std(AHF)/2/pi