%3pi/4
clear all

IntT = 42e-9*[35 49 75 90 100];

name{1} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-23-212840Uno3pi435cycles.mat';
name{2} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-18-113453Uno3pi449cycles.mat';
name{3} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-19-145847Uno3pip475cycles.mat';
name{4} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-18-072131Uno3pi490cycles.mat';
name{5} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-18-155453Uno3pi4100cyclesnotverygood.mat';

a = 1e-6*[1.5, 2, 3, 3.5, 4];
b = 1e-6*[0.5, 1, 2, 3, 3.5];
c = 1e-6*[2, 3, 4, 4, 4.5];

for aux = 1:1:5
    
    clear Stadevi diffe signal rabi lowref oneref zeroref w array_pos arrp Stadevibestfitsignal deriva pbest
    load(name{aux});

rabi = Scan.ExperimentData{1}{1};
lowref = Scan.ExperimentData{1}{2};
oneref = Scan.ExperimentData{1}{3};
zeroref = Scan.ExperimentData{1}{4};
experrorONE = Scan.ExperimentDataError{1}{3}/sqrt(100000); %doing moving avg for repetitions as well; formula probably not entirely correct, but order of magnitude ok
experrorZERO = Scan.ExperimentDataError{1}{4}/sqrt(100000);  %doing moving avg for repetitions as well
experrorSIG = Scan.ExperimentDataError{1}{1}/sqrt(100000);  %doing moving avg for repetitions as well
w = Scan.vary_begin:((Scan.vary_end -Scan.vary_begin)/(Scan.vary_points-1)):Scan.vary_end;
rabinorm = (rabi - oneref)./(zeroref - oneref); %rabi normalized

figure(13)
plot(w,rabi,w,zeroref,w,oneref)
mean(zeroref)
mean(oneref)

% with moving average
%comment if not wanted
lagpts = 4; %needs to be an even number, >=2

one1ma = tsmovavg(oneref, 's', lagpts);
one1ma = one1ma(lagpts:1:end);
zero1ma = tsmovavg(zeroref, 's', lagpts);
zero1ma = zero1ma(lagpts:1:end);
rabima = tsmovavg(rabi, 's', lagpts);
rabima = rabima(lagpts:1:end);

%all vectors -lagpts long
%%%%% Experimenta do get moved too????
errONE = tsmovavg(experrorONE, 's', lagpts);
errONE = errONE(lagpts:1:end);
errZERO = tsmovavg(experrorZERO, 's', lagpts);
errZERO = errZERO(lagpts:1:end);
errSIG = tsmovavg(experrorSIG, 's', lagpts);
errSIG = errSIG(lagpts:1:end);

rabinorm2 = (rabima - one1ma)./(zero1ma - one1ma); %rabi normalized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltaA = errZERO; %already corrected for std dev of mean for repetitions
deltaB = errONE;
deltaC = errSIG;
deltaDbarNEW = sqrt(deltaA.^2.*((one1ma - rabima)./(zero1ma - one1ma).^2).^2 + deltaB.^2.*((1./(one1ma - zero1ma) + (rabima - one1ma)./(zero1ma - one1ma).^2)).^2 + deltaC.^2.*(1./(zero1ma - one1ma)).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(120); 
plot(w,rabinorm,'b',w(floor(lagpts/2):1:end-floor(lagpts/2)),rabinorm2,'r*')

rabinorm = rabinorm2;
w = w(floor(lagpts/2):1:end-floor(lagpts/2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(12); 
plot(w,rabinorm)
hold on
myfun = @(p, w) p(1) * cos(2*pi*p(2)*(p(3) - w)) + p(4) + p(5)*cos(2*pi*p(2)*(p(6) + (p(3) - w))) + p(7)*cos(2*pi*p(2)*(p(6) - (p(3) - w)));
%p(1) amp
%p(2) freq of osci
%p(3) reson freq
%p(4) offset
%p(5) amp
%p(6) hyperfine
%p(7) amp
% initial values 
pinit = [0.2, a(aux), 3.157e9,mean(rabinorm), 0.2, 2.1e6, 0.2];
% bounds for fitting parameters 
LB = [0.1,   b(aux),      3.1565e9,     0.1,  0.1, 1.9e6, 0.1];
UB = [0.4,    c(aux),       3.1575e9,     0.9,  0.4, 2.4e6, 0.4];
[pbest,delta_p]=easyfit(w, rabinorm, pinit, myfun, LB, UB);
hold off

%ANALYTICAL
bestfitsignal = pbest(1) * cos(2*pi*pbest(2)*(pbest(3) - w)) + pbest(4) + pbest(5)*cos(2*pi*pbest(2)*(pbest(6) + (pbest(3) - w))) + pbest(7)*cos(2*pi*pbest(2)*(pbest(6) - (pbest(3) - w)));
Stadevibestfitsignal = sqrt(bestfitsignal.*(1-bestfitsignal));
deriva = pbest(1)*2*pi*pbest(2)*sin(2*pi*pbest(2)*(pbest(3)-w)) + pbest(5)*2*pi*pbest(2)*sin(2*pi*pbest(2)*(pbest(6) + (pbest(3)-w))) - pbest(7)*2*pi*pbest(2)*sin(2*pi*pbest(2)*(pbest(6) - (pbest(3)-w)));
[fulletaA(aux), array_pos] = min((Stadevibestfitsignal./abs(deriva))*sqrt(IntT(aux)));
maxvectorfulletaA(aux) = w(array_pos);  
  
figure(34)

%how many wavelengths from fitted resonance? (should want 1)
ene = 1;

plot(w,(((Stadevibestfitsignal./abs(deriva))*sqrt(IntT(aux)))))
hold on
hline1 = line([pbest(3) ;pbest(3)],[0 ;max((((Stadevibestfitsignal./abs(deriva))*sqrt(IntT(aux)))))]);
set(hline1,'Color','r');
hline2 = line([pbest(3)+ene/2/pbest(2) ;pbest(3)+ene/2/pbest(2)],[0 ;max((((Stadevibestfitsignal./abs(deriva))*sqrt(IntT(aux)))))]);
set(hline2,'Color','r');
set(hline2,'LineStyle','--');
hline3 = line([pbest(3)-ene/2/pbest(2) ;pbest(3)-ene/2/pbest(2)],[0 ;max((((Stadevibestfitsignal./abs(deriva))*sqrt(IntT(aux)))))]);
set(hline3,'Color','r');
set(hline3,'LineStyle','--');

[min_diff,arrayposresw] = min(abs(pbest(3)-w)); %position of resonance
[min_diff,array_posmaxw] = min(abs((pbest(3)+ene/2/pbest(2))-w));
[min_diff,array_posminw] = min(abs((pbest(3)-ene/2/pbest(2))-w));
ARR = (((Stadevibestfitsignal./abs(deriva))*sqrt(IntT(aux))));
if array_posmaxw > length(ARR)
    array_posmaxw = length(ARR);
end
[myminANA(aux),array_poschosenw] = min(ARR(array_posminw:1:array_posmaxw));
plot(w(array_poschosenw+array_posminw-1),ARR(array_poschosenw+array_posminw-1),'*')
hold off

%%%%%%%%%%%%%%%%%%%%%%
%standard deviation for the sensitivity, calculated at the chosen w
stanANA = sqrt( (  ( (1 - 2*bestfitsignal)./(2*sqrt(bestfitsignal.*(1-bestfitsignal)))  )  .* (1./(deriva) ) ).^2   .*    (deltaDbarNEW).^2  );
vectstANA(aux) = stanANA(array_poschosenw+array_posminw-1);

%NUMERICAL
Stadevi = deltaDbarNEW; %STD signal                 %S(S-1), shot noise limited,probably not correct ------ sqrt(rabinorm.*(1-rabinorm));
Stadevi = Stadevi(2:1:end); %Stadevi(1:1:end-1);
diffe = diff(rabinorm)./diff(w);
[fulletaN(aux), array_pos] = min((Stadevi./abs(diffe))*sqrt(IntT(aux)));
maxvectorfulletaN(aux) = w(array_pos); 

figure(33)
plot(w(1:1:end-1),(((Stadevi./abs(diffe))*sqrt(IntT(aux)))))
hold on
hline1 = line([pbest(3) ;pbest(3)],[0 ;max((((Stadevi./abs(diffe))*sqrt(IntT(aux)))))]);
set(hline1,'Color','r');
hline2 = line([pbest(3)+ene/2/pbest(2) ;pbest(3)+ene/2/pbest(2)],[0 ;max((((Stadevi./abs(diffe))*sqrt(IntT(aux)))))]);
set(hline2,'Color','r');
set(hline2,'LineStyle','--');
hline3 = line([pbest(3)-ene/2/pbest(2) ;pbest(3)-ene/2/pbest(2)],[0 ;max((((Stadevi./abs(diffe))*sqrt(IntT(aux)))))]);
set(hline3,'Color','r');
set(hline3,'LineStyle','--');

[min_diff,arrayposresw] = min(abs(pbest(3)-w)); %position of resonance
[min_diff,array_posmaxw] = min(abs((pbest(3)+ene/2/pbest(2))-w)); 
[min_diff,array_posminw] = min(abs((pbest(3)-ene/2/pbest(2))-w));  
ARR = (((Stadevi./abs(diffe))*sqrt(IntT(aux))));
if array_posmaxw > length(ARR)
    array_posmaxw = length(ARR);
end
[myminNUM(aux),array_poschosenw] = min(ARR(array_posminw:1:array_posmaxw));
plot(w(array_poschosenw+array_posminw-1),ARR(array_poschosenw+array_posminw-1),'*')
hold off

%%%%%%%%%%%%%%%%%%%%%%
%with std as vector, not number
%stanNUMNEW = sqrt( (  ( (1 - 2*rabinorm(2:1:end))./(2*sqrt(rabinorm(2:1:end).*(1-rabinorm(2:1:end))))  )  .* (1./(diffe) ) ).^2   .*    (deltaDbarNEW(2:1:end)).^2  );
stanNUMNEW = sqrt( (  ( (1 - 2*rabinorm(1:1:end-1))./(2*sqrt(rabinorm(1:1:end-1).*(1-rabinorm(1:1:end-1))))  )  .* (1./(diffe) ) ).^2   .*    (deltaDbarNEW(1:1:end-1)).^2  );
vectstNUMNEW(aux) = stanNUMNEW(array_poschosenw+array_posminw-1);

aa = diff(w);
delw = aa(1);
stanNUMMasashi = [];
for m = 1:1:length(rabinorm)-1
stanNUMMasashi(m) = delw*sqrt((rabinorm(m) + rabinorm(m+1) - 2*rabinorm(m)*rabinorm(m+1))^2/(4*rabinorm(m)*(1-rabinorm(m))*(rabinorm(m+1) - rabinorm(m))^4)*deltaDbarNEW(m)^2 + (rabinorm(m)*(1-rabinorm(m)))/(rabinorm(m+1) - rabinorm(m))^4*deltaDbarNEW(m+1)^2);
end
vectstNUMMasashi(aux) = stanNUMMasashi(array_poschosenw+array_posminw-1);

end

%those plots are mostly wrong bc they do not multiply by sqrt(Reps X Avgs)
% figure(67)
% subplot(3,1,1)
% plot(IntT/1e-6,fulletaA,'ro',IntT/1e-6,fulletaN,'bo',IntT/1e-6,myminANA,'r*',IntT/1e-6,myminNUM,'b*')
% title('full sensitivity eta = stdev / derivative * sqrt(T) [sqrt(Hz)], red fit, blue num // o whole range, * 1 wavelength')
% xlabel('\mus')
% %axis([44 170 0 50])
% subplot(3,1,2)
% plot(IntT/1e-6,fulletaA/(28024.9540*10^6)*10^9,'ro',IntT/1e-6,fulletaN/(28024.9540*10^6)*10^9,'bo',IntT/1e-6,myminANA/(28024.9540*10^6)*10^9,'r*',IntT/1e-6,myminNUM/(28024.9540*10^6)*10^9,'b*')
% title('full sensitivity eta = stdev / derivative * sqrt(T) [nT/sqrt(Hz)], red fit, blue num // o whole range, * 1 wavelength')
% subplot(3,1,3)
% errorbar(IntT/1e-6,myminANA/(28024.9540*10^6)*10^9,vectstANA/(28024.9540*10^6)*10^9,'r')
% hold on
% errorbar(IntT/1e-6,myminNUM/(28024.9540*10^6)*10^9,vectstNUMNEW/(28024.9540*10^6)*10^9,'g')
% hold off
% title('full sensitivity eta = stdev / derivative * sqrt(T) [nT/sqrt(Hz)], red fit, blue num, with std dev, 1 wavelength')
% xlabel('Rotary pi // 20 lagpts // error propag calculated after moving avg of raw signal // @ chosen w of min eta')

figure(78)
%times sqrt(RepsXAvg)
errorbar(IntT/1e-6,sqrt(100000*40)*myminNUM/(28024.9540*10^6)*10^9,vectstNUMNEW/(28024.9540*10^6)*10^9)%[nT/sqrt(Hz)]
%errorbar(IntT/1e-6,sqrt(100000*40)*myminNUM,vectstNUMNEW) %in Hz
%plot(IntT/1e-6,sqrt(100000*40)*myminNUM/(28024.9540*10^6)*10^9)
%plot(IntT/1e-6,sqrt(100000*40)*myminNUM/1e-6)
title('full sensitivity eta = stdev / derivative * sqrt(T) [nT/sqrt(Hz)], 4 lagpts, 1 wavelength')

%Masashi's formula is better
figure(790)
%times sqrt(RepsXAvg)
errorbar(IntT/1e-6,sqrt(100000*40)*myminNUM/(28024.9540*10^6)*10^9,vectstNUMMasashi/(28024.9540*10^6)*10^9,'b*')%[nT/sqrt(Hz)]
%errorbar(IntT/1e-6,sqrt(100000*40)*myminNUM,vectstNUMMasashi,'b*') %in Hz
%plot(IntT/1e-6,sqrt(100000*40)*myminNUM/(28024.9540*10^6)*10^9)
%plot(IntT/1e-6,sqrt(100000*40)*myminNUM/1e-6)
title('full sensitivity eta = stdev / derivative * sqrt(T) [nT/sqrt(Hz)], 4 lagpts, 1 wavelength')

% 
  time3pi4old= IntT/1e-6;
   sig3pi4old = sqrt(100000*40)*myminNUM/(28024.9540*10^6)*10^9;
   err3pi4old = sqrt(100000*40)*vectstNUMNEW/(28024.9540*10^6)*10^9;
  
  save('3piov4old.mat','time3pi4old','sig3pi4old','err3pi4old')

% error propagation:
% error propagation for the normalizing; 
% error propagation for the function eta = std(signal)/deriv(S,w) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
