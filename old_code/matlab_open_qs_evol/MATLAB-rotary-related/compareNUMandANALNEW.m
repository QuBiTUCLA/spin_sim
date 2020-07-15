IntT = 56e-9*[50 75 87 117 164];

name{1} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-11-201924Uno50cycles.mat';
name{2} ='1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-17-202727Uno75cycles.mat';
name{3} ='1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-17-163452Uno87cyclesincomplete.mat';%alternative
name{4} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-14-174935Uno117cycles.mat';
name{5} = '1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-15-034143Uno-164cycles.mat';

a = 1e-6*[2, 3, 3.5, 4.5, 6.5];
b = 1e-6*[1, 2, 3,   4,   6];
c = 1e-6*[3, 4, 4,   5,   7];

cvector = [];

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

% for the calculation of C
% figure(12)
% plot(w,rabi,w,zeroref,w,oneref)
zz = mean(zeroref)*1000*100e-9;
umum = mean(oneref)*1000*100e-9;

cvector(aux) = sqrt(1/(1 + 2*(zz + umum)/(zz - umum)^2));

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
%%%%% Experimental errors do get moved too, for consistency
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
deltaDbarNEW =  sqrt(deltaA.^2.*((one1ma - rabima)./(zero1ma - one1ma).^2).^2 + deltaB.^2.*((1./(one1ma - zero1ma) + (rabima - one1ma)./(zero1ma - one1ma)).^2)+ deltaC.^2.*(1./(zero1ma - one1ma)).^2);
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
UB = [0.4,    c(aux),       3.1575e9,     0.7,  0.4, 2.4e6, 0.4];
[pbest,delta_p]=easyfit(w, rabinorm, pinit, myfun, LB, UB);
hold off

%ANALYTICAL
bestfitsignal = pbest(1) * cos(2*pi*pbest(2)*(pbest(3) - w)) + pbest(4) + pbest(5)*cos(2*pi*pbest(2)*(pbest(6) + (pbest(3) - w))) + pbest(7)*cos(2*pi*pbest(2)*(pbest(6) - (pbest(3) - w)));

Stadevibestfitsignal = sqrt(bestfitsignal.*(1-bestfitsignal)); %%%%%%%%%%%%%%%%%%%%%%%%S(S-1) is shot-noise limited, so that maybe not best option; only use when fitting is bad
%full treatment below
f1 = cos(2* pbest(2) *pi * (pbest(3) - w));
f2 = -2 *pi* pbest(1)* (pbest(3) - w) .*sin(2* pbest(2)*pi* (pbest(3) - w)) - 2*pbest(5)*pi*(pbest(3) + pbest(6) - w).*sin(2*pbest(2)*pi*(pbest(3) + pbest(6) - w)) - 2*pbest(7)*pi*(-pbest(3) + pbest(6) + w).*sin(2*pbest(2)*pi*(-pbest(3) + pbest(6) + w));
f3 = -2*pbest(2)*pi*pbest(1)*sin(2*pbest(2)*pi*(pbest(3) - w)) - 2*pbest(2)*pbest(5)*pi*sin(2*pbest(2)*pi*(pbest(3) + pbest(6) - w)) + 2*pbest(2)*pbest(7)*pi*sin(2*pbest(2)*pi*(-pbest(3)+ pbest(6) + w));
f5 = cos(2*pbest(2)*pi*(pbest(3) + pbest(6) - w));
f6 = -2*pbest(2)*pbest(5)*pi*sin(2*pbest(2)*pi*(pbest(3) + pbest(6) - w)) - 2*pbest(2)*pbest(7)*pi*sin(2*pbest(2)*pi*(-pbest(3) + pbest(6) + w));
f7 = cos(2*pbest(2)*pi*(-pbest(3) + pbest(6) + w));

%Stadevibestfitsignal = sqrt(delta_p(4)^2*ones(1,length(f1)) + delta_p(1)^2*f1.^2 + delta_p(2)^2*f2.^2 + delta_p(3)^2*f3.^2 + delta_p(5)^2*f5.^2 + delta_p(6)^2*f6.^2 + delta_p(7)^2*f7.^2);

deriva = abs(pbest(1)*2*pi*pbest(2)*sin(2*pi*pbest(2)*(pbest(3)-w)) + pbest(5)*2*pi*pbest(2)*sin(2*pi*pbest(2)*(pbest(6) + (pbest(3)-w))) - pbest(7)*2*pi*pbest(2)*sin(2*pi*pbest(2)*(pbest(6) - (pbest(3)-w))));

%how many wavelengths from fitted resonance? (should want 1)
ene = 1;

%%%%% min of eta or max of sensitivity?
%%% here we find the minimum of the sensitivity
[fulletaA(aux), array_pos] = min((Stadevibestfitsignal./abs(deriva))*sqrt(IntT(aux)));
maxvectorfulletaA(aux) = w(array_pos);  

figure(34)

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

%%%here we find the maximum of the derivative
% [value,array_pos] = max(deriva);
% fulletaA(aux) = Stadevibestfitsignal(array_pos)/abs(deriva(array_pos))*sqrt(IntT(aux));%here find the maximum of the derivative only
% maxvectorfulletaA(aux) = w(array_pos);  
% 
% figure(34)
% 
% plot(w,deriva)
% hold on
% hline1 = line([pbest(3) ;pbest(3)],[0 ;max(deriva)]);
% set(hline1,'Color','r');
% hline2 = line([pbest(3)+ene/2/pbest(2) ;pbest(3)+ene/2/pbest(2)],[0 ;max(deriva)]);
% set(hline2,'Color','r');
% set(hline2,'LineStyle','--');
% hline3 = line([pbest(3)-ene/2/pbest(2) ;pbest(3)-ene/2/pbest(2)],[0 ;max(deriva)]);
% set(hline3,'Color','r');
% set(hline3,'LineStyle','--');
% 
% [min_diff,arrayposresw] = min(abs(pbest(3)-w)); %position of resonance
% [min_diff,array_posmaxw] = min(abs((pbest(3)+ene/2/pbest(2))-w));
% [min_diff,array_posminw] = min(abs((pbest(3)-ene/2/pbest(2))-w));
% ARR = deriva;
% if array_posmaxw > length(ARR)
%     array_posmaxw = length(ARR);
% end
% [value,array_poschosenw] = max(ARR(array_posminw:1:array_posmaxw));
% myminANA(aux) = (((Stadevibestfitsignal(array_poschosenw)/deriva(array_poschosenw))*sqrt(IntT(aux))));
% plot(w(array_poschosenw+array_posminw-1),ARR(array_poschosenw+array_posminw-1),'*')
% hold off

%%%%%%%%

%standard deviation for the sensitivity, calculated at the chosen w
stanANA = sqrt( (  ( (1 - 2*bestfitsignal)./(2*sqrt(bestfitsignal.*(1-bestfitsignal)))  )  .* (1./(deriva) ) ).^2   .*    (Stadevibestfitsignal).^2  );
vectstANA(aux) = stanANA(array_poschosenw+array_posminw-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NUMERICAL
Stadevi = deltaDbarNEW; %STD signal  ///S(S-1), shot noise limited,probably not correct ------ sqrt(rabinorm.*(1-rabinorm));
Stadevi = Stadevi(2:1:end); %Stadevi(1:1:end-1);
diffe = abs(diff(rabinorm)./diff(w));

%%%%% min of eta or max of sensitivity?

%%%here find the minimum of the sensitivity
[fulletaN(aux), array_pos] = min((Stadevi./diffe)*sqrt(IntT(aux))*sqrt(40*100000));

maxvectorfulletaN(aux) = w(array_pos); 

figure(33)
plot(w(1:1:end-1),(((Stadevi./abs(diffe))*sqrt(IntT(aux))*sqrt(40*100000))))
hold on
hline1 = line([pbest(3) ;pbest(3)],[0 ;max((((Stadevi./abs(diffe))*sqrt(IntT(aux))*sqrt(40*100000))))]);
set(hline1,'Color','r');
hline2 = line([pbest(3)+ene/2/pbest(2) ;pbest(3)+ene/2/pbest(2)],[0 ;max((((Stadevi./abs(diffe))*sqrt(IntT(aux))*sqrt(40*100000))))]);
set(hline2,'Color','r');
set(hline2,'LineStyle','--');
hline3 = line([pbest(3)-ene/2/pbest(2) ;pbest(3)-ene/2/pbest(2)],[0 ;max((((Stadevi./abs(diffe))*sqrt(IntT(aux))*sqrt(40*100000))))]);
set(hline3,'Color','r');
set(hline3,'LineStyle','--');

[min_diff,arrayposresw] = min(abs(pbest(3)-w)); %position of resonance
[min_diff,array_posmaxw] = min(abs((pbest(3)+ene/2/pbest(2))-w)); 
[min_diff,array_posminw] = min(abs((pbest(3)-ene/2/pbest(2))-w));  
ARR = (((Stadevi./abs(diffe))*sqrt(IntT(aux)))*sqrt(40*100000));
if array_posmaxw > length(ARR)
    array_posmaxw = length(ARR);
end
[myminNUM(aux),array_poschosenw] = min(ARR(array_posminw:1:array_posmaxw));
plot(w(array_poschosenw+array_posminw-1),ARR(array_poschosenw+array_posminw-1),'*')
hold off

%%%here find the maximum of the derivative only
% [value,array_pos] = max(diffe);
% fulletaN(aux) = Stadevi(array_pos)/abs(diffe(array_pos))*sqrt(IntT(aux))*sqrt(40*100000);
% maxvectorfulletaN(aux) = w(array_pos); 
% 
% figure(33)
% plot(w(1:1:end-1),diffe)
% hold on
% hline1 = line([pbest(3) ;pbest(3)],[0 ;max(diffe)]);
% set(hline1,'Color','r');
% hline2 = line([pbest(3)+ene/2/pbest(2) ;pbest(3)+ene/2/pbest(2)],[0 ;max(diffe)]);
% set(hline2,'Color','r');
% set(hline2,'LineStyle','--');
% hline3 = line([pbest(3)-ene/2/pbest(2) ;pbest(3)-ene/2/pbest(2)],[0 ;max(diffe)]);
% set(hline3,'Color','r');
% set(hline3,'LineStyle','--');
% 
% [min_diff,arrayposresw] = min(abs(pbest(3)-w)); %position of resonance
% [min_diff,array_posmaxw] = min(abs((pbest(3)+ene/2/pbest(2))-w)); 
% [min_diff,array_posminw] = min(abs((pbest(3)-ene/2/pbest(2))-w));  
% ARR =diffe;
% if array_posmaxw > length(ARR)
%     array_posmaxw = length(ARR);
% end
% [number,array_poschosen] = max(ARR(array_posminw:1:array_posmaxw));
% myminNUM(aux) = Stadevi(array_poschosen)/diffe(array_poschosen)*sqrt(IntT(aux))*sqrt(40*100000);
% plot(w(array_poschosenw+array_posminw-1),ARR(array_poschosenw+array_posminw-1),'*')
% hold off

%%%%%%%%%%%%%%%%%%%%%%
%with std as vector, not number
%stadevi of eta
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

% figure(11)
% errorbar(IntT/1e-6,2*pi*myminANA,2*pi*vectstANA,'r')
% hold on
% errorbar(IntT/1e-6,2*pi*myminNUM,2*pi*vectstNUMNEW,'g')
% hold on
% errorbar(IntT/1e-6,2*pi*myminNUM,2*pi*vectstNUMMasashi,'b')
% hold off
% title('full sensitivity eta = stdev / derivative * sqrt(T) [1/sqrt(s)], red fit, blue num, with std dev, 1 wavelength')
% xlabel('Rotary pi // 4 lagpts // error propag calculated after moving avg of raw signal // @ chosen w of min eta')


figure(80)
%subplot(3,1,1)
%plot(IntT/1e-6,fulletaA,'ro',IntT/1e-6,fulletaN,'bo',IntT/1e-6,myminANA,'r*',IntT/1e-6,myminNUM,'b*')
%title('full sensitivity eta = stdev / derivative * sqrt(T) [sqrt(Hz)], red fit, blue num // o whole range, * 1 wavelength')
%xlabel('\mus')
%axis([44 170 0 50])
%subplot(3,1,2)
%plot(IntT/1e-6,fulletaA/(28024.9540*10^6)*10^9,'ro',IntT/1e-6,fulletaN/(28024.9540*10^6)*10^9,'bo',IntT/1e-6,myminANA/(28024.9540*10^6)*10^9,'r*',IntT/1e-6,myminNUM/(28024.9540*10^6)*10^9,'b*')
%title('full sensitivity eta = stdev / derivative * sqrt(T) [nT/sqrt(Hz)], red fit, blue num // o whole range, * 1 wavelength')
%subplot(3,1,3)

plot(IntT/1e-6,myminANA/(28024.9540*10^6)*10^9/sqrt(2*pi),'r',IntT/1e-6,myminNUM/(28024.9540*10^6)*10^9/sqrt(2*pi),'b')

klklkl
%plot(IntT/1e-6,myminNUM/(28024.9540*10^6)*10^9,'g')
 
%errorbar(IntT/1e-6,myminANA/(28024.9540*10^6)*10^9,vectstANA/(28024.9540*10^6)*10^9,'r')
%hold on
%errorbar(IntT/1e-6,myminNUM/(28024.9540*10^6)*10^9,vectstNUMNEW/(28024.9540*10^6)*10^9,'g')
%hold on
errorbar(IntT/1e-6,2*pi*myminNUM/(28024.9540*10^6)*10^9,2*pi*vectstNUMMasashi/(28024.9540*10^6)*10^9,'b')
hold off
title('full sensitivity eta = stdev / derivative * sqrt(T) [nT/sqrt(Hz)], red fit, blue num, with std dev, 1 wavelength')
xlabel('Rotary pi // 4 lagpts // error propag calculated after moving avg of raw signal // @ chosen w of min eta')
 
  timepi= IntT/1e-6;
  sigpi = myminNUM/(28024.9540*10^6)*10^9/sqrt(2*pi)
  errpiold = vectstNUMNEW/(28024.9540*10^6)*10^9; 
  errpinew = vectstNUMMasashi/(28024.9540*10^6)*10^9; 
 %%%%% save('pisensitivityBIS','timepi','sigpi','errpiold','errpinew')

cvector

mean(cvector)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the S(w) X w plot
% if aux == 1
%  w50 = w;
%  sig50 = rabinorm;
%  fit50 = bestfitsignal;
%  res50 = arrayposresw;
%  del50 = delta_p(3);
%  save('50cyclespi.mat','w50','sig50','fit50','res50','del50')
% end
% if aux == 2
%  w75 = w;
%  sig75 = rabinorm;
%  fit75 = bestfitsignal;
%  res75 = arrayposresw;
%  del75 = delta_p(3);
%  save('75cyclespi.mat','w75','sig75','fit75','res75','del75')
% end
% if aux == 3
%  w87 = w;
%  sig87 = rabinorm;
%  fit87 = bestfitsignal;
%  res87 = arrayposresw;
%  del87 = delta_p(3);
%  save('87cyclespi.mat','w87','sig87','fit87','res87','del87')
% end
% if aux == 4
%  w117 = w;
%  sig117 = rabinorm;
%  fit117 = bestfitsignal;
%  res117 = arrayposresw;
%  del117 = delta_p(3);
%  save('117cyclespi.mat','w117','sig117','fit117','res117','del117')
% end
% if aux == 5
%  w164 = w;
%  sig164 = rabinorm;
%  fit164 = bestfitsignal;
%  res164 = arrayposresw;
%  del164= delta_p(3);
%  save('164cyclespi.mat','w164','sig164','fit164','res164','del164')
% end

%%%%%%%%%%%%%%%%%%%%%%
