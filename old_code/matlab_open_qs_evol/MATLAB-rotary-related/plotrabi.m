
% UNINDENTIFIED NV
%load('Rabi2011-09-01-080936') %%% Frequencies found: 9.924e6 MHz, 0.35MHz
%load('Rabi2011-09-02-094328') %%%% Frequencies 7.813, 8.24 and 8.362MHz
%load('Rabi2011-09-04-090055') %%% Bad S/N FFT
%load('Rabi2011-09-09-020632')

% 1st NV stripline
%load('StripTrialRabi2011-10-05-143523')
%load('StripRabi2011-10-06-063631');
%load('Rabi20mhz2011-10-28-132429'); %%%% Frequencies 10.98MHz whole rabi

% NV Nara
%load('NaraRabi-100G-10dBm');
%load('Nara1stRabi-0field-0dBm');
%load('NaraRabi-100G-minus5dBm')
%load('NaraRabi-100G-0dBm-10micros');
%load('NaraRabi-100G-0dBm-10micros-higherfreq');
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2011-11-11-093637Nara')
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2011-11-12-125332Nara')
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2011-11-14-215611Nara') %Rabi with 2dBm
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2011-11-15-151908Nara.') %Rot with 2dBm, pi-rotary
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2011-11-16-022757Nara') % New Rabi with 0dBm
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2011-11-13-041652Nara4dBm')
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2011-11-14-002708Nara2nd')
%load('NaraRabi-100G-0dBm-10to15micros');

%NV Zorba
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2011-11-17-024544Zorba-nowait')
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2011-11-17-235916Zorba-withwait') %120s wait
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2011-11-18-111314Zorba-withwait')
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2011-11-18-145833Zorba-withwait')
%all below are with 120s wait
load('1DExp-seq-Rabi-vary-length_rabi_pulse-2011-11-18-164828Zorba-0to2-5dBm') %very good
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2011-11-18-183434Zorba-2to4-5dBm')
%rotary by parts 1
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2011-11-19-024344Zorba')
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2011-11-19-104600Zorba')
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2011-11-19-132459Zorba')
%rotary by parts 2
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2011-11-19-181234Zorba-1to2-5dBm')
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2011-11-19-211823Zorba-2to3-5dBm') %3055
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2011-11-20-015512Zorba-1to2-3065')
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2011-11-20-200724Zorba-3to4')
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2011-11-21-124915Zorba-3to4-bis')
%load('1DExp-seq-ramsey-vary-tau-2011-11-21-095811Zorba')
%load('1DExp-seq-ramsey-vary-tau-2011-11-22-050833Zorba-5MHzdet')
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2011-11-24-105016Zorba-5to6') % 5 to 6 E
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2011-11-24-221641Zorba-3pi4-0to1') %3pi4 0to1 F
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2011-11-20-114120-Zorba-0to1')
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2011-11-28-113818Zorba')
%dont understand the saved exp above
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2011-11-28-153209Zorba-5pi-0to1')
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2011-11-29-183841Zorba-5pi-5to6');
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2011-12-07-105900Uno-1to6-0dBm.mat')
%echo
%load('1DExp-seq-echo_envelope-vary-tau-2011-11-23-094620Zorba-1of2.mat');

%%% NV Uno
%S(w) X w
%load('1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-09-181648Uno50cycles');
%load('1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-10-124433Uno75cycles'); 
%load('1DExp-seq-RotaryEchoCompleteEvol-SwXw-vary-mw_freq-2012-01-13-023143Uno87cycles.mat') %1.75 oscillations expected
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2012-02-07-154902Rabitocheck.mat') %up to 2mus, gives freq of 16.1MHz
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2012-02-07-181536Rabicheck5mus.mat')
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2012-02-13-150434Beli-Rabicheck5mus.mat')
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2012-02-13-220218Beli-checkRabi5mus_8dBm.mat')
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2012-02-14-051442Beli-Rabicheck5mus-10dBm.mat')

%NV Fuji
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2012-02-17-050003FujiRabicheck5mus.mat') %12.57MHz Rabi

%NV noname
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2012-04-18-213341NoNameusualRabi.mat')
%load('1DExp-seq-OUnoise_Rabi-vary-length_rabi_pulse-2012-04-18-185736NoNameModDepth1commented.mat')
%load('1DExp-seq-OUnoise_Rabi-vary-length_rabi_pulse-2012-04-19-102314NoNameModDepth08commented.mat')
%load('1DExp-seq-OUnoise_Rabi-vary-length_rabi_pulse-2012-04-19-124431NoNameModDepth08UNcommentedNoise0.mat')

%NV Mindlin
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2012-05-23-055352.mat');
%load('1DExp-seq-ramsey-vary-tau-2012-05-23-174810Ramsey5MHZMindlin.mat');
%load('1DExp-seq-ramsey-vary-tau-2012-05-24-070636Ramsey5MHzMindlin2.mat');
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2012-05-24-080032.mat');
%load('1DExp-seq-ramsey-vary-tau-2012-05-24-130916Ramsey3rd.mat');
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2012-05-24-160601-Rabi1mus.mat');
%load('1DExp-seq-ramsey-vary-tau-2012-05-25-003809Ramsey.mat');
%load('1DExp-seq-Rabimodulated-vary-length_rabi_pulse-2012-05-25-121808.mat');

%NV Hana
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2012-05-30-155146Hana.mat');

%NV Manuelzao
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2012-05-30-184854Manuelzao5dBmRabi.mat')
%load('1DExp-seq-ramsey-vary-tau-2012-05-30-234349Manuelzao5MHz.mat')
%load('1DExp-seq-Rabimodulated-vary-length_rabi_pulse-2012-05-31-114659Manuelzao1-4Rabi.mat');
%load('1DExp-seq-Rabimodulated-vary-length_rabi_pulse-2012-05-31-173958ManuelzaoStNoise10.mat');
%load('1DExp-seq-Rabi-vary-length_rabi_pulse-2012-05-31-190308Manuelzao.mat')
%load('1DExp-seq-ramsey-vary-tau-2012-06-01-002910Manuelzao5MHz.mat');
%load('1DExp-seq-Rabimodulated-vary-length_rabi_pulse-2012-06-01-112037Manuelzao-Rabinonoise4mus.mat'); %this is the standard
%load('1DExp-seq-Rabimodulated-vary-length_rabi_pulse-2012-06-03-195620Manuelzao-5pcstatic.mat');
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2012-06-08-224101ManuelzaoREPi0to310avgs.mat');
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2012-06-09-045635ManuelzaoREPi0to339avgs.mat');
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2012-06-09-005044ManuelzaoREPi0to320avgs.mat');
%load('1DExp-seq-RotaryEchoCompleteEvolution-vary-length_rotary_pulse-2012-06-10-013743ManuelzaoREpi5to6mus.mat'); 
%load('1DExp-seq-Rabimodulated-vary-length_rabi_pulse-2012-06-01-112037Manuelzao-Rabinonoise4mus.mat');

%x = x(1:1:400);
%rabinorm = rabinorm(1:1:400);
% only part above gives 10.74MHz, very good oscillations

rabi = Scan.ExperimentData{1}{1};
lowref = Scan.ExperimentData{1}{2};
oneref = Scan.ExperimentData{1}{3};
zeroref = Scan.ExperimentData{1}{4};

%error propagation
% experrorONE = Scan.ExperimentDataError{1}{3};
% experrorZERO = Scan.ExperimentDataError{1}{4};
% experrorSIG = Scan.ExperimentDataError{1}{1};
% deltaA = errZERO;
% deltaB = errONE;
% deltaDbarNEW = sqrt(deltaA.^2.*((oneref - rabi)./(zeroref - oneref).^2).^2 + deltaB.^2.*((1./(oneref - zeroref) + (rabi - oneref)./(zeroref - oneref).^2)).^2);
%%%%%%%%%%%%%

x = Scan.vary_begin:((Scan.vary_end -Scan.vary_begin)/(Scan.vary_points-1)):Scan.vary_end;
 figure(11);
 plot(x,rabi,'k',x,lowref,'r',x,oneref,'b',x,zeroref,'g') %,x,15+2*exp(-x/(5*1e-6)),'r');
 
 %average then normalize
rabinorm = (rabi - oneref)./(zeroref - oneref); %rabi normalized

%cut undesirable
%x = x(1:1:800);
%rabinorm = rabinorm(1:1:800);

hold on
%figure(100); plot(x/52e-9,rabinorm)
figure(100); plot(x,rabinorm,'b')

klklk
%figure(12);
%errorbar(x,rabinorm,deltaDbarNEW,'r')

%  save('50cycles-Uno-3157.txt','x','rabinorm','-ascii');
%  fid = fopen('50cycles-Uno-3157.txt', 'w');
%  for k = 1:length(x)
%      fprintf(fid, '%d %d\n', x(k), rabinorm(k));
%       end;
%  fclose(fid);
 
 

periodorecipe(rabinorm-mean(rabinorm),x)



NFFT = 2^nextpow2(length(x));
Fs = 1/abs(x(1)-x(2));                    % Sampling frequency
fftx = Fs/2*linspace(0,1,NFFT/2+1);
ffty = fft(rabinorm-mean(rabinorm),NFFT)/length(x);
ffty = ffty(1:NFFT/2+1);
%%% PLOT FFT
figure(1002); plot(fftx,abs(ffty));
%[pks,locs] = max(abs(ffty));

lklklk

%  figure(12);
%  plot(x,rabinorm,'r')
% title('rabi: avg then norm')

%%%% ABOVE: good code
%%%% BELOW: tests

%%%%%%%%%%%%%%%%%%%%%%%%%%%normalize then average 
rabinew = {};
for aux=1:1:length(Scan.ExperimentDataEachAvg{1})
    rabinew{aux} = (Scan.ExperimentDataEachAvg{1}{aux}{1} - Scan.ExperimentDataEachAvg{1}{aux}{3})./(Scan.ExperimentDataEachAvg{1}{aux}{4} - Scan.ExperimentDataEachAvg{1}{aux}{3});
end

dim = ndims(rabinew{1});          %get the number of dimensions for your arrays
M = cat(dim+1,rabinew{:});        %convert to a (dim+1)-dimensional matrix
meanArray = mean(M,dim+1);        %get the mean across arrays

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% normalize than average with moving average
rabinewMA = {};
lagpts =50;
for aux=1:1:length(Scan.ExperimentDataEachAvg{1})
 
    
   ONE = tsmovavg(Scan.ExperimentDataEachAvg{1}{aux}{3}, 's', lagpts);
   ONE(1:1:lagpts) = Scan.ExperimentDataEachAvg{1}{aux}{3}(1:1:lagpts);
   ZERO = tsmovavg(Scan.ExperimentDataEachAvg{1}{aux}{4}, 's', lagpts);
   ZERO(1:1:lagpts) = Scan.ExperimentDataEachAvg{1}{aux}{4}(1:1:lagpts);
   
    rabinewMA{aux} = (Scan.ExperimentDataEachAvg{1}{aux}{1} - ONE)./(ZERO - ONE);
    
end

dim = ndims(rabinewMA{1});          %get the number of dimensions for your arrays
MA = cat(dim+1,rabinewMA{:});        %convert to a (dim+1)-dimensional matrix
meanArrayMA = mean(MA,dim+1) ;    %get the mean across arrays

 %figure(13);

 %plot(x,meanArray,'r')
 %title('rabi: norm then avg')

%%%%%% teste
for aux=1:1:length(Scan.ExperimentDataEachAvg{1})
rabitest{aux} = Scan.ExperimentDataEachAvg{1}{aux}{1};
onereftest{aux} = Scan.ExperimentDataEachAvg{1}{aux}{3};
zeroreftest{aux} = Scan.ExperimentDataEachAvg{1}{aux}{4};
end

dim1 = ndims(rabitest{1});          %get the number of dimensions for your arrays
M1 = cat(dim1+1,rabitest{:});        %convert to a (dim+1)-dimensional matrix
rabi1 = mean(M1,dim1+1);        %get the mean across arrays

dim2 = ndims(onereftest{1});          %get the number of dimensions for your arrays
M2 = cat(dim2+1,onereftest{:});        %convert to a (dim+1)-dimensional matrix
one1 = mean(M2,dim2+1);        %get the mean across arrays

dim3 = ndims(zeroreftest{1});          %get the number of dimensions for your arrays
M3 = cat(dim3+1,zeroreftest{:});        %convert to a (dim+1)-dimensional matrix
zero1 = mean(M3,dim3+1);        %get the mean across arrays

rabinorm2 = (rabi1 - one1)./(zero1 - one1); %rabi normalized
 %figure(15);

  %plot(x,rabinorm2,'r')
 %title('rabi: avg then norm')


%Moving average
% % if want to use moving average
%  lagpts defined above
   one1ma = tsmovavg(one1, 's', lagpts);
   one1ma(1:1:lagpts) = one1(1:1:lagpts);
   zero1ma = tsmovavg(zero1, 's', lagpts);
   zero1ma(1:1:lagpts) = zero1(1:1:lagpts);
   rabinormMA = (rabi1 - one1ma)./(zero1ma - one1ma); %rabi normalized

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT COMPARISON
   
figure(16)
%rabinorm2 -> usual avg then norm, 'b'
%meanArray -> norm than avg
%rabinormMA -> avg then norm, moving avg, 'k'
%meanArrayMA -> norm than avg, moving avg, 'm'
plot(x,rabinorm2,'b',x,rabinormMA,'k',x,meanArrayMA,'m')
title('rabi: blue is AVG -> NORM, black is AVG -> NORM with moving avg, magenta is NORM -> AVG with moving avg')
axis([min(x) max(x) -2 2])


%%%% FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NFFT = 2^nextpow2(length(x));
Fs = 1/abs(x(1)-x(2));                    % Sampling frequency
fftx = Fs/2*linspace(0,1,NFFT/2+1);
ffty = fft(rabinorm-mean(rabinorm),NFFT)/length(x);
ffty = ffty(1:NFFT/2+1);
%%% PLOT FFT
%figure(1002); plot(fftx,abs(ffty));
%[pks,locs] = max(abs(ffty));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %%% Fit
% 
% figure(1000);
% myfun = @(p, x) p(1)*sin(2*pi*p(2)*x + p(3)).* sin(2*pi*p(4)*x + p(5)) + p(6);      
% % initial values 
% pinit = [1, 7e6, pi, 3e5, pi,  mean(rabinorm)];
% %pinit = [1, 7e6, pi, mean(rabinorm)];
% 
% 
% % bounds for fitting parameters 
% LB = [0.3, 5e6,0 , 1e5 , 0,-0.2];
% UB = [1, 10e6,2*pi, 2e6, 2*pi,1.2];
% %LB = [0.3, 5e6,0, 0];
% %UB = [1.5, 10e6,2*pi,1.7];
% 
% 
% [pbest,delta_p]=easyfit(x, rabinorm, pinit, myfun, LB, UB);
% 
% %%close all;
% figure(1001);
% subplot(2,1,1);
% plot(x, rabinorm, 'ro', x, myfun(pbest, x), 'k-');
% title('fit rabi')
% subplot(2,1,2);
% plot(x,rabinorm-myfun(pbest, x));
% title('residuals')
% text(0,-2,[func2str(myfun) sprintf('\n %d, %d, %d, %d, %d,%d', pbest(1), pbest(2),pbest(3),pbest(4),pbest(5),pbest(6))])
% %text(1,-2,[func2str(myfun) sprintf('\n %d, %d, %d , %d', pbest(1), pbest(2),pbest(3),pbest(4))])