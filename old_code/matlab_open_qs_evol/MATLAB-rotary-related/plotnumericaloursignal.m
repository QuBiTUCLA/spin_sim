load('1.mat') %delta 0
load('2.mat') %+2.1
load('3.mat') %-2.1

TOTS = (signal1 + signal2 + signal3)/3;
TOTT = tempo1;

figure(34)
plot(TOTT/1e-6,TOTS,'b')
axis([0 2 0 1])
xlabel('\mus')
title('For Rabi 17.857, Delta = 0,pm2.1MHz, up to 500cycles, 40pts per cycle')

%Calculating FT

timearray = TOTT;
NFFT = 2^nextpow2(length(timearray));
Fs = 1/abs(timearray(1)-timearray(2));                    % Sampling frequency
fftx = Fs/2*linspace(0,1,NFFT/2+1);

ffty = fft(TOTS-mean(TOTS),NFFT)/length(timearray);
ffty = ffty(1:NFFT/2+1);

%normalize ffty
mais5 = max(abs(ffty));
menos5 = min(abs(ffty));
vetorS = (abs(ffty) - menos5)/(mais5 - menos5);

figure(345)
plot(fftx/1e6,vetorS, 'g')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%try to see if adding noise adds peak
clear all

load('noisy1.mat')
load('noisy2.mat')
load('noisy3.mat')

TOTSN = (signaln1 + signaln2 + signaln3)/3;
TOTTN = tempon1;

figure(340)
plot(TOTTN/1e-6,TOTSN,'b')
axis([0 2 0 1])
xlabel('\mus')
title('For Rabi 17.857, Delta = 0,pm2.1MHz, up to 500cycles, 40 pts per cycle, avg 200 times over static noise 1% Rabi')

%Calculating FT

clear timearray NFFT Fs fftx ffty mais5 menos5 vetorS

timearray = TOTTN;
NFFT = 2^nextpow2(length(timearray));
Fs = 1/abs(timearray(1)-timearray(2));                    % Sampling frequency
fftx = Fs/2*linspace(0,1,NFFT/2+1);

ffty = fft(TOTSN-mean(TOTSN),NFFT)/length(timearray);
ffty = ffty(1:NFFT/2+1);

%normalize ffty
mais5 = max(abs(ffty));
menos5 = min(abs(ffty));
vetorS = (abs(ffty) - menos5)/(mais5 - menos5);

figure(3450)
plot(fftx/1e6,vetorS, 'g')