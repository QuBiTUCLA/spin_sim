clear all

%file names: Rabi-T2-Delta-nbcycles-ptspercycle
%file = '73-10-21-50-20';
%file = '20-no-21-1000-20';
%file = '10-no-0-50-20';
%file = '9-no-23-4mus-20';
%file = 'test';
file = '10-no-21-100-20';

Rabi = 10e6; %in Hz
flipangle = pi/2;

openfig(file,'new');
h = gcf;
line = findall(h, 'Type', 'Line'); 
x = get(line, 'xdata');
y = get(line, 'ydata');
close(gcf);

t = x{1};
%t = x{1}/Rabi;
rabi = y{2};
rot = y{1};

plot(t,rabi,'r',t,rot,'b')
title('red: Rabi, blue: Rotary. Time in s.')

%%%% FFT Rabi
NFFT = 2^nextpow2(length(t));
Fs = 1/abs(t(1)-t(2));                    % Sampling frequency
fftx = Fs/2*linspace(0,1,NFFT/2+1);
ffty = fft(rabi-mean(rabi),NFFT)/length(t);
ffty = ffty(1:NFFT/2+1);
%%% PLOT FFT
figure(1002); plot(fftx,abs(ffty));
title('FFT Rabi')
%[pks,locs] = max(abs(ffty));

clear NFFT Fs fftx ffty 

%%%% FFT Rot
NFFT = 2^nextpow2(length(t));
Fs = 1/abs(t(1)-t(2));                    % Sampling frequency
fftx = Fs/2*linspace(0,1,NFFT/2+1);
ffty = fft(rot-mean(rot),NFFT)/length(t);
ffty = ffty(1:NFFT/2+1);
%%% PLOT FFT
figure(1003); plot(fftx,abs(ffty));
title('FFT Rot')
%[pks,locs] = max(abs(ffty));