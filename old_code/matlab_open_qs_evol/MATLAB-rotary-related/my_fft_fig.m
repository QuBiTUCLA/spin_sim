function [fftx, ffty] = my_fft_fig(file)

openfig(file,'new');
h = gcf;
line = findall(h, 'Type', 'Line'); 
x = get(line(1), 'xdata');
y = get(line(1), 'ydata');
close(gcf);

x = x*2*pi/(2*pi*1e7);
figure(1); plot(x,y);

NFFT = 2^nextpow2(length(x));

Fs = 1/abs(x(1)-x(2));                    % Sampling frequency

fftx = Fs/2*linspace(0,1,NFFT/2+1);
ffty = fft(y-mean(y),NFFT)/length(x);

ffty = ffty(1:NFFT/2+1);

figure(2); plot(fftx,abs(ffty));

