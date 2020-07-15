%function [fftx, ffty] = my_fft_mat(matfile)
function my_fft_mat(matfile)

load(matfile);

%yref =Scan.ExperimentData{1}{1};
y = Scan.ExperimentData{1}{1};
%y = yreal./yref;
x = Scan.vary_begin:((Scan.vary_end -Scan.vary_begin)/(Scan.vary_points-1)):Scan.vary_end;

% % if want to use moving average
%   lagpts = 20;
%   ymvavg = tsmovavg(yref, 's', lagpts);
%   y(lagpts:1:end) = yreal(lagpts:1:end)./ymvavg(lagpts:1:end);
% % if not just comment

%figure(1); plot(x,y); %with moving avg
figure(3); plot(x,y); %without moving avg
%figure(5);plot(x,yreal,x,yref);

NFFT = 2^nextpow2(length(x));

Fs = 1/abs(x(1)-x(2));                    % Sampling frequency

fftx = Fs/2*linspace(0,1,NFFT/2+1);
ffty = fft(y-mean(y),NFFT)/length(x);

ffty = ffty(1:NFFT/2+1);

figure(2); plot(fftx,abs(ffty));

[pks,locs] = max(abs(ffty));

Rabi = fftx(locs)/1e6 %in MHz
pipulse = 1e9/fftx(locs)/2 %in ns
piovertwopulse = pipulse/2 %in ns


