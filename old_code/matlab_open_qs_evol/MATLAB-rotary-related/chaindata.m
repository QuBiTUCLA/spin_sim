clear all

load('2to4-Rot2011-09-17-053055');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% All avgs
rot = Scan.ExperimentData{1}{1};
lowref = Scan.ExperimentData{1}{2};
oneref = Scan.ExperimentData{1}{3};
zeroref = Scan.ExperimentData{1}{4};
rabi = Scan.ExperimentData{1}{5};

x = Scan.vary_begin:((Scan.vary_end -Scan.vary_begin)/(Scan.vary_points-1)):Scan.vary_end;

rabinorm = (rabi - oneref)./(zeroref - oneref); %rabi normalized
rotnorm = (rot - oneref)./(zeroref - oneref); %rotary normalized

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Some Avgs

%%
% AUTO
%%

%array to discard
%this discarded array was obtained by function maximizepeak
%disc = [2:1:4 6 8 9 12:1:26 43:1:50]; %4to6-Rot2011-09-18-051509 
disc = [1:1:17 19 24 28 30 34 36:1:48 50]; %2to4-Rot2011-09-17-053055 
% [2:1:5 8 11 13 17 28 29 46]; %Rot2011-09-14-100503

aux2=1;
for aux=1:1:length(Scan.ExperimentDataEachAvg{1})
    if ~any(aux==disc)
rotCl{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{1};
oneCl{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{3};
zeroCl{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{4};
rabiCl{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{5};
aux2 = aux2+1;
    end
end

dim = ndims(rotCl{1});          %get the number of dimensions for your arrays

M1 = cat(dim+1,rotCl{:});        %convert to a (dim+1)-dimensional matrix
rot1 = mean(M1,dim+1);        %get the mean across arrays

M2 = cat(dim+1,oneCl{:});        %convert to a (dim+1)-dimensional matrix
one1 = mean(M2,dim+1);        %get the mean across arrays

M3 = cat(dim+1,zeroCl{:});        %convert to a (dim+1)-dimensional matrix
zero1 = mean(M3,dim+1);        %get the mean across arrays

M4 = cat(dim+1,rabiCl{:});        %convert to a (dim+1)-dimensional matrix
rabi1 = mean(M4,dim+1);        %get the mean across arrays

rotClean = (rot1 - one1)./(zero1 - one1); %rabi normalized
rabiClean = (rabi1 - one1)./(zero1 - one1); %rabi normalized

%%
% MANUAL
%%

%array to discard
%this discarded array was obtained manually
%array to discard
%this discarded array was obtained by looking
%disc = [1:1:3 9 5 10 11 16 21 22 32];
disc = [1:1:15 24 37 40 47]; %2to4-Rot2011-09-17-053055

aux2=1;
for aux=1:1:length(Scan.ExperimentDataEachAvg{1})
    if ~any(aux==disc)
rotClm{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{1};
oneClm{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{3};
zeroClm{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{4};
rabiClm{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{5};
aux2 = aux2+1;
    end
end

dim = ndims(rotClm{1});          %get the number of dimensions for your arrays

M1 = cat(dim+1,rotClm{:});        %convert to a (dim+1)-dimensional matrix
rot1m = mean(M1,dim+1);        %get the mean across arrays

M2 = cat(dim+1,oneClm{:});        %convert to a (dim+1)-dimensional matrix
one1m = mean(M2,dim+1);        %get the mean across arrays

M3 = cat(dim+1,zeroClm{:});        %convert to a (dim+1)-dimensional matrix
zero1m = mean(M3,dim+1);        %get the mean across arrays

M4 = cat(dim+1,rabiClm{:});        %convert to a (dim+1)-dimensional matrix
rabi1m = mean(M4,dim+1);        %get the mean across arrays

rotCleanm = (rot1m - one1m)./(zero1m - one1m); %rabi normalized
rabiCleanm = (rabi1m - one1m)./(zero1m - one1m); %rabi normalized

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Scan rot lowref oneref zeroref rabi;
load('4to6-Rot2011-09-18-051509');

%%%% All avgs
rot = Scan.ExperimentData{1}{1};
lowref = Scan.ExperimentData{1}{2};
oneref = Scan.ExperimentData{1}{3};
zeroref = Scan.ExperimentData{1}{4};
rabi = Scan.ExperimentData{1}{5};

x = [x Scan.vary_begin:((Scan.vary_end -Scan.vary_begin)/(Scan.vary_points-1)):Scan.vary_end];

rabinorm = [rabinorm (rabi - oneref)./(zeroref - oneref)]; %rabi normalized
rotnorm = [rotnorm (rot - oneref)./(zeroref - oneref)]; %rotary normalized

%%%% Disc Avgs AUTO
%array to discard
%this discarded array was obtained by function maximizepeak
disc = [2:1:4 6 8 9 12:1:26 43:1:50]; %4to6-Rot2011-09-18-051509 
%disc = [1:1:17 19 24 28 30 34 36:1:48 50]; %2to4-Rot2011-09-17-053055 
% [2:1:5 8 11 13 17 28 29 46]; %Rot2011-09-14-100503

aux2=1;
for aux=1:1:length(Scan.ExperimentDataEachAvg{1})
    if ~any(aux==disc)
rotCl{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{1};
oneCl{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{3};
zeroCl{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{4};
rabiCl{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{5};
aux2 = aux2+1;
    end
end

dim = ndims(rotCl{1});          %get the number of dimensions for your arrays

M1 = cat(dim+1,rotCl{:});        %convert to a (dim+1)-dimensional matrix
rot1 = mean(M1,dim+1);        %get the mean across arrays

M2 = cat(dim+1,oneCl{:});        %convert to a (dim+1)-dimensional matrix
one1 = mean(M2,dim+1);        %get the mean across arrays

M3 = cat(dim+1,zeroCl{:});        %convert to a (dim+1)-dimensional matrix
zero1 = mean(M3,dim+1);        %get the mean across arrays

M4 = cat(dim+1,rabiCl{:});        %convert to a (dim+1)-dimensional matrix
rabi1 = mean(M4,dim+1);        %get the mean across arrays

rotClean = [rotClean (rot1 - one1)./(zero1 - one1)]; %rabi normalized
rabiClean = [rabiClean (rabi1 - one1)./(zero1 - one1)]; %rabi normalized

%%%% Disc Avgs MANUALLY
clear rotClm oneClm zeroClm rabiClm

%array to discard
%this discarded array was obtained by looking
disc = [1:1:3 9 5 10 11 16 21 22 32];
%[1:1:15 24 37 40 47]; %2to4-Rot2011-09-17-053055

aux2=1;
for aux=1:1:length(Scan.ExperimentDataEachAvg{1})
    if ~any(aux==disc)
rotClm{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{1};
oneClm{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{3};
zeroClm{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{4};
rabiClm{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{5};
aux2 = aux2+1;
    end
end

dim = ndims(rotClm{1});          %get the number of dimensions for your arrays

M1 = cat(dim+1,rotClm{:});        %convert to a (dim+1)-dimensional matrix
rot1m = mean(M1,dim+1);        %get the mean across arrays

M2 = cat(dim+1,oneClm{:});        %convert to a (dim+1)-dimensional matrix
one1m = mean(M2,dim+1);        %get the mean across arrays

M3 = cat(dim+1,zeroClm{:});        %convert to a (dim+1)-dimensional matrix
zero1m = mean(M3,dim+1);        %get the mean across arrays

M4 = cat(dim+1,rabiClm{:});        %convert to a (dim+1)-dimensional matrix
rabi1m = mean(M4,dim+1);        %get the mean across arrays

rotCleanm = [rotCleanm (rot1m - one1m)./(zero1m - one1m)]; %rabi normalized
rabiCleanm = [rabiCleanm (rabi1m - one1m)./(zero1m - one1m)]; %rabi normalized

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Scan rot lowref oneref zeroref rabi;
clear rotClm oneClm zeroClm rabiClm;
clear rotCl oneCl zeroCl rabiCl;
load('0to2-Rot2011-09-18-220847');

%%%% All avgs
rot = Scan.ExperimentData{1}{1};
lowref = Scan.ExperimentData{1}{2};
oneref = Scan.ExperimentData{1}{3};
zeroref = Scan.ExperimentData{1}{4};
rabi = Scan.ExperimentData{1}{5};

x = [Scan.vary_begin:((Scan.vary_end -Scan.vary_begin)/(Scan.vary_points-1)):Scan.vary_end x];

rabinorm = [(rabi - oneref)./(zeroref - oneref) rabinorm]; %rabi normalized
rotnorm = [(rot - oneref)./(zeroref - oneref) rotnorm]; %rotary normalized

%%%% Disc Avgs AUTO
%array to discard
%this discarded array was obtained by function maximizepeak
disc = [1 3 5:1:7 9:1:11 16:1:19 21:1:29 31:1:34 37:1:39 41:1:43 45:1:48 50];  % Rot0to2 
%[2:1:4 6 8 9 12:1:26 43:1:50]; %4to6-Rot2011-09-18-051509 
%disc = [1:1:17 19 24 28 30 34 36:1:48 50]; %2to4-Rot2011-09-17-053055 
% [2:1:5 8 11 13 17 28 29 46]; %Rot2011-09-14-100503

aux2=1;
for aux=1:1:length(Scan.ExperimentDataEachAvg{1})
    if ~any(aux==disc)
rotCl{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{1};
oneCl{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{3};
zeroCl{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{4};
rabiCl{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{5};
aux2 = aux2+1;
    end
end

dim = ndims(rotCl{1});          %get the number of dimensions for your arrays

M1 = cat(dim+1,rotCl{:});        %convert to a (dim+1)-dimensional matrix
rot1 = mean(M1,dim+1);        %get the mean across arrays

M2 = cat(dim+1,oneCl{:});        %convert to a (dim+1)-dimensional matrix
one1 = mean(M2,dim+1);        %get the mean across arrays

M3 = cat(dim+1,zeroCl{:});        %convert to a (dim+1)-dimensional matrix
zero1 = mean(M3,dim+1);        %get the mean across arrays

M4 = cat(dim+1,rabiCl{:});        %convert to a (dim+1)-dimensional matrix
rabi1 = mean(M4,dim+1);        %get the mean across arrays

rotClean = [(rot1 - one1)./(zero1 - one1) rotClean ]; %rabi normalized
rabiClean = [(rabi1 - one1)./(zero1 - one1) rabiClean ]; %rabi normalized

%%%% Disc Avgs MANUALLY
clear rotClm oneClm zeroClm rabiClm

%array to discard
%this discarded array was obtained by looking
disc = [1]; %0to6
%[1:1:3 9 5 10 11 16 21 22 32]; %Rot4to6
%[1:1:15 24 37 40 47]; %2to4-Rot2011-09-17-053055

aux2=1;
for aux=1:1:length(Scan.ExperimentDataEachAvg{1})
    if ~any(aux==disc)
rotClm{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{1};
oneClm{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{3};
zeroClm{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{4};
rabiClm{aux2} = Scan.ExperimentDataEachAvg{1}{aux}{5};
aux2 = aux2+1;
    end
end

dim = ndims(rotClm{1});          %get the number of dimensions for your arrays

M1 = cat(dim+1,rotClm{:});        %convert to a (dim+1)-dimensional matrix
rot1m = mean(M1,dim+1);        %get the mean across arrays

M2 = cat(dim+1,oneClm{:});        %convert to a (dim+1)-dimensional matrix
one1m = mean(M2,dim+1);        %get the mean across arrays

M3 = cat(dim+1,zeroClm{:});        %convert to a (dim+1)-dimensional matrix
zero1m = mean(M3,dim+1);        %get the mean across arrays

M4 = cat(dim+1,rabiClm{:});        %convert to a (dim+1)-dimensional matrix
rabi1m = mean(M4,dim+1);        %get the mean across arrays

rotCleanm = [(rot1m - one1m)./(zero1m - one1m) rotCleanm]; %rabi normalized
rabiCleanm = [(rabi1m - one1m)./(zero1m - one1m) rabiCleanm ]; %rabi normalized



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT

% FFT all avgs
%%%% FFT Rabi
NFFT = 2^nextpow2(length(x));
Fs = 1/abs(x(1)-x(2));                    % Sampling frequency
fftx = Fs/2*linspace(0,1,NFFT/2+1);
ffty = fft(rabinorm-mean(rabinorm),NFFT)/length(x);
ffty = ffty(1:NFFT/2+1);

%%%% FFT Rot
NFFT = 2^nextpow2(length(x));
Fs = 1/abs(x(1)-x(2));                    % Sampling frequency
fftx2 = Fs/2*linspace(0,1,NFFT/2+1);
ffty2 = fft(rotnorm-mean(rotnorm),NFFT)/length(x);
ffty2 = ffty2(1:NFFT/2+1);

% FFT AUTO Clean
%%%% FFT Rabi
NFFT = 2^nextpow2(length(x));
Fs = 1/abs(x(1)-x(2));                    % Sampling frequency
fftxc = Fs/2*linspace(0,1,NFFT/2+1);
fftyc = fft(rabiClean-mean(rabiClean),NFFT)/length(x);
fftyc = fftyc(1:NFFT/2+1);

%%%% FFT Rot
NFFT = 2^nextpow2(length(x));
Fs = 1/abs(x(1)-x(2));                    % Sampling frequency
fftx2c = Fs/2*linspace(0,1,NFFT/2+1);
ffty2c = fft(rotClean-mean(rotClean),NFFT)/length(x);
ffty2c = ffty2c(1:NFFT/2+1);

% FFT MANUAL Clean
%%%% FFT Rabi
NFFT = 2^nextpow2(length(x));
Fs = 1/abs(x(1)-x(2));                    % Sampling frequency
fftxcm = Fs/2*linspace(0,1,NFFT/2+1);
fftycm = fft(rabiCleanm-mean(rabiCleanm),NFFT)/length(x);
fftycm = fftycm(1:NFFT/2+1);
%%% PLOT FFT
%figure(50020); plot(fftxc,abs(fftyc));
%[pks,locs] = max(abs(ffty));

%%%% FFT Rot
NFFT = 2^nextpow2(length(x));
Fs = 1/abs(x(1)-x(2));                    % Sampling frequency
fftx2cm = Fs/2*linspace(0,1,NFFT/2+1);
ffty2cm = fft(rotCleanm-mean(rotCleanm),NFFT)/length(x);
ffty2cm = ffty2cm(1:NFFT/2+1);

%%%% PLOT RABI
figure(99);
subplot(3,2,1);
plot(x,rabinorm,'r')
title('rabi all avgs')
axis([min(x) max(x) -0.5 1.5])
subplot(3,2,2);
plot(fftx,abs(ffty),'r')
axis([min(fftx) max(fftx) 0 0.08])
subplot(3,2,3);
plot(x,rabiClean,'r')
title('rabi cleaned auto')
axis([min(x) max(x) -0.5 1.5])
subplot(3,2,4);
plot(fftxc,abs(fftyc),'r')
axis([min(fftxc) max(fftxc) 0 0.08])
subplot(3,2,5);
plot(x,rabiCleanm,'r')
title('rabi cleaned manually')
axis([min(x) max(x) -0.5 1.5])
subplot(3,2,6);
plot(fftxcm,abs(fftycm),'r')
axis([min(fftxcm) max(fftxcm) 0 0.08])

%%%% PLOT RABI
figure(100);
subplot(3,2,1);
plot(x,rotnorm,'b')
title('rot all avgs')
axis([min(x) max(x) -0.5 1.5])
subplot(3,2,2);
plot(fftx2,abs(ffty2),'b')
axis([min(fftx2) max(fftx2) 0 0.04])
subplot(3,2,3);
plot(x,rotClean,'b')
title('rot cleaned auto')
axis([min(x) max(x) -0.5 1.5])
subplot(3,2,4);
plot(fftx2c,abs(ffty2c),'b')
axis([min(fftx2c) max(fftx2c) 0 0.04])
subplot(3,2,5);
plot(x,rotCleanm,'b')
title('rot cleaned manually')
axis([min(x) max(x) -0.5 1.5])
subplot(3,2,6);
plot(fftx2cm,abs(ffty2cm),'b')
axis([min(fftx2cm) max(fftx2cm) 0 0.04])

