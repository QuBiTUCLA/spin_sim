function [] = abcd()

clear all

t = 50e-9:5e-9:2000e-9;

%static noise in Delta
Distr = 'Normal';
AD = 0; %A param for normal is mean
BD = 0.1; %B param for normal is std dev
stanoiseD = random(Distr,AD,BD,[length(t),1]);
noi = stanoiseD';

sig = 1*cos(2*pi*10*10^6*t) + 0.6*cos(2*pi*11*10^6*t) + noi; %+ 0.23;

%figure(1); plot(t,sig)

NFFT = 2^nextpow2(length(t));
Fs = 1/abs(t(1)-t(2));                    % Sampling frequency
fftx = Fs/2*linspace(0,1,NFFT/2+1);
ffty = fft(sig-mean(sig),NFFT)/length(t);
ffty = ffty(1:NFFT/2+1);
%%% PLOT FFT
%figure(1002); plot(fftx,abs(ffty));

Mmax = 2; %max number of presumed lines

Amax = sum(ffty);
xmin = fftx(1);
xmax = fftx(end);

w = 0.1*10^6;
W = 0.1*10^6;

%%%
proba = [];

for M = 1:1:Mmax
        
%%%%%%%%%%%%%%%%%%
fun = @(P,fftx) myfun(P,W,w,fftx,abs(ffty),BD,Mmax);     
pinit = [0.8, 10*10^6, 0.5, 11.2*10^6];
LB = [0,    9*10^6,     0, 10*10^6];
UB = [1,    11*10^6,    1, 13*10^6];
[pbest,delta_p]=easyfit(fftx, abs(ffty), pinit, fun, LB, UB);
%%%%%%%%%%%%%%%%%%%%%

%fun = @(P) myfun(P,W,w,fftx,ffty,BD,Mmax);     
    
%chisqmin  = fminsearch
%a = lsqnonlin(fun,[1, 10*10^6, 0.6, 11*10^6]);

end

l;l;l;l


%%%%%%%%%%%%%%%%%%%    


%deter =  111;



proba(M) = 1/((Amax*(xmax-xmin))^M)*exp(-result/2)*factorial(M)*(4*pi)^M/Sqrt(deter);


end

function chisq=myfun(P,W,w,fftx,ffty,BD,Mmax)

    chisq = 0;
    
    for k=1:1:length(fftx)
        
        y(k) = 0;

	for j = 1:1:Mmax
       
      y(k) = y(k) + P(2*j-1)*exp(((P(2*j) - fftx(k)))^2/(2*(W^2 + w^2)));  
        
    end 
    
    
    end
    
    chisq = chisq + ((y - abs(ffty))/(1/BD)).^2;
    chisq = chisq.';
   
end


