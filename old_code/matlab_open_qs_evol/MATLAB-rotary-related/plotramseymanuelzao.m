load('1DExp-seq-ramsey-vary-tau-2012-06-01-002910Manuelzao5MHz.mat');

rabi = Scan.ExperimentData{1}{1};
lowref = Scan.ExperimentData{1}{2};
oneref = Scan.ExperimentData{1}{3};
zeroref = Scan.ExperimentData{1}{4};
x = Scan.vary_begin:((Scan.vary_end -Scan.vary_begin)/(Scan.vary_points-1)):Scan.vary_end;

%rabi = (rabi - oneref)./(zeroref - oneref); 

%%%%% LAG PTS
lagpts = 10;
rabima = tsmovavg(rabi,'s',lagpts) ;
rabima = rabima(lagpts:1:end);
rabi = rabima;

%xlag = x(1:1:end-lagpts);
xlag = x(lagpts:1:end);
%xlag = x(floor(lagpts/2):1:end-floor(lagpts/2));
x = xlag;
%%%%%% LAG PTS

%%%%%% SHORTEN SIGNAL
rabi = rabi(1:1:350); %600 works well, tau = 1.9
x = x(1:1:350);
%%%%%% SHORTEN SIGNAL

init = [0,0,0,0,0,0];
Lo = [0,0,0,0,0,0];
Hi = [0,0,0,0,0,0];

ramseyfit = @(p,t) p(1) - (p(2)*cos(t*2*pi*1e6*p(3) + p(6)*pi) + p(2)*cos(t*2*pi*1e6*(p(3)+p(4))+ p(6)*pi) + p(2)*cos(t*2*pi*1e6*(p(3)-p(4)))+ p(6)*pi).*exp(-(t/p(5)/1e-6).^2);

%offset
init(1) = 12;
Lo(1) = 11.5; 
Hi(1) = 12.5; 

%amp 1
init(2) = 0.32;
Lo(2) = 0.2; 
Hi(2) = 0.55;

%freq delta
init(3) = 4.6;
Lo(3) = 4.3;
Hi(3) = 6;

%freq HF
init(4) = 2.17;
Lo(4) = 1.9;
Hi(4) = 2.3;

%T2*
init(5) = 2.4;
Lo(5) = 2;
Hi(5) = 4;

%phase phi -> phi is really close to 0
init(6) =0;
Lo(6) = -0.1;
Hi(6) =0.1;

% initial values 
pinit = init;

% bounds for fitting parameters 
LB = Lo;
UB = Hi;

[pbest,delta_p]=easyfit(x, rabi, pinit, ramseyfit, LB, UB);
pbest
delta_p
pbest(5)


figure(78)
plot(x,rabi,'ob',x,ramseyfit(pbest,x),'r')

%for 10 lagpts

 
 
