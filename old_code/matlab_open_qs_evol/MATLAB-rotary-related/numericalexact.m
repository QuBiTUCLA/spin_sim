function numericalexact()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Rabi g1 g2 anglerot Delta DabEq aux nb_aux noise;

do_rabi =0;
do_rotary = 1;

Rabi = 2*pi*20e6;        %10e6; % bare Rabi in circular freq
g1 = 0;%1/1e-3; %1/100e-3; %Right now T1 infinite
g2 = 0; %1/1e-6;

%For rotary
anglerot =pi;
Delta=2*pi*1e6*(1)*2.1;   %2.1;     %2.1e6;        %3/100*Rabi; % By definition, w0 - w, but here parametrized by % of Rabi;

DabEq =0; %all pop in ground state
rhoggEq = (-DabEq + 1)/2;

u_init = 0;
v_init = 0;
w_init = -1;

nbcycle =30;    % ceil(Rabi/anglerot/g2); 
%nbcycle is Rotary cycles
% total time: nbcycle*2*anglerot/Rabi
timeunit =anglerot/Rabi;
period = 2*timeunit;

N = 20; % nb of points per cycle/2  ie nb of points per cycle is 2*N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NUMERICAL

%original
%options = odeset('RelTol',1e-12,'AbsTol',ones(1,3)*1e-12);
%modified to see if improves convergence
options = odeset('RelTol',1e-3,'AbsTol',ones(1,3)*1e-6);
%used (engouh) -3,-6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Numerical no noise

Wnn = [];
Unn = [];
Vnn = [];

WRotnn = [];
URotnn = [];
VRotnn = [];

u_new = u_init;
v_new = v_init;
w_new = w_init;


%% Usual Rabi
if do_rabi
for q = 0:1:nbcycle-1
    
[T,Y] = ode45(@Propagate_phase0nn,[q*period (2*q+1)/2*period],[u_new v_new w_new],options);

t0 = q*period:period/2/N:(2*q+1)/2*period;
y0=interp1(T,Y,t0);

if q ~= 0
   Wnn = Wnn(1:1:end-1);
   Unn = Unn(1:1:end-1);
   Vnn = Vnn(1:1:end-1);
end
Wnn = [Wnn y0(:,3).'];
Unn = [Unn y0(:,1).'];
Vnn = [Vnn y0(:,2).'];

u_new = y0(end,1);
v_new = y0(end,2);
w_new = y0(end,3);

[TT,YY] = ode45(@Propagate_phase0nn,[(2*q+1)/2*period (q+1)*period],[u_new v_new w_new],options);

tt0 = (2*q+1)/2*period:period/2/N:(q+1)*period;

yy0=interp1(TT,YY,tt0);

Wnn = Wnn(1:1:end-1);
Unn = Unn(1:1:end-1);
Vnn = Vnn(1:1:end-1);

Wnn = [Wnn yy0(:,3).'];
Unn = [Unn yy0(:,1).'];
Vnn = [Vnn yy0(:,2).'];

u_new = yy0(end,1);
v_new = yy0(end,2);
w_new = yy0(end,3);

end
end

%% Rotary
if do_rotary
u_new = 0;
v_new = 0;
w_new = -1;
for q = 0:1:nbcycle-1
    
[T,Y] = ode45(@Propagate_phase0nn,[q*period (2*q+1)/2*period],[u_new v_new w_new],options);

t0 = q*period:period/2/N:(2*q+1)/2*period;
y0=interp1(T,Y,t0);

if q ~= 0
   WRotnn = WRotnn(1:1:end-1);
   URotnn = URotnn(1:1:end-1);
   VRotnn = VRotnn(1:1:end-1);
end
WRotnn = [WRotnn y0(:,3).'];
URotnn = [URotnn y0(:,1).'];
VRotnn = [VRotnn y0(:,2).'];

u_new = y0(end,1);
v_new = y0(end,2);
w_new = y0(end,3);

[TT,YY] = ode45(@Propagate_phase180nn,[(2*q+1)/2*period (q+1)*period],[u_new v_new w_new],options);

tt0 = (2*q+1)/2*period:period/2/N:(q+1)*period;

yy0=interp1(TT,YY,tt0);

WRotnn = WRotnn(1:1:end-1);
URotnn = URotnn(1:1:end-1);
VRotnn = VRotnn(1:1:end-1);

WRotnn = [WRotnn yy0(:,3).'];
URotnn = [URotnn yy0(:,1).'];
VRotnn = [VRotnn yy0(:,2).'];

u_new = yy0(end,1);
v_new = yy0(end,2);
w_new = yy0(end,3);
 
end
end

%FITS
% %PRB modified results without noise for rotary
 n=0:1/2/N:nbcycle;
 c = (3*g1 + 2*g2)/4;
C = sqrt(1 + c^2/Rabi^2);
 theta = atan(Rabi/c);
% 
 prb = 0.5*(sin(anglerot/2)^2*C*sin(2*n*pi + theta).*cos(4*Delta*n/Rabi*sin(anglerot/2)).*exp(-(c*2*n*anglerot/Rabi)*(3/2)*pi/anglerot*sin(anglerot/2))) + 0.5*(1 + cos(anglerot/2)^2*exp(-(c*2*n*anglerot/Rabi)*anglerot/pi*sin(anglerot/2)));

%%%% No relax
%norelaxrot = 0.5*(1 + cos(anglerot/2)^2) + 0.5*sin(anglerot/2)^2.*cos(4*n*Delta/Rabi*sin(anglerot/2)).*cos((2*pi*n/Rabi)*Rabi);
norelaxrot =  0.5*(1 + cos(anglerot/2)^2) + 0.5*sin(anglerot/2)^2.*cos((2*pi*n/Rabi)*(2*Delta/pi)*sin(anglerot/2)).*cos((2*pi*n/Rabi)*Rabi);

%%%%% in time
tunit =  (2*anglerot*nbcycle/Rabi)/(nbcycle*N*2); %1e-9; %make same number per cycle as in the simulation
tfinal = 2*anglerot*nbcycle/Rabi;
t = 0:tunit:tfinal;
%norelaxrot2 =  0.5*(1 + cos(anglerot/2)^2) + 0.5*sin(anglerot/2)^2.*cos((t/anglerot)*(2*Delta)*sin(anglerot/2)).*cos((t*pi/anglerot)*Rabi);

%%correct, without mod, valid for theta < 2*pi
norelaxrot3 =  0.5*(1 + cos(anglerot/2)^2) + 0.5*sin(anglerot/2)^2.*cos((t*pi/anglerot)*(2*Delta/pi)*sin(anglerot/2)).*cos((t*pi/mod(anglerot,2*pi))*Rabi); %correct, with mod
relaxrot3 = 0.5*C*sin(anglerot/2)^2.*cos((t*pi/anglerot)*(2*Delta/pi)*sin(anglerot/2)).*sin((t*pi/mod(anglerot,2*pi))*Rabi + theta).*exp(-c*t*3*pi/2/anglerot*sin(anglerot/2)) + 0.5*(1 + cos(anglerot/2)^2.*exp(-c*t*anglerot/pi*sin(anglerot/2))); %correct
norelaxrot44 = 1 - sin(anglerot/2)^2*sin(2*Delta*(t*Rabi/2/anglerot)/Rabi*sin(anglerot/2)).^2.*(sin((t*pi/mod(anglerot,2*pi))*Rabi));



%Paola's round formula
Proundrot =  0.5*(1 + cos(anglerot/2)^2) + 0.5*sin(anglerot/2)^2.*cos(4*Delta/Rabi*round(Rabi*t/2/anglerot)).*cos(t*Rabi); %ORIGINAL

roundrot =  0.5*(1 + cos(anglerot/2)^2) + 0.5*sin(anglerot/2)^2.*(cos(4*Delta/Rabi*floor(Rabi*t/2/anglerot + 1/2)).*cos(t*Rabi/2).^2); %WORKIGN

%addi =  0.5*sin(anglerot/2)^2.*(-cos(4*Delta/Rabi*floor(Rabi*t/2/anglerot + 32)).*sin(t*Rabi/2).^2) + 0.5; %works for Rabi/Delta = 20;
addi =  0.5*sin(anglerot/2)^2.*(-cos(4*Delta/Rabi*floor(Rabi*t/2/anglerot) + 4*Delta/Rabi*0.5).*sin(t*Rabi/2).^2) + 0.5; %ALSO WORKS for Rabi/Delta = 20;
%addi =  0.5*sin(anglerot/2)^2.*(cos(4*Delta/Rabi*floor(-Rabi*t/2/anglerot + 16)).*sin(t*Rabi/2).^2) + 0.5; %TRY NOT WORKING Rabi/Delta = 20;
%addi =  0.5*sin(anglerot/2)^2.*(cos(4*Delta/Rabi*floor(-Rabi*t/2/anglerot + 24)).*sin(t*Rabi/2).^2) + 0.5; %works for Rabi/Delta = 30;
%%%%%%%%%%addi =  0.5*sin(anglerot/2)^2.*(cos(4*Delta/Rabi*floor(Rabi*t/2/anglerot - 31)).*sin(t*Rabi/2).^2) + 0.5; %works for Rabi/Delta = 40;
%addi =  0.5*sin(anglerot/2)^2.*(cos(4*Delta/Rabi*floor(-Rabi*t/2/anglerot + 32)).*sin(t*Rabi/2).^2) + 0.5; %NEW TRY Rabi/Delta = 40;
%addi =  0.5*sin(anglerot/2)^2.*(cos(4*Delta/Rabi*floor(-Rabi*t/2/anglerot +40)).*sin(t*Rabi/2).^2) + 0.5; %works for Rabi/Delta = 50;


norelaxrottest =  0.5*(1 + cos(anglerot/2)^2) + 0.5*sin(anglerot/2)^2.*cos((t*pi/anglerot)*(2*Delta/pi)*sin(anglerot/2)).*cos((t*pi/mod(anglerot,2*pi))*Rabi); %TEST

rabitest =(1 - Rabi^2/(Rabi^2+Delta^2)*sin(t/2*sqrt(Delta^2+Rabi^2)).^2);
%For anglerot > 2*pi, both with and without mod pick up the echo peaks but
%are off for the other parts of the echo: the with mod option oscillates
%quicker, the without mod option oscillates slower

%%%%% Test of frequencies only, forget about normalization
% noddc -> DOUBLE PEAKS
%auxodd = zeros(1,length(t));    %cos(t*Rabi*pi/mod(anglerot,2*pi));
auxodd = cos(t*Rabi*pi/mod(anglerot,2*pi));
nodd_array = 1:2:31;
numro = 0;
for nodd = nodd_array
   %if (nodd - int64(anglerot/pi) ~= 0)
   %auxodd = auxodd + 1*(Delta/Rabi/nodd)^2*cos(nodd*t*Rabi*pi/anglerot);
   %%so far the best fit is above
   %auxodd = auxodd + 1*(Delta/Rabi)/nodd^2*cos(nodd*t*Rabi*pi/anglerot);
   %above not too bad a fit either
   auxodd = auxodd + (1/nodd)*cos(nodd*t*Rabi*pi/anglerot);
   numro = numro + 1/nodd;
   %end   
end
%auxodd = auxodd/(1 + (1*Delta/Rabi)*numro^2);
auxodd = auxodd/(1 + numro);

%neven -> SINGLE PEAKS
auxeven  = zeros(1,length(auxodd));
neven_array = 2:2:12;
numm = 0;
for neven = neven_array
    %auxeven = auxeven + (Delta/Rabi/neven)*cos(neven*t*Rabi*pi/anglerot); 
    auxeven = auxeven + (1/neven^2)*cos(neven*t*Rabi*pi/anglerot); 
    numm = numm + 1/neven;
 end
%auxeven = auxeven/(Delta/Rabi*numm);
%auxeven = auxeven/(numm);

% auxodd = cos(t*Rabi*pi/mod(anglerot,2*pi));
% n_array = 1:1:10;
% no = 0;
% for indn = n_array
%     auxodd = auxodd + 1*(Delta/Rabi/indn)*cos(indn*t*Rabi*pi/anglerot);
%     no = no + 1/indn;
% end
% auxodd = auxodd/(1 + 1*Delta/Rabi*no);

%signal = 0.5*(1 + cos(anglerot/2)^2) + 0.5*sin(anglerot/2)^2.*(cos((t/anglerot)*(2*Delta)*sin(anglerot/2)).*auxodd); % + Delta/Rabi*auxeven);
signal = 0.5*(1 + cos(anglerot/2)^2) + 0.5*sin(anglerot/2)^2.*(cos(4*Delta/Rabi*round(Rabi*t/2/anglerot)).*auxodd); % + Delta/Rabi*auxeven);
%signal = cos((t/anglerot)*(2*Delta)*sin(anglerot/2)).*auxodd + ((1/20)*auxeven); 
% cut below to test
%
%up to auxodd array, formula above is same formula as in norelaxrot3
%factor 2 inside modulation is not correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting

figure(100);
hold on
if do_rabi
%%plot(0:1/20:nbcycle,Zavg,'r');  %numerical with noise
%plot([0:1/2/N:nbcycle]*2*anglerot/Rabi,(1-Wnn)/2,'r'); %numerical no noise
plot([0:1/2/N:nbcycle],(1-Wnn)/2,'k'); %numerical no noise
%plot(t,rabitest,'m');
hold on
end

if do_rotary
%plot(0:1/20:nbcycle,Zrotavg,'b');


%plot([0:1/2/N:nbcycle],(1-WRotnn)/2,'r'); %USED plot in time
% 
% %working more or less
% %triangle = (cos(Delta*t*(4*anglerot*Delta/Rabi))).*(cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))))/2 + 0.5; 
% 
% triangle = (cos(anglerot*Delta*t/5) ).*(cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))))/2 + 0.5 ; %working more or less
% 
% %triangle2 = (cos(anglerot*Delta/5*(t+5*((floor(-sin(2*pi*1/5*t)))+1)) )).*(cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))))/2 + 0.5 ;
% %triangle2 = anglerot*Delta/5*(1*t+2*4*Rabi/5/Delta*((floor(-sin(4*Rabi*t/5/Delta)))+1));
% 
% %triangle2 = (1 * t*Rabi/2/anglerot + 4*Rabi/5/Delta * ((floor(-sin(t * Rabi/2/anglerot)))+1));
 %time_shifting_func = @(t, a) -2*a*( (t)/a  - floor( (t)/a + 1/2 ))+(t);
% 
% %triangle2 = (0.5 * cos(anglerot*4*Delta/Rabi*time_shifting_func(t*Rabi/2/an%glerot,4*Rabi/Delta/5)/5) + 0.5 );%.*(cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))))/2 + 0.5 ; %working more or less
%  
% %triangle2 = (0.5 * cos(anglerot * Delta/5*time_shifting_func(t*Rabi/2/anglerot,4*Rabi/Delta/5)) + 0.5 );%.*(cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))))/2 + 0.5 ; %working more or less
%  
%triangle2 = (cos(anglerot * Delta *t*(2*pi*0.0322))).*(cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))))/2 + 0.5 ; %working more or less

%triangle2 = 0.5*sin(anglerot/2)^2*(cos(t*2*Delta*sin(anglerot/2)/anglerot)).*(cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5)))) + 0.5 + 0.5*cos(anglerot/2)^2; %working more or less

%triangle4 = 0.5*sin(anglerot/2)^2*cos(Delta*t*(pi^2*2*0.0322)/sin(anglerot/2)).*(cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5)))) +0.5 +0.5*cos(anglerot/2)^2;

%%%%%%%%%%%%%working last one
%triangle4 =0.5*sin(anglerot/2)^2*cos(Delta*t*0.6356/sin(anglerot/2)).*(cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5)))) +0.5 + 0.5*cos(anglerot/2)^2; 

%triangle4 =0.5*(cos(anglerot * Delta *t*(2*pi*0.0322))).*cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))) +0.5; 
%triangle5 =0.5*(cos(anglerot * Delta *t*(2*pi*0.0322)+pi)).*cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot)).*cos(Delta*t/3)) +0.5; 
%triangle4 =0.5*cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5)).*cos(Delta*t/3)) +0.5; 


%used%%%triangle4 =0.5*cos(2*anglerot*abs(t*Delta/anglerot^2 - floor(t*Delta/anglerot^2 +0.5))).*cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))) +0.5; 
%resi = 1.5*Delta/Rabi/anglerot.*sin(2*anglerot*abs(t*Delta/anglerot^2 - floor(-pi/4+ t*Delta/anglerot^2 +0.5))).*cos(2*anglerot*abs(pi/4+ t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))); 


triangle4 =0.5*cos((2*pi*0.0322)*t*Delta*anglerot).*cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))) +0.5; 
resi =Delta/Rabi/2*(1).*sin((2*pi*0.0322)*t*Delta*anglerot).*cos(2*anglerot*abs(pi/4 + t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))); 



%resi = 1.5*Delta/Rabi/anglerot*sin(2*anglerot*abs(t*Delta/anglerot^2 - floor(-pi/4 + t*Rabi/2/anglerot +0.5))).*cos(2*anglerot*abs(pi/4 + t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))); 
%resi = 2*Delta/Rabi/anglerot*sin(anglerot * Delta *t*(2*pi*0.0322)).*cos(2*anglerot*abs(pi/4 + t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))); 

%triangle4 =0.5*(cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5)))) +0.5;   

%triangle4 =0.5*sin(anglerot/2)^2*cos(Delta*0.6356/sin(anglerot/2)*time_shifting_func(t,5*Delta*t/8/anglerot)).*(cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5)))) +0.5 + 0.5*cos(anglerot/2)^2;

%triangle4 =0.5*cos(Delta*t*0.6356/sin(anglerot/2)).*(cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))))+0.5;

%formNoRWA = 0.5 - 0.5*cos(abs(2*(Rabi^2*t/Delta/anglerot - floor(Rabi^2*t/Delta/anglerot + 0.5)))); 
%formNoRWA = 0.5 - 0.5*cos(abs(2*(t/anglerot*Delta - floor(t/anglerot*Delta + 0.5)))); 
%formNoRWA = 0.5 + 0.5*cos(2*anglerot/Delta*(1 - cos(Delta/2/anglerot*t)));

%formNoRWA = 0.5 - 0.5*cos(0 + anglerot/Rabi*abs(2*(t*Rabi^2/Delta/2/anglerot - floor(t*Rabi^2/Delta/2/anglerot + 0.5))));
%formNoRWA = 0.5 - 0.5*cos(0 + anglerot*abs(2*(t/Delta*Rabi^2/2/anglerot - floor(t/Delta*Rabi^2/2/anglerot + 0.5))));

% 
% %triangle3 =  (-0.5 * cos(anglerot * Delta *t/5 ) + 0.5);
% 
% % timeincycle = t/(2*anglerot/Rabi);
% % for aux = 0:2:(timeincycle(end)/8)
% % [min_diff, arrayposmin] = min(abs(timeincycle-Rabi/Delta/5*2*(aux-1)));
% % [min_diff, arrayposmax] = min(abs(timeincycle-Rabi/Delta/5*2*(aux+1)));
% % triangle2(arrayposmin:1:arrayposmax) = cos(anglerot*Delta/5*(4*Rabi/5/Delta*aux) - timeincycle(arrayposmin:1:arrayposmax));     
% % %triangle2 = (0.5 * cos(1/5*anglerot * Delta *(-t + 4*Rabi/5/Delta*(floor(t/(2*Rabi/anglerot)) + 1/2))) + 0.5)
% % end
%


%used%newf = 0.5*sin(anglerot/2)^2*cos(2*anglerot/sin(anglerot/2)^2*abs(t*Delta/anglerot/pi*sin(anglerot/2) - floor(t*Delta/anglerot/pi*sin(anglerot/2)))).*cos(2*anglerot/sin(anglerot/2)*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))) +0.5 + 0.5*cos(anglerot/2)^2;
%newf = 0.5*sin(anglerot/2)^2*cos(0.6356*t*Delta).*(cos(2*anglerot/sin(anglerot/2)^1*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))))+0.5 +0.5*cos(anglerot/2)^2;

%%%%%%%%%%%%%%%%%%%%newf = 0.5*cos(2*t*Delta/pi -0.5*t*Delta/pi*Delta^2/Rabi^2-0.5*t*Delta/pi*Delta^3/Rabi^3).*(cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))))+0.5;
newf = 0.5*cos(2*t*Delta/pi -0.5*t*Delta/pi*Delta^2/Rabi^2-0.5*t*Delta/pi*Delta^3/Rabi^3).*(cos((sqrt(Rabi^2+Delta^2))/Rabi*2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))))+0.5;

%newf = 0.5*cos(2*pi*t*Delta/pi^2 -0.5*pi*t*Delta/pi^2*(1/(1-Delta/Rabi)-1-Delta/Rabi)).*(cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))))+0.5;
%newf = 0.5*cos(2*t*Delta/pi*sqrt(1 - Delta^2/2/Rabi^2)).*(cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))))+0.5;
%resi = Delta/Rabi/2*floor(t.*Rabi/2/anglerot/16 + 1).*cos(-pi/2 + 2*t*Delta/anglerot*sin(anglerot/2)).*(cos(2*anglerot/sin(anglerot/2)^0*abs(pi/4 + t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))));

%used%%%resi = tempo.*(Delta/Rabi/2).*sin(2*pi*abs(t*Delta./pi^2 - floor(t*Delta./pi^2 + 0.5))).*(cos(2*anglerot*abs(pi/4 + t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))));
%resi = (Delta/Rabi/2).*sin(2*pi*t*Delta/pi^2  -0.5*pi*t*Delta/pi^2*Delta^2/Rabi^2-0.5*pi*t*Delta/pi^2*Delta^3/Rabi^3);%.*(cos(2*anglerot*abs(pi/4 + t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))));
resi = (Delta/Rabi/2).*sin(2*t*Delta/pi -0.5*t*Delta/pi*Delta^2/Rabi^2-0.5*t*Delta/pi*Delta^3/Rabi^3).*cos(sqrt(Rabi^2+Delta^2)/Rabi*(2*anglerot*abs(pi/4 + t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))));% - (Delta^2/Rabi^2/4);%

%resi2 = Delta^2/Rabi^2/2*cos(-pi/2 + 2*t*Delta/anglerot*sin(anglerot/2)).*(cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))));
%del = Delta/Rabi/2*ones(1,length(newf));
%del2 = Delta/Rabi/2 + Delta^2/Rabi^2/2*t.*Rabi/2/anglerot/16;
%newf = 0.5*sin(anglerot/2)^2*(besselj(1,2*t*Delta/anglerot)).*cos(2*anglerot/sin(anglerot/2)*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))) +0.5 + 0.5*cos(anglerot/2)^2;

omega = sqrt(Rabi^2 + Delta^2);
triangle  = 2*anglerot/Rabi*abs(t./(2*anglerot/Rabi) - floor(t./(2*anglerot/Rabi) + 0.5));
arg = omega*triangle/2;
sq = square(pi*Rabi/anglerot*t); %both square functions a bit delayed, i think its only numerics
%sq = ones(1,length(t));
%for aux=1:1:length(t)
 %  if mod(floor(t(aux)*Rabi/anglerot),2) == 1
  %    sq(aux) = -1; 
  % end
%end
%mathefun = cos(arg).^2 + Delta^2*sin(arg).^2/omega^2;
%mathefun = sin(arg).^2./((sq).^2);

%exa = 0.5 + 0.5*cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5)));

%figure(33)
%plot([0:1/2/N:nbcycle],mathefun,'r',[0:1/2/N:nbcycle],sq,'b')
%ll;l
% close all
%plot([0:1/2/N:nbcycle],triangle4,'g-') 



%plot([0:1/2/N:nbcycle],norelaxrottest,'r')



tri = 2*anglerot/Rabi*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5));
michaels = cos(tri*Rabi/2).^2 - Delta*Rabi/2*sin(tri*Rabi) + Delta^2*Rabi^2/4*sin(tri*Rabi/2).^2;
%plot([0:1/2/N:nbcycle],michaels,'r')

plot([0:1/2/N:nbcycle],(1-WRotnn)/2,'b') %,[0:1/2/N:nbcycle],newf + 1*resi,'r',[0:1/2/N:nbcycle], (1-WRotnn)/2 - newf,'b',[0:1/2/N:nbcycle],resi,'g',[0:1/2/N:nbcycle],(1-WRotnn)/2 - newf - resi-0.05,'r') %,[0:1/2/N:nbcycle],resi2-0.1,'r');


%%%%Analytical formula AHT rotated back
Om = sqrt(4*Delta^2*sin(anglerot/2)^2 + anglerot/2*Rabi^2);
aht = cos(t*Om/2/anglerot).^2 + 4*Delta^2/Om^2*cos(anglerot/2 + Rabi*tri).^2*sin(anglerot/2)^2.*sin(t*Om/2/anglerot).^2;
%%%%
hold on;
%plot([0:1/2/N:nbcycle],aht,'r')
plot([0:1/2/N:nbcycle],norelaxrot3,'r')
hold on;
plot([0:1/2/N:nbcycle],norelaxrot44,'g')

l;l;l
%figure(17)
%plot([0:1/2/N:nbcycle], (1-WRotnn)/2 - newf,'b',[0:1/2/N:nbcycle], (1-WRotnn)/2 - newf2,'r')

l;ll
plot([0:1/2/N:nbcycle],(1-WRotnn)/2,'b',[0:1/2/N:nbcycle],triangle4,'g',[0:1/2/N:nbcycle], (1-WRotnn)/2 - triangle4,'r',[0:1/2/N:nbcycle],resi,'b') %,[0:1/2/N:nbcycle],norelaxrot3, 'k'); %USED plot in time
 axis([0 nbcycle -0.2 1])
 
figure(56)
plot([0:1/2/N:nbcycle],(1-WRotnn)/2,'b',[0:1/2/N:nbcycle],triangle4+resi,'g',[0:1/2/N:nbcycle], (1-WRotnn)/2 - triangle4-resi,'r') %,[0:1/2/N:nbcycle],resi,'b') %,[0:1/2/N:nbcycle],norelaxrot3, 'k'); %USED plot in time
 axis([0 nbcycle -0.2 1])
 
 
% 
% klklklkkl
% 
% close all
% plot([0:1/2/N:nbcycle],(1-WRotnn)/2 ,'r',[0:1/2/N:nbcycle],triangle2,'b'); %USED plot in time
% axis([0 nbcycle 0 1])
% klklklk
% 
%  
% plot([0:1/2/N:nbcycle],(1-WRotnn)/2 ,'r',[0:1/2/N:nbcycle],0.01*time_shifting_func(t*Rabi/2/anglerot, 4*Rabi/Delta/5),'g*'); %USED plot in time
% axis([0 nbcycle 0 1])
% klklklk
% 
% plot([0:1/2/N:nbcycle],(1-WRotnn)/2 - triangle,'g')
% 
% resgreen = 0.1*(cos(Delta*t*(4*anglerot*Delta/Rabi) + pi/2)).*(cos(2*anglerot*abs(t*Rabi/2/anglerot - floor(t*Rabi/2/anglerot + 0.5))))/2; %+ 0.5;
% plot([0:1/2/N:nbcycle]*2*anglerot/Rabi,resgreen,'m')
% 
% %resdores = (1-WRotnn)/2 - triangle - resgreen - 0.1; %-0.1 just to plot below
% %plot([0:1/2/N:nbcycle]*2*anglerot/Rabi,resdores,'k')

%to save rotary data
% tempo3 = [0:1/2/N:nbcycle]*2*anglerot/Rabi;
% signal3= (1-WRotnn)/2;
% save('3.mat','tempo3','signal3')

m,m,,m,m,

%plot([0:1/20:nbcycle],(1-WRotnn)/2,'b'); %USED plot in cycles
%models
%plot([0:1/20:nbcycle]*2*anglerot/Rabi,norelaxrot,'+k');

%plot(t,norelaxrot3,'g+');
%plot(t,relaxrot3,'g+');


%plot(t,Proundrot,'g')
%plot(t,roundrot + addi - 0.5,'b*')
%plot(t,signal,'g')

diffe = (1-WRotnn)/2 - roundrot + 0.5; % difference num and roundrot
plot(t/(2*anglerot/Rabi),diffe, 'r')
plot(t/(2*anglerot/Rabi),addi,'k*')


plot(t/(2*anglerot/Rabi),(diffe-addi),'g')
%trying to add that slow osci fiven by diffe-addi
%newosci = (Rabi.*t/2/anglerot)*(Delta/500/Rabi).*cos(Rabi*t).*sin(4*Delta/Rabi*floor(Rabi*t/2/anglerot + 0.5))    ;
%newosci = (t*Delta^2/8/anglerot/Rabi/2/pi).*cos(Rabi*t).*sin(4*Delta/Rabi*floor(Rabi*t/2/anglerot + 0.5))    ;
%newosci = 0.01*cos(Rabi*t*pi/anglerot).*sin(2*Delta*t/anglerot*sin(anglerot/2)); %miniecho
newosci = (Delta/Rabi).*(Rabi*t/2/anglerot/505).*(sin(4*Delta/Rabi*floor(Rabi*t/2/anglerot+ 0.5)).*cos(t*Rabi/2).^2 - sin(4*Delta/Rabi*(floor(Rabi*t/2/anglerot) + 0.5)).*sin(t*Rabi/2).^2);
plot(t/(2*anglerot/Rabi),newosci,'m*')

%plot([0:1/20:nbcycle]*2*anglerot/Rabi,prb,'k'); %USED plot in time
%plot(t/(2*anglerot/Rabi),norelaxrottest,'m');

%test plot time of maximum signal, approximately 
%n=0:1:2;
%plot(n*pi*anglerot/2/Delta/sin(anglerot/2),zeros(1,length(n)),'s')
%coseno = cos(t*2*Delta/anglerot*sin(anglerot/2))/2 + 0.5;
%plot(t,coseno,'r')
end

hold off
%title('pop in 0, red is rabi, blue is rotary, line is num, points is mm')
%title('Red: Numerics, Green: Analy'); %, Green: Without mod (original form)')
title('Red numerics, blue new formula')

% %mm
% plot(0:1/2:(length(Zmmavg)-1)/2,real(Zmmavg),'sr')
% plot(0:1/2:(length(Zmmrotavg)-1)/2,real(Zmmrotavg),'sb')
% 
% %mm no noise
% plot(0:1/2:(length(Zmmnn)-1)/2,real(Zmmnn),'dk')
% plot(0:1/2:(length(Zmmrotnn)-1)/2,real(Zmmrotnn),'dk')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating FT

timearray = [0:1/2/N:nbcycle]*2*anglerot/Rabi;
NFFT = 2^nextpow2(length(timearray));
Fs = 1/abs(timearray(1)-timearray(2));                    % Sampling frequency
fftx = Fs/2*linspace(0,1,NFFT/2+1);

ffty = fft((1-WRotnn)/2-mean((1-WRotnn)/2),NFFT)/length(timearray);
ffty = ffty(1:NFFT/2+1);

fftysignal = fft(signal-mean(signal),NFFT)/length(signal);
fftysignal = fftysignal(1:NFFT/2+1);

fftyP = fft(roundrot-mean(roundrot),NFFT)/length(roundrot);
fftyP = fftyP(1:NFFT/2+1);

fftyD = fft(diffe - mean(diffe),NFFT)/length(diffe);
fftyD = fftyD(1:NFFT/2+1);

fftyDA = fft((diffe-addi)-mean(diffe - addi),NFFT)/length(diffe);
fftyDA = fftyDA(1:NFFT/2+1);

figure(345)
plot(fftx/1e6-Rabi/2/pi/1e6,abs(fftyDA), 'g')


%%% PLOT FFT
figure(1002);
hold on
mais = max(abs(ffty));
menos = min(abs(ffty));
vetor = (abs(ffty) - menos)/(mais - menos);
mais2 = max(abs(fftysignal));
menos2 = min(abs(fftysignal));
vetor2 = (abs(fftysignal) -menos2)/(mais2-menos2);

mais3 = max(abs(fftyP));
menos3 = min(abs(fftyP));
vetorP = (abs(fftyP) - menos3)/(mais3 - menos3);

mais4 = max(abs(fftyD));
menos4 = min(abs(fftyD));
vetorD = (abs(fftyD) - menos4)/(mais4 - menos4);

soma = roundrot + addi - 0.5;
fftyS = fft(soma - mean(soma),NFFT)/length(soma);
fftyS = fftyS(1:NFFT/2+1);
mais5 = max(abs(fftyS));
menos5 = min(abs(fftyS));
vetorS = (abs(fftyS) - menos5)/(mais5 - menos5);

%plot(fftx/1e6-Rabi/2/pi/1e6,vetor,'r')
%plot(fftx/1e6-Rabi/2/pi/1e6,vetor,'r',fftx/1e6-Rabi/2/pi/1e6,vetorP,'b',fftx/1e6-Rabi/2/pi/1e6,vetorD,'g'); %normalized
plot(fftx/1e6-Rabi/2/pi/1e6,vetor,'r',fftx/1e6-Rabi/2/pi/1e6,vetorS,'b'); %normalized
%plot(fftx/1e6-Rabi/2/pi/1e6,abs(ffty),'r',fftx/1e6-Rabi/2/pi/1e6,abs(fftysignal),'g');
%red is numeric fft,green is formula "signal" FFT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% integration function
function dy = Propagate_phase0(t,y)

global Rabi g1 g2 Delta DabEq aux nb_aux noise;

Eps =noise(aux,nb_aux);

% Optical Bloch Equations

dy = zeros(3,1);    % a column vector

%y(1) u
%y(2) v
%y(3) w

dy(1) =  +y(2)*(Delta+Eps) -y(1)*g2;
dy(2) =  -y(1)*(Delta+Eps) -y(3)*Rabi - y(2)*g2;
dy(3) = -(y(3) - DabEq)*g1 + Rabi*y(2);

function dy = Propagate_phase180(t,y)

global Rabi g1 g2 Delta DabEq aux nb_aux noise;

Eps =noise(aux,nb_aux);

% Optical Bloch Equations

dy = zeros(3,1);    % a column vector

%y(1) u
%y(2) v
%y(3) w

dy(1) =  +y(2)*(Delta+Eps) -y(1)*g2;
dy(2) =  -y(1)*(Delta+Eps) +y(3)*Rabi - y(2)*g2;
dy(3) = -(y(3) - DabEq)*g1 - Rabi*y(2);

% integration function
function dy = Propagate_phase0nn(t,y)

global Rabi g1 g2 Delta DabEq aux nb_aux noise;

% Optical Bloch Equations

dy = zeros(3,1);    % a column vector

%y(1) u
%y(2) v
%y(3) w

dy(1) =  +y(2)*(Delta) -y(1)*g2;
dy(2) =  -y(1)*(Delta) -y(3)*Rabi - y(2)*g2;
dy(3) = -(y(3) - DabEq)*g1 + Rabi*y(2);

function dy = Propagate_phase180nn(t,y)

global Rabi g1 g2 Delta DabEq aux nb_aux noise;

% Optical Bloch Equations

dy = zeros(3,1);    % a column vector

%y(1) u
%y(2) v
%y(3) w

dy(1) =  +y(2)*(Delta) -y(1)*g2;
dy(2) =  -y(1)*(Delta) +y(3)*Rabi - y(2)*g2;
dy(3) = -(y(3) - DabEq)*g1 - Rabi*y(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

