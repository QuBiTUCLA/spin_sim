t=0:1e-9:50e-6;

%%%

Omega = 12e6;

Rabi = 5e6;
eta = sqrt(3/2)*sqrt(1-4*Rabi^2/Omega^2)
T1 = 10; %100e-3; 
T2 = 20e-6;
%%%

R = (4*Rabi^2 -(4- 2/3*eta^2*Omega^2) + 2*1i*Omega^2*(1/T1 + 1/T2) + 1/(T1*T2));
S = (4*Rabi^2 + 1/(T1*T2) + eta^2/2/T1/T2 +(eta^4*Omega/8/T1/(1/T1 + 1/T2)));
sat = 1/T1/S/T2;
entsat = 1/T1/S*eta^2*Omega/2*imag(1 - eta^2*Omega^2/2/R);
bewegung = -1/T1/S*eta^2*Omega/2*imag((1- (1-8/(eta^2))*eta^2*Omega^2/2/R)*exp(2*1i*Omega*t)) + 1/T1/S*2/T2*real(eta^2*Omega^2/2/R*exp(2*1i*Omega*t));
analytpop =[sat + entsat+ bewegung];

figure(15); plot(bewegung./(bewegung+sat+entsat));
figure(55); plot(t,analytpop);
