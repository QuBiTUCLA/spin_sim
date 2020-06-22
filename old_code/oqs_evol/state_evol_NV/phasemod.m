%% Definitions

% hbar = 1 

% those are spin operators for spin 1
I3=speye(3);
Z3=sparse(diag([1 0 -1])); 
X3=sparse([0 1 0;1 0 1;0 1 0])/sqrt(2);
Y3=sparse([0 -1i 0;1i 0 -1i;0 1i 0])/sqrt(2);

% those are spin operators for spin 1/2
I2=speye(2);
Z2=(1/2)*sparse(diag([1 -1])); 
X2=(1/2)*sparse([0 1;1 0]);
Y2=(1/2)*sparse([0 -1i;1i 0]);

%% NV params

De=1420; % ZFS excited state, MHz
Dg=2870; % ZFS ground state, MHz
ge=17.588/2/pi; % (-1*) gyromagnetic ratio electron/(2*pi), MHz/G 

% Transition strengths
G=77; % direct decay by fluorescence, sublevel preserving, for 0 and \pm 1 sublevels, MHz 
G1m=30; % decay to 0 via singlet, for \pm 1 sublevels, MHz
Gm0=3.333; % decay to 0 via singlet, for 0 sublevel, MHz  
Gx=1.5; % direct decay to 0 by fluorescence, for \pm 1 sublevels, MHz 

% Rabi frequency as given by the excitation laser 
W=G*1.1;
Wx=Gx*W/G; 

%% Hamiltonian

Ne=3; % number of electron levels in the ground state

Ng=Ne; % total number of levels in the e in the ground state 

N=7; % 7 total number of states

% possible values for magnetic field
% We are aiming at the LACs in the excited or gd state of the electron
%B=(Dg)/ge; 
B=800;

% Zeeman and ZFS terms
Hg=kron(Dg*Z3*Z3+ge*B*Z3*pi,speye(Nn));
He=kron(De*Z3*Z3+ge*B*Z3*pi,speye(Nn));

% build block-diagonal Hamiltonian
H=blkdiag(Hg,He); 
clear Hg He

%% Master equation

% build dissipation operators

for k=1:Ne
    Ls{k}=sqrt(G)*sparse(k,k+Ne,1,2*Ne+1,2*Ne+1);
end

Ls{Ne+1}=sqrt(G1m)*sparse(2*Ne+1,Ne+1,1,2*Ne+1,2*Ne+1);
Ls{Ne+2}=sqrt(Gm0)*sparse(Ne-1,2*Ne+1,1,2*Ne+1,2*Ne+1);
if Ne>2
    Ls{Ne+3}=sqrt(G1m)*sparse(2*Ne+1,2*Ne,1,2*Ne+1,2*Ne+1);
end

% build excitation operators
Ld=Ls;

for k=1:Ne
    Ls{length(Ld)+k}=sqrt(W)*sparse(k+Ne,k,1,2*Ne+1,2*Ne+1),speye(Nn));
end

% build master equation

Glight=-1i*(kron(H,speye(N))-kron(speye(N),H))-Lind(Ls); % for Ne=3, Ls(1) until Ls(9); d(rho)/dt with laser
Gdark =-1i*(kron(H,speye(N))-kron(speye(N),H))-Lind(Ld); % for Ne=3, Ld(1) until Ld(6) (corresponding to Ls(1) to Ls(6)); d(rho)/dt w/o laser

clear Ls Ld

%% time evolution
tdark=1.5; %us % discrete evolution for rho w/o illumination
t1ight=.5; %us % discrete evolution for rho under illumination
dt=t1+td; % used?
S=expms(Glight*t1ight); % for rho under illumination
Sd=expms(Gdark*tdark); % for rho w/o illumination

tmax=100;
for t=1:tmax
    for h=1:2
        r0=X{h}*S1*S*r0; %first illuminate for dt, then evolve in dark for 0.5 mus
    end
  pop(:,t+1)=pops(Sd*r0,N);
  if mod(t*10,tmax)==0, disp([num2str(t) '%']); end % show how much of computation is already done
end
 
% Version ssroC33
% also-result-ssroC33.jpg
time=[(0:t)*dt];
figure
subplot(1,2,1) %|0>|000>
plot(time,real(pop0(1:Ng,:))')
hold on
plot(time,real(pop0(1+Ng:2*Ng,:))','--')
plot(time,real(pop0(1+2*Ng:N,:))',':')
title('|0>|000>')
axis tight 

subplot(1,2,2) % |+1>|111>
plot(time,real(pop1(1:Ng,:))')
hold on
plot(time,real(pop1(1+Ng:2*Ng,:))','--')
plot(time,real(pop1(1+2*Ng:N,:))',':')
title('|+1>|111>')
axis tight

% only building legend
for k=1:Ne
    se=['|',num2str(k-2),' '];
    for h=1:Nn
        sleg{(k-1)*Nn+h}=[se, dec2bin(h-1,Nc),'>'];
        sleg{+Ne*Nn+(k-1)*Nn+h}=sleg{(k-1)*Nn+h};
    end
end

for h=1:Nn
        sleg{2*Ne*Nn+h}=['|s ', dec2bin(h-1,Nc),'>'];
end
legend(sleg)
 
% fluorescence intensity (at LAC) with DFS encoding
% READOUT
figure %result-ssroC33.jpg / fig 5C
plot(time,real(pe0))
hold on
plot(time, real(pe1),'r')
disp(['Total differential emission: ' num2str(abs(sum(pe0-pe1)))])

% Version polarizeC33

time=[(0:t)*dt*3]; % why times 3?
figure
plot(time,real(pop(1+Nn:2*Nn,:))')
%legend(sleg)
hold on
plot(time,real(pop(1+2*Nn:3*Nn,:))','--')
plot(time,real(pop(1:Nn,:))',':')

% plot(time,real(pop(1:Ng,:))')
% hold on
% plot(time,real(pop(1+Ng:2*Ng,:))','--')
% plot(time,real(pop(1+2*Ng:N,:))',':')
title('Polarize')
axis tight 
for k=1:Ne
    se=['|',num2str(k-2),' '];
    for h=1:Nn
        sleg{(k-1)*Nn+h}=[se, dec2bin(h-1,Nc),'>'];
        sleg{+Ne*Nn+(k-1)*Nn+h}=sleg{(k-1)*Nn+h};
    end
end

for h=1:Nn
        sleg{2*Ne*Nn+h}=['|s ', dec2bin(h-1,Nc),'>'];
end
legend(sleg)
