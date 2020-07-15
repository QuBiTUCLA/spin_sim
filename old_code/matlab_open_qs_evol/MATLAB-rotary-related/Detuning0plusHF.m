
anglerot = pi;
Rabi = 2*pi*18e6;
%have to be the same as defined in the 2 cleancodeWITH... below

%Plot expected signal for zero detuning PLUS HF detuning

DeltaArr = [2*pi*2e6 0];

MeanZmmrotavg = [];
MeanZmmrot = [];
MeanZmmavg = [];
MeanZmm = [];

for Delta = DeltaArr

[A,B,C,D] = cleancodeGENERALIZED(Delta);

if Delta == DeltaArr(1)
MeanZmmrotavg = real(A);
MeanZmmrot = real(B);

MeanZmmavg =  real(C);
MeanZmm = real(D);
    
else
MeanZmmrotavg = MeanZmmrotavg + real(A);
MeanZmmrot = MeanZmmrot + real(B);

MeanZmmavg = MeanZmmavg + real(C);
MeanZmm = MeanZmm + real(D);

end

end

MeanZmmrotavg = MeanZmmrotavg/length(DeltaArr);
MeanZmmrot = MeanZmmrot/length(DeltaArr);

MeanZmmavg = MeanZmmavg/length(DeltaArr);
MeanZmm = MeanZmm/length(DeltaArr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting

figure(999); 

hold on
%matrix multiplication result rotary
plot([0:1/2:(length(MeanZmmrotavg)-1)/2]*2*anglerot/Rabi,real(MeanZmmrotavg),'b')
%matrix multiplication result rotary no noise
plot([0:1/2:(length(MeanZmmrot)-1)/2]*2*anglerot/Rabi,real(MeanZmmrot),'b*')
%matrix multiplication result rabi
plot([0:1/2:(length(MeanZmmavg)-1)/2]*2*anglerot/Rabi,real(MeanZmmavg),'r')
%matrix multiplication result rabi no noise
plot([0:1/2:(length(MeanZmm)-1)/2]*2*anglerot/Rabi,real(MeanZmm),'r*')
%prb results
%plot(0:1/20:nbcycle,prb,'k*')
hold off

%title('Rotary: pop in 0, blue is rot matrix mult, black is rot PRB approx, red is rabi matrix mult')
title('blue is rot matrix mult, red is rabi matrix mult for sum of 0 + HF detunings')
xlabel(['time, for angle' sprintf('%0.3f * pi', anglerot/pi)])

