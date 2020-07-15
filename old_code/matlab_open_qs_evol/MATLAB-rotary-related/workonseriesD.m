clear all
load('signalD.mat') % time and signal vectors

data = signal;
datacomplete = signal;

myfun0 = @(p, t) p(1) * cos(2*pi*p(2)*t + p(3)*pi) + p(4);

%remove peak 2Rabi
pinit = [0.5, 35.7e6, 0, 0.5];
% bounds for fitting parameters 
LB = [  0.05,   35.6e6,      0,       0.1];
UB = [  1,   35.8e6,      2,       1];
figure(999)
[pbest1,delta_p]=easyfit(time, data, pinit, myfun0, LB, UB);

best1 = myfun0(pbest1,time);
data = data - best1;

%remove peak 4Rabi
myfun = @(p, t) p(1) * cos(2*pi*p(2)*t + p(3)*pi);
pinit = [0.5, 71.43e6, pbest1(3)];
% bounds for fitting parameters 
LB = [  0.01,   71.4e6,      0];
UB = [  1,   71.45e6,     2];
figure(999)
[pbest2,delta_p]=easyfit(time, data, pinit, myfun, LB, UB);

best2 = myfun(pbest2,time);
data = data - best2;

periodorecipe(data,time)