clear all;

open('FactorXRatio.fig');

h = gcf;

line = findall(h, 'Type', 'Line');

x = get(line(16), 'xdata');

onepi = get(line(20), 'ydata');
fivepi = get(line(18), 'ydata');
ninepi = get(line(16), 'ydata');


close all;

figure(42);


% threepi = get(line(19), 'ydata');
% 

% 
% sevenpi = get(line(17), 'ydata');
% 

%elevenpi = get(line(7), 'ydata');

x = x(30:1:end);
onepi = onepi(30:1:end);
fivepi = fivepi(30:1:end);
ninepi = ninepi(30:1:end);
%plot(x,onepi,'b',x,fivepi,'r')

save('bluuup.mat','x','fivepi');


% this is the fitting function
%proposed

myfun = @(p, x) p(1)*pi*x + p(2) + pi*p(3)*x.^(p(4)).*cos(pi*p(5)*x.^(p(6))).^2;  
% initial values
pinit = [0.000495, 0.454, 42e-6, 1.59, 1.22, 0.5];
LB =    [0.000490, 0.450, 30e-6,    1,   -3, 0.45];
UB =    [0.000500, 0.460, 60e-6,    2,    3, 0.55];


[pbest,delta_p]=easyfit(x,fivepi, pinit, myfun, LB, UB)
