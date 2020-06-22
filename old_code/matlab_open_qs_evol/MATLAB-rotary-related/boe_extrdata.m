clear all;

open('FactorXRatio.fig');

h = gcf;

line = findall(h, 'Type', 'Line');

x = get(line(16), 'xdata');

realonepi = get(line(21), 'ydata'); %Real 1pi % fit does not work
onepi = get(line(20), 'ydata'); %actually 3pi
fivepi = get(line(18), 'ydata'); %actually 7pi
ninepi = get(line(16), 'ydata'); %actually 11pi
sevenpi = get(line(17), 'ydata'); %actually 9pi
threepi = get(line(19), 'ydata'); %actually 5pi
elevenpi = get(line(15), 'ydata'); %actually 13pi

close all;

figure(42);



x = x(30:1:end);
onepi = onepi(30:1:end);
fivepi = fivepi(30:1:end);
ninepi = ninepi(30:1:end);
sevenpi = sevenpi(30:1:end);
realonepi = realonepi(30:1:end);
threepi = threepi(30:1:end);
elevenpi = elevenpi(30:1:end);
%plot(x,onepi,'b',x,fivepi,'r')

%save('bluuup.mat','x','fivepi');


% this is the fitting function
%proposed


x = x - min(x);

myfun = @(p, x) p(1)*x + p(2);
% initial values
pinit = [1, 0];

fivepi = onepi;

[pbest,delta_p]=easyfit(x,fivepi, pinit, myfun, [], [])

figure(2);
plot(x, fivepi - myfun(pbest, x), 'o');

y = fivepi - myfun(pbest, x);
y = y;


% myfun2 = @(p, x) -p(1)*x.^(p(2)).*cos(p(3)*x.^(p(4)));
% 
% pinit = [0.004, 0.05, -6.08, 0.3];
% LB =    [0.003, 0.05, -10.0, 0];
% UB =    [0.005, 0.05,  10.0, 1];

myfun2 = @(p, x) -0.00253*x.^(0.5).*cos(p(1)*x.^(p(2)));

pinit = [ -5.92, 0.38];
LB =    [ -15.0,  0];
UB =    [ 15.0,   1];

%below, BOes's
%myfun2 = @(p, x) p(1) * x.^(0.5) .* (       -cos( (p(3) - p(2) * x) .* x )       );

% initial values boe's
% pinit = [0.005, 0.01, 0.8];
% LB =    [0.001, 0.001, -pi];
% UB =    [0.010, 0.2, pi];

[pbest2,delta_p]=easyfit(x,y, pinit, myfun2, LB, UB)

figure(100);
plot(x, myfun2(pbest2, x) + myfun(pbest, x), x, fivepi, 'o');

