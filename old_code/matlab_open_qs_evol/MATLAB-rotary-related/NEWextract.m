clear all;

%OLD FIG

open('FactorXRatio.fig');

h = gcf;

line = findall(h, 'Type', 'Line');

x1 = get(line(16), 'xdata');
onepi = get(line(21), 'ydata'); %Real 1pi % fit does not work
threepi = get(line(20), 'ydata'); %actually 3pi
fivepi = get(line(19), 'ydata'); %actually 5pi

x1 = x1(1:1:20);
onepi = onepi(1:1:20);
threepi = threepi(1:1:20);
fivepi = fivepi(1:1:20);

close all;

open('NEWRatio.fig');

h = gcf;

line = findall(h, 'Type', 'Line');

x = get(line(8), 'xdata');

cincopi = get(line(7), 'ydata'); %Real 1pi
trespi = get(line(8), 'ydata');
umpi = get(line(9), 'ydata');

close all;

open('Ratio-end.fig');

h = gcf;

line = findall(h, 'Type', 'Line');

x2 = get(line(8), 'xdata');

cincopi2 = get(line(7), 'ydata'); %Real 1pi
trespi2 = get(line(8), 'ydata');
umpi2 = get(line(9), 'ydata');

close all;

figure(42);

ONE = [onepi umpi umpi2];
THREE = [threepi trespi trespi2];
FIVE = [fivepi cincopi cincopi2];
XIS = [x1 x x2];

plot(XIS,ONE,'b',XIS,THREE,'r',XIS,FIVE,'k')
axis([XIS(1) XIS(end) min(FIVE) max(FIVE)])

save('bluuup.mat','XIS','ONE','THREE','FIVE');

figure(43);


%DATA TO BE FITTED
fitdata = THREE;
x = XIS;

%x = x - min(x);

%try expon
%myfun = @(p, x) p(1) + p(2) * x;
% initial values
%pinit = [0.45, 0.0017];
%LB = [0.4, 0. ];
%UB = [0.5, 0.002 ];

% myfun = @(p, x) p(1)*x + p(2);
% % initial values
% pinit = [1, 0];
% 

%[pbest,delta_p]=easyfit(x,fitdata, pinit, myfun, LB, UB)

%figure(2);
%plot(XIS,fitdata - myfun(pbest, x), 'o');

%y =fitdata - myfun(pbest, x);
%y = y;

%%%

new = fitdata -( -exp(-0.281*x) + 0.00169*x + 0.461 );

figure(2);
plot(XIS, new,'o')

a



%%%

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
plot(x, myfun2(pbest2, x) + myfun(pbest, x), x, fitdata, 'o');

