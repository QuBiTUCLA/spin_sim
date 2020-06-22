%function modsinfit(file)

file = 'long.fig';

openfig(file,'new');
h = gcf;
line = findall(h, 'Type', 'Line'); 
x = get(line(1), 'xdata');
rabi = get(line(5), 'ydata');
zeroref = get(line(2),'ydata');
oneref = get(line(3),'ydata');
close(gcf);

%close all;
%figure(1000);
%plot(x, rabi,'k',x,zeroref,'b',x,oneref,'r');

y = (rabi - oneref)./(zeroref - oneref); %y is rabinorm

figure(2000);
plot(x,y);

figure(2003);
%fit decaying sin
myfun = @(p, x) p(1)*sin(2*pi*p(2)*x + p(3)) .* sin(2*pi*p(4)*x + p(5)) + p(6)

% initial values 
pinit = [1, 7e6, pi, 0.4e6, pi, mean(y)];

% bounds for fitting parameters 
LB = [0.3, 5e6,0 , 9e4 , 0,0];
UB = [1.5, 10e6,2*pi, 5e5, 2*pi,1.2];

[pbest,delta_p]=easyfit(x, y, pinit, myfun, LB, UB);

%%close all;
figure(1001);
subplot(2,1,1);
plot(x, y, 'ro', x, myfun(pbest, x), 'k-');
title('fit')
subplot(2,1,2);
plot(x,y-myfun(pbest, x));
title('residuals')
text(0,-1.5,[func2str(myfun) sprintf('\n %d, %d, %d, %d, %d, %d', pbest(1), pbest(2),pbest(3),pbest(4),pbest(5),pbest(6))])