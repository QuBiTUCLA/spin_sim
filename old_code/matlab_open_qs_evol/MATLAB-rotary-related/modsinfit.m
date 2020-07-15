%function modsinfit(file)

%file = '1DExp-seq-Rabi-vary-length_rabi_pulse-2011-08-28-095700';
file = 'Rabi2011-09-24-072744';

load(file);

rabi = Scan.ExperimentData{1}{1};
lowref = Scan.ExperimentData{1}{2};
oneref = Scan.ExperimentData{1}{3};
zeroref = Scan.ExperimentData{1}{4};
%rabi = Scan.ExperimentData{1}{5}; %%%%%%% 5

x = Scan.vary_begin:((Scan.vary_end -Scan.vary_begin)/(Scan.vary_points-1)):Scan.vary_end;
 
rabinorm = (rabi - oneref)./(zeroref - oneref); %rabi normalized


y =rabinorm;
%yreal = Scan.ExperimentData{1}{2};
%y = yreal./yref;

%if want to use moving average
%lagpts = 200;
%ymvavg = tsmovavg(yref, 's', lagpts);
%y(lagpts:1:end) = yreal(lagpts:1:end)./ymvavg(lagpts:1:end);
%%if not just comment

%close all;
%plot(x, y);


%fit modulated
myfun = @(p, x) p(1)*sin(2*pi*p(2)*x + p(3)) .* sin(2*pi*p(4)*x + p(5)) + p(6);

% initial values 
pinit = [0.6, 10e6, 0, 1.5e5, pi, mean(y)];

% bounds for fitting parameters 
LB = [0.5, 9e6,0 , 0.8e5 , 0,0];
UB = [1, 11e6,2*pi, 2e5, 2*pi,1];

[pbest,delta_p]=easyfit(x, y, pinit, myfun, LB, UB);

%%close all;
figure(1001);
subplot(3,1,1);
plot(x, y);
subplot(3,1,2);
plot(x, y, 'ro', x, myfun(pbest, x), 'k-');
title('fit')
subplot(3,1,3);
plot(x,y-myfun(pbest, x));
title('residuals')
text(0,-3,[func2str(myfun) sprintf('\n %d, %d, %d, %d, %d, %d', pbest(1), pbest(2),pbest(3),pbest(4),pbest(5),pbest(6))])