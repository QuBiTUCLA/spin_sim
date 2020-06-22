function sinexp_fit(file)

% I assume your data variables are x and y, change accordingly
% 
% %fit sin not decaying
%  myfun = @(p,x) p(1)*sin(2*pi*p(2)*x + p(3)) + p(4)
%  pinit = [0.5, 4.8e6, -pi/2, 0.7];
%  LB = [0, 1e6, -pi, 0];
%  UB = [1, 7e6, pi, 1];
%  [pbest,delta_p]=easyfit(x, y, pinit, myfun, LB, UB);

% 
% openfig(file,'new');
% h = gcf;
% line = findall(h, 'Type', 'Line'); 
% x = get(line(1), 'xdata');
% y = get(line(1), 'ydata');
% close(gcf);
% % 
% %close all;
% figure(1);
% plot(x, y);



load(file);

yref =Scan.ExperimentData{1}{1};
yreal = Scan.ExperimentData{1}{2};
%y = yreal./yref;

y = yreal;
% if want to use moving average
%lagpts = 200;
%ymvavg = tsmovavg(yref, 's', lagpts);
%y(lagpts:1:end) = yreal(lagpts:1:end)./ymvavg(lagpts:1:end);
% if not just comment


x = Scan.vary_begin:((Scan.vary_end -Scan.vary_begin)/(Scan.vary_points-1)):Scan.vary_end;




%close all;
figure(1000);
plot(x, y);

% 
% %fit decaying sin
% myfun = @(p, x) p(1) * sin(2*pi*p(2)*x + p(3)).*exp(-x/p(4)) + p(5)
% 
% % initial values 
% pinit = [0.5, 2e6, 0, 1e-6, 0.5];
% 
% % bounds for fitting parameters 
% LB = [0, 1e6, 0, 1e-9, 0];
% UB = [1, 3e6, 2*pi, 5e-3, 1];
% 
% [pbest,delta_p]=easyfit(x, y, pinit, myfun, LB, UB);


%fit decaying sin
myfun = @(p, x) p(1) * sin(2*pi*p(2)*x + p(3)) + p(4)

% initial values 
pinit = [4, 10e6, 0,mean(y)];

% bounds for fitting parameters 
LB = [1, 0.1e6, 0,8];
UB = [10, 20e6,2*pi,16];

[pbest,delta_p]=easyfit(x, y, pinit, myfun, LB, UB);

%close all;
figure(1001);
plot(x, yreal);
figure(1002);
plot(x, yreal./yref);

