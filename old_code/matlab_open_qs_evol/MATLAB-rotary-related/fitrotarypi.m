% experiments were done using 28ns as pi time, 21ns as 3pi/4 time
% == rabi is 17.86MHz

load('3pi4sig.mat'); %vectors timered and signalred

 save('3pi4RotSignal.txt','timen','signaln','-ascii');
 fid = fopen('3pi4RotSignal.txt', 'w');
 for k = 1:length(timen)
     fprintf(fid, '%d %d\n', timen(k), signaln(k));
 end;
 fclose(fid);


figure(120)
plot(timered,signalred,'k.')
axis([min(timered) max(timered) 0 max(signalred)])

%p(1) angle in pi units
%p(2) rabi in Hz
%p(3) amp nuc 0
%p(4) amp nuc +1
%p(5) amp nuc -1
%p(6) small detuning delta
%p(7) HF

rotaryfit = @(p,t) 0.5*(1 + cos(pi*p(1)/2)^2 + sin(pi*p(1)/2)^2*cos(t*(2*pi*p(2))*pi/(p(1)*pi)).*(p(3)*cos(t*2*p(6)/(pi*p(1))*sin(pi*p(1)/2))  + p(4)*cos(t*2*(p(6)+p(7))/(pi*p(1))*sin(pi*p(1)/2)) +  p(5)*cos(t*2*(p(6)-p(7))/(pi*p(1))*sin(pi*p(1)/2))));

% initial values 
pinit = [1, 17.5e6, 0.2, 0.2, 0.2, 0.1e6, 2e6];

% bounds for fitting parameters 
LB = [0.8,   17e6,      0.1,    0.1, 0.1, 0, 2e6];
UB = [1.2,   18e6,      0.5,    0.5, 0.5, 1.6e6, 2.2e6];

[pbest,delta_p]=easyfit(timered, signalred, pinit, rotaryfit, LB, UB);

% %simplified
% rotaryfit = @(p,t) cos(t*(2*pi*17.5e6)*pi).*(p(1)*cos(t*2*p(4)/pi)  + p(2)*cos(t*2*(p(4)+p(5))/pi) +  p(3)*cos(t*2*(p(4)-p(5))/pi));
% 
% %p(1) amp nuc 0
% %p(2) amp nuc +1
% %p(3) amp nuc -1
% %p(4) small detuning delta
% %p(5) HF
% 
% 
% % initial values 
% pinit = [0.2, 0.2, 0.2, 0.1e6, 2e6];
% 
% % bounds for fitting parameters 
% LB = [ 0.1,    0.1, 0.1, 0, 1.9e6];
% UB = [ 0.5,    0.5, 0.5, 1.6e6, 2.2e6];
% 
% [pbest,delta_p]=easyfit(timered, signalred, pinit, rotaryfit, LB, UB);