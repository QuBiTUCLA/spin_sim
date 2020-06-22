function [] = fitSwXwRabi() 

%close all
clear all
load('1DExp-seq-Rabi-SwXw-vary-mw_freq-2012-01-16-011901Uno-Rabi164cycles.mat')
%load('1DExp-seq-Rabi-SwXw-vary-mw_freq-2012-01-21-043844Rabi164cycles.mat');
%load('rabi_sensitivity_162cycleMAsashi.mat')

rabi = Scan.ExperimentData{1}{1};
lowref = Scan.ExperimentData{1}{2};
oneref = Scan.ExperimentData{1}{3};
zeroref = Scan.ExperimentData{1}{4};

w = Scan.vary_begin:((Scan.vary_end -Scan.vary_begin)/(Scan.vary_points-1)):Scan.vary_end;

%average then normalize
rabinorm = (rabi - oneref)./(zeroref - oneref); %rabi normalized

figure(13)
plot(w,rabinorm)

figure(12); 
%subplot(2,1,1)
plot(w,rabinorm)
hold on
%fit decaying sin
myfun = @(p, w) p(7) - p(1)*(p(4)^2./((w - p(6)).^2 + p(4)^2)).*sin(pi*9.184e-6*sqrt(p(4)^2 + (p(6) - w).^2)).^2 - p(2)*(p(4)^2./((w - p(6) + p(5)).^2 + p(4)^2)).*sin(pi*9.184e-6*sqrt(p(4)^2 + (p(6) - w + p(5)).^2)).^2 -p(3)*(p(4)^2./((w - p(6) - p(5)).^2 + p(4)^2)).*sin(pi*9.184e-6*sqrt(p(4)^2 + (p(6) - w - p(5)).^2)).^2;

%p(1) amp
%p(2) amp
%p(3) amp
%p(4) Rabi
%p(5) hyperfine
%p(6) reson freq
%p(7) offset
%p(8) 'frequency' pi*time = pi*9184e-6;

pinit = [0.3, 0.3, 0.3, 17.5e6, 2.1e6,3.1569e9, 1];
% bounds for fitting parameters 
LB = [0.1,    0.1,    0.1,      17e6,1.9e6,3.15e9, 0.6];
UB = [0.5,   0.5,     0.5,  19e6,2.2e6,3.16e9, 1.2];
[pbest,delta_p]=easyfit(w, rabinorm, pinit, myfun, LB, UB);
hold off


% sigsig{aux} = rabinorm;
% freqw{aux} = w;

% figure(98)
% plot(w,rabinorm,'r.')
% 
% color = {'r', 'b', 'k', 'b', 'k'};
% color2 = {'r', 'b', 'k', 'b', 'k'};
% figure(99)
% subplot(2,1,1)
% hold on
% for kk=[1 2 5]
% plot(freqw{kk},plotfc(pbest{kk},freqw{kk}),color{kk})
% hold on
% end
% for  kk=[1 2 5]
% plot(freqw{kk},sigsig{kk},'LineStyle',':','Linewidth',0.3,'Color',color2{kk}) 
% hold on
% end
% hold off
% legend('2.8\mus','4.2\mus','9.2\mus')
% axis([3.1565e9 3.1575e9 0.45 0.8])
% %subplot(2,1,2)
end

% function curve = plotfc(p,w)
%  curve = p(1) * cos(2*pi*p(2)*(p(3) - w)) + p(4) + p(5)*cos(2*pi*p(2)*(p(6) + (p(3) - w))) + p(7)*cos(2*pi*p(2)*(p(6) - (p(3) - w)));
% 
% end
