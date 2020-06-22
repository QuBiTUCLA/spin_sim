close all

fighandle=openfig('Rotarynoisetofit3.fig');
ax=findall(fighandle,'Type','line');
x=get(ax,'Xdata');
y=get(ax,'YData');

% %some params Rotarynoisetofit
% %no T1, no T2
% Rabi =  2*pi*10e6;
% Delta = 5/100*Rabi;
% stdev = (1/100*Rabi)^2;
% anglerot = 1*pi;
% timeunit =anglerot/Rabi;
% tauC = 20*timeunit;

% %some params Rotarynoisetofit2
% %no T1
% g2 = 1/10e-6;
% Rabi =  2*pi*10e6;
% Delta = 5/100*Rabi;
% stdev = (1/100*Rabi)^2;
% anglerot = 1*pi;
% timeunit =anglerot/Rabi;
% tauC = 20*timeunit;

%some params Rotarynoisetofit3
%no T1
g2 = 1/10e-6;
Rabi =  2*pi*10e6;
Delta = 3/100*Rabi;
stdev = (1/100*Rabi)^2;
anglerot = 1*pi;
timeunit =anglerot/Rabi;
tauC = 20*timeunit;

% %some params Rotarynoisetofit4
% %no T1
% g2 = 0; %no T2
% Rabi =  2*pi*10e6;
% Delta = 3/100*Rabi;
% stdev = (1/100*Rabi)^2;
% anglerot = 1*pi;
% timeunit =anglerot/Rabi;
% tauC = 20*timeunit;

expfitnn = 0.5*(1 + exp(-x{2}*timeunit/(1/g2/1.5)));
% in the absence of noise, rotary decays as exp(-3t/(2T_2))

%with noise
expfitnoise =0.5 + 0.5*exp(-x{2}*timeunit*3/(2/g2)).*(exp(-x{2}*timeunit/(3*tauC)));
%can be that factor that multiplies tauC is Delta/sqrt(stdev)
%not really

figure(11);
plot(x{1},y{1},':b',x{2},y{2},'b',x{2},expfitnn,'r',x{2},expfitnoise,'g');