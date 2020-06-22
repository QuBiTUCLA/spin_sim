close all

fighandle=openfig('Rabinoisetofit.fig');
ax=findall(fighandle,'Type','line');
x=get(ax,'Xdata');
y=get(ax,'YData');

%some params Rotarynoisetofit3
%no T1
g2 = 1/10e-6;
Rabi =  2*pi*10e6;
Delta = 3/100*Rabi;
stdev = (1/100*Rabi)^2;
anglerot = 1*pi;
timeunit =anglerot/Rabi;
tauC = 20*timeunit;

factor = 4*Rabi^2/stdev^2/tauC;

expfitnn = 0.5*(1 + exp(-x{2}*timeunit/(1/g2)));
% in the absence of noise, rotary decays as exp(-t/T_2)

%with noise
%expfitnoise =0.5 + 0.5*exp(-x{2}*timeunit/(1/g2)).*(exp(-x{2}*timeunit/(250*tauC)));
expfitnoise = 0.5*(1 + exp(-x{2}*timeunit/(factor)));
%can be that factor that multiplies tauC is Delta/sqrt(stdev)
%not really

figure(11);
plot(x{1},y{1},':b',x{2},y{2},'b',x{2},expfitnn,'r',x{2},expfitnoise,'g');