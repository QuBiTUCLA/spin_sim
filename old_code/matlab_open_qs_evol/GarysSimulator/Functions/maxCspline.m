function maxX = maxCspline(interval,Xvalues,Yvalues,Gvalues)
%MAXCSPLINE
%This function evaluate the maximizer of a cubic hermite spline interpolation
%given an interval=[max min] with values at Xvalues=[Xa Xb] equal to Yvalues=[Ya Yb]
%and their respective gradient values Gvalues=[Ga Gb]
%May 14th 2011 Gary Wolfowicz

%Check if interval is type [min max]
Dint = Xvalues(2) - Xvalues(1);
% err = 1e2*eps(max(Xvalues));
% if(Dint <= err) %use of err because division by Dint later
%     error('Interval must be [min max] with min < max')
% end

%Evaluation of coefficient of equation: derivative of interpolation = 0
B = [Yvalues(1) Dint*Gvalues(1) Yvalues(2) Dint*Gvalues(2)];
A = B*[6 3 -6 3;-6 -4 6 -2;0 1 0 0].';

%Determinant of equation
delta = A(2)^2 - 4*A(1)*A(3);

%Evaluation of interpolation at extrema
Tint = (interval -  Xvalues(1))/Dint;
Yint = (1+2*Tint).*(1-Tint).^2*Yvalues(1) + Tint.*(1-Tint).^2*Dint*Gvalues(1) + ...
        Tint.^2.*(3-2*Tint)*Yvalues(2) + Tint.^2.*(Tint-1)*Dint*Gvalues(2);

if(delta < 0)
    [~,maxIndex] = max(Yint);
    maxX = interval(maxIndex);
    return;
end

%Evaluation of extrema for delta >= 0
err = 1e2*eps(max(B));
if(abs(A(1)) > err) %prevent numerical error (round off to zero)
    tM = (-A(2) - sqrt(delta))/(2*A(1));
    tP = (-A(2) + sqrt(delta))/(2*A(1));
elseif(abs(A(2)) > err)
    tM = -A(3)/A(2);
    tP = -A(3)/A(2);
else
    tM = 0;
    tP = 0;
end

%Evaluation of interpolation at extrema (CODE TO IMPROVE)
if(tM >= Tint(1) && tM <= Tint(2))
    YM = (1+2*tM)*(1-tM)^2*Yvalues(1) + tM*(1-tM)^2*Dint*Gvalues(1) + ...
        tM^2*(3-2*tM)*Yvalues(2) + tM^2*(tM-1)*Dint*Gvalues(2);
else 
    tM = Tint(1); %This assignment will have no influence in the max below
    YM = Yint(1);
end
if(tM ~= tP)
    if(tP >= Tint(1) && tP <= Tint(2))
        YP = (1+2*tP)*(1-tP)^2*Yvalues(1) + tP*(1-tP)^2*Dint*Gvalues(1) + ...
            tP^2*(3-2*tP)*Yvalues(2) + tP^2*(tP-1)*Dint*Gvalues(2);
    else 
        tP = Tint(1);
        YP = Yint(1);
    end
else
    YP = YM;
end

maxXs = [interval(1) interval(2) Dint*tM+Xvalues(1) Dint*tP+Xvalues(1)];
[~,maxIndex] = max([Yint(1) Yint(2) YM YP]);
maxX = maxXs(maxIndex);

% %For plot if necessary
% x=linspace(interval(1),interval(2),30);
% t=(x-Xvalues(1))/Dint;
% y = (1+2.*t).*(1-t).^2.*Yvalues(1) + t.*(1-t).^2.*Dint.*Gvalues(1) + ...
%         t.^2.*(3-2.*t).*Yvalues(2) + t.^2.*(t-1).*Dint.*Gvalues(2);
% plot(x,y)

end

