

Rabi = 2*pi*20e6;
sigma = 20/100*Rabi;
tauc = 200e-9;
anglerot = pi;

cond1 = 0;
cond2 = 0;

isfast = 0;
isslow = 0;

if tauc*sigma < anglerot/2
   cond1 = 1;
end

if tauc > anglerot/2/Rabi
   cond2 = 1;
end

if 1/tauc > sigma^2/Rabi
   isfast = 1; 
end

if 1/tauc < sigma^2/Rabi
   isslow = 1; 
end

candoOUZnoiseslow = 0;
 if cond1 + cond2 + isslow == 3
    candoOUZnoiseslow = 1 
    cond1
    cond2
    isslow
 end
candoOUZnoiseslow
cond1
    cond2
    isslow
    tauc*sigma/(anglerot/2)

klklk
 
 
um = cond1
dois = cond2

fast = isfast
slow = isslow

joint3 = isslow + cond1 + cond2 

sigma*tauc

%cond3 = 0;
% if sigma*tauc < 1 %kubo number
%     cond3 = 1;
% end
%tres = cond3