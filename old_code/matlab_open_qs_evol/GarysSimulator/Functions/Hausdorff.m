function out = Hausdorff(A,B,dt,order)
%HAUSDORFF series expansion according to 2nd order grape paper
%the exp part is not included

out = B;

nestCom = B;
for ctOrd = 1:order
    nestCom = A*nestCom-nestCom*A;
    out = out + (1i*dt).^(ctOrd)./prod(1:ctOrd+1)*nestCom;
end

out = -1i.*dt.*out;

end

