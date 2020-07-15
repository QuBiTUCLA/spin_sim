% V=zeros(length(U),size(U{1},1)*size(U{1},2));
% for i = length(U)
%     V(i,:) = U{i}(:);
% end
% V=sparse(V);
% tic
% csvwrite('Ulist.dat',V)
% toc
% 
% tic
% M = csvread('Ulist.dat');
% toc
% % dlmWrite?
tic
U=propagator.unitary();
toc
