%close all;
%clear all;
% 
% fit_par_arr = [-0.1:0.05:0.5]; %[-0.1:0.005:0.5];
% rabi_arr = 2*pi*1e6*[10:20:30];        %[5:1:40];
% angle_par = pi*[1/2 5];             %[3 5 7 13];
% 
% color = {'--b*', '--r+','--ko', '--g.', '--ys','--cx','--md','--b'};
% 
% figure(40)
% hold on
% 
% ARRD = {};
% 
% for q = 1:1:length(angle_par)
% 
% ARRD{q} = zeros(1,length(rabi_arr));
% 
% end
% 
% nn = 1;
% for n = angle_par
% 
%         %test
%         tryarr = zeros(1,length(rabi_arr));
%    
%         qq = 1;
% for m = rabi_arr
%    
% 
%     
%    m/(2*pi*1e6)
%    arr = [];
%    for k = fit_par_arr
%       
%       [a, b] = numericalexactBOE(m, k,n);
%       arr(end+1) = sum((a-b).^2);
%    end
%    
%    [c,d] = min(arr);
%    
%    fit_par_arr(d)
%    ARRD{nn}(end) = fit_par_arr(d);
%    %test
%    tryarr(qq) = fit_par_arr(d);
%    
%    qq = qq + 1;
% end
% 
% hold on
% %plot(rabi_arr/(2*pi*1e6),ARRD{nn},color{nn})
% plot(rabi_arr/(2*pi*1e6),tryarr,color{nn})
% hold on
% 
% nn = nn + 1;
% 
% end
% 
% %hold off
% ylabel('factor that needs to be added')
% xlabel('Ratio Rabi/Delta')
% %legend('pi/2','pi','3pi', '5pi','7pi','9pi','11pi','13pi')
% legend('3pi', '5pi','7pi','13pi')

%fit_par_arr = 16+[-0.05:0.005:0.05];
%fit_par_arr = 32+[-0.05:0.005:0.05];

% res = [];
% 
% for m = rabi_arr
% 	hlp = [];
% 	tic
% for k = fit_par_arr
% 
% 	k
% 	[a, b] = numericalexactBOE(m, k);
% 	hlp(end+1) = sum((a-b).^2);
% end;
% 	res = [res.' hlp.'].';
% 	save('res');
% 	toc
% end;
% 
% imagesc(fit_par_arr, rabi_arr/2/pi/1e6, res);
% xlabel('Fit Parameter');
% ylabel('Rabi Frequency');
% 
fit_par_arr = [-1:0.05:30];
arr = [];

for xx = fit_par_arr
xx
[a, b] = numericalexactBOE(xx);

arr(end+1) = sum((a-b).^2);

end


figure(12)
plot(fit_par_arr,arr)
% 
% [c,d] = min(arr);
%     
%   fit_par_arr(d)
