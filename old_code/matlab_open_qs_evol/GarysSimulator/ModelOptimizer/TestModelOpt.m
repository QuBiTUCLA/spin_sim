%Test script of model opt
clear all
ct=1;
for i=1:3
settings{ct} = 'params.timeStep = 1e-8;';ct=ct+1;
settings{ct} = 'params.timeStep = 5e-9;';ct=ct+1;
settings{ct} = 'params.timeStep = 3e-9;';ct=ct+1;
settings{ct} = 'params.timeStep = 1e-9;';ct=ct+1;
end

for i=1:length(settings)
    settings{i}
    [results{i},parameters{i}] = LauncherModelOpt2(settings{i});
end
% save('NoiseFieldTest1','results','parameters');

%% plot
% resultIdx = 1;
% map=mapping{resultIdx};
% mapCost=mappingCost{resultIdx};
% bestPts=bestPoints{resultIdx};
% bestCts=bestCosts{resultIdx};
% pms=params{resultIdx};
% 
% %Scatter plot of all points (3D map)
% figure(2)
% [map, mapIdx] = unique(map.','rows');
% map = map.';
% mapCost = mapCost(mapIdx);
% 
% showVal = log(mapCost); %Cost normalization to see more clearly 
% showVal = showVal - min(showVal) + 1;
% 
% scatter3(pms.NVC*map(1,:),pms.NVC*map(2,:),...
%     pms.NVC*map(3,:),showVal,showVal,'filled');
% text(pms.NVC*bestPts(1,1),pms.NVC*bestPts(2,1),...
%     pms.NVC*bestPts(3,1),'\leftarrow\fontsize{16}{\color{red}Best}');
% xlabel('SzIx (Hz)');
% ylabel('SzIy (Hz)');
% zlabel('SzIz (Hz)');
% title('NV-Carbon 1 dipole tensor search map');
