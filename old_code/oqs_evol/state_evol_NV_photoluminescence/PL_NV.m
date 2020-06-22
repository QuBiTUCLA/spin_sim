close all;

% calculate stuff using Liou_NV

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% States sorted by increasing energy
N = 7;

bs = cell(N, 1);
for k = 1:N
	bs{k} = basis(N, k);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time vector
tinit = 0;
tstep = 1e-9;
tfinal = 2500e-9;
tlist = tinit:tstep:tfinal;

%%% init zero_g
rho0 = bs{1}*bs{1}';
rho = ode2es(Liou_NV, rho0);

tpops = cell(N);

for k = 1:N
		tpops{k} = real(esval(expect(bs{k}*bs{k}', rho), tlist));
end

figure(1);
plot(tlist,tpops{1},tlist,tpops{2},tlist,tpops{3},tlist,tpops{4},tlist,tpops{5},tlist,tpops{6},tlist,tpops{7})
xlabel('Time (s)');
ylabel('Population, starting at zero_g')
legend('zero_g', 'm_one_g', 'one_g', 'meta', 'zero_e', 'm_one_e', 'one_e');

%%%% init m_one_g
rho0 = bs{2}*bs{2}';
rho = ode2es(Liou_NV, rho0);

tpops2 = cell(N);

for k = 1:N
		tpops2{k} = real(esval(expect(bs{k}*bs{k}', rho), tlist));
end


figure(11);
plot(tlist,tpops2{1},tlist,tpops2{2},tlist,tpops2{3},tlist,tpops2{4},tlist,tpops2{5},tlist,tpops2{6},tlist,tpops2{7})
xlabel('Time (s)');
ylabel('Population, starting at m_one_g')
legend('zero_g', 'm_one_g', 'one_g', 'meta', 'zero_e', 'm_one_e', 'one_e');

%%%% init superposition zero_g and m_one_g
rho0 = 0.5*bs{1}*bs{1}' + 0.5*bs{2}*bs{2}';
rho = ode2es(Liou_NV, rho0);

tpops3 = cell(N);

for k = 1:N
		tpops3{k} = real(esval(expect(bs{k}*bs{k}', rho), tlist));
end

figure(111);
plot(tlist,tpops3{1},tlist,tpops3{2},tlist,tpops3{3},tlist,tpops3{4},tlist,tpops3{5},tlist,tpops3{6},tlist,tpops3{7})
xlabel('Time (s)');
ylabel('Population, starting at sup zero_g/m_one_g')
legend('zero_g', 'm_one_g', 'one_g', 'meta', 'zero_e', 'm_one_e', 'one_e');

%%% plot
figure(2);
plot(tlist, tpops{5} + tpops{6} + tpops{7}, tlist, tpops2{5} + tpops2{6} + tpops2{7},tlist, tpops3{5} + tpops3{6} + tpops3{7})
xlabel('Time (s)');
ylabel('Instantaneous PL collected')
legend('init zero_g', 'init m_one_g', 'init superp zero_g/one_g');

% Now I want the integrated PL collected
int_PL_zero_g = zeros(1,length(tlist));
int_PL_m_one_g = zeros(1,length(tlist));
int_PL_sup = zeros(1,length(tlist));
for q = 1:1:length(tlist)
   int_PL_zero_g(q) = sum(tpops{5}(1:1:q)) + sum(tpops{6}(1:1:q)) + sum(tpops{7}(1:1:q));
   int_PL_m_one_g(q) = sum(tpops2{5}(1:1:q)) + sum(tpops2{6}(1:1:q)) + sum(tpops2{7}(1:1:q));
   int_PL_sup(q) = sum(tpops3{5}(1:1:q)) + sum(tpops3{6}(1:1:q)) + sum(tpops3{7}(1:1:q));
end

figure(3);
plot(tlist,int_PL_zero_g,tlist,int_PL_m_one_g,tlist,int_PL_sup)
xlabel('Time (s)');
ylabel('Integrated PL collected up to this time')
legend('init zero_g', 'init m_one_g', 'init superp zero_g/one_g');

figure(4);
plot(tlist,int_PL_zero_g-int_PL_m_one_g)
xlabel('Time (s)');
ylabel('Integrated PL collected up to this time')
legend('Difference init zero_g - init m_one_g');

%reproduce measured photoluminescence curve
acq_time = 100e-9;
delta = round(acq_time/tstep);
PL_zero_g = zeros(1,length(tlist)-delta);
PL_m_one_g = zeros(1,length(tlist)-delta);
PL_sup = zeros(1,length(tlist)-delta);
for p = 1:1:(length(tlist) - delta)
   PL_zero_g(p) = sum(tpops{5}(p:1:p+delta)) + sum(tpops{6}(p:1:p+delta)) + sum(tpops{7}(p:1:p+delta));
   PL_m_one_g(p) = sum(tpops2{5}(p:1:p+delta)) + sum(tpops2{6}(p:1:p+delta)) + sum(tpops2{7}(p:1:p+delta));
   PL_sup(p) = sum(tpops3{5}(p:1:p+delta)) + sum(tpops3{6}(p:1:p+delta)) + sum(tpops3{7}(p:1:p+delta));
end

figure(5);
plot(tlist(1:1:end-delta),PL_zero_g,tlist(1:1:end-delta),PL_m_one_g,tlist(1:1:end-delta),PL_sup)
xlabel('Time (s)');
ylabel(['PL curve for acq duration ' num2str(acq_time)  ', acq starts at time'])
legend('PL init zero_g', 'PL init m_one_g', 'PL init sup zero_g/m_one_g');

figure(6);
plot(tlist(1:1:end-delta),PL_zero_g -PL_m_one_g)
xlabel('Time (s)');
ylabel(['PL curve for acq duration ' num2str(acq_time)  ', acq starts at time'])
legend('Difference PL init zero_g - PL init m_one_g');

%try to find optimum acq length
l = 1;
PL = [];

acq_time_arr = [100 200 300 500 1000 1500] * 1e-9;

for acq_time = acq_time_arr   
delta = round(acq_time/tstep);
PL(l,:) = zeros(1,length(tlist));
for p = 1:1:(length(tlist) - delta)
   PL(l,p) = sum(tpops{5}(p:1:p+delta)) + sum(tpops{6}(p:1:p+delta)) + sum(tpops{7}(p:1:p+delta)) - sum(tpops2{5}(p:1:p+delta)) - sum(tpops2{6}(p:1:p+delta)) - sum(tpops2{7}(p:1:p+delta));
end
    l = l + 1;
end

figure(7);
plot(tlist,PL')
xlabel('Time (s)');
ylabel(['Difference PL init zero_g - PL init m_one_g, acq starts at time'])
entries = cellfun(@num2str, num2cell(acq_time_arr/1e-9), 'Uniformoutput', false);

legend(entries);
title('For different datection lengths (in ns)');
