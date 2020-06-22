function L = Liou_NV

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% States sorted by increasing energy
N = 7;
zero_g = basis(N,1);
m_one_g = basis(N,2);
one_g = basis(N,3);
meta = basis(N,4);
zero_e = basis(N,5);
m_one_e = basis(N,6);
one_e = basis(N,7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants

ZFS_g = 2.870e9;
ZFS_e = 1.420e9;
giro = 2.8e6; %Hz/G
sat = 1; %saturation param for laser power; varies from 0 to 1

H0 = 0 * zero_g*zero_g'; % in the interaction frame;

H1 = 0 * zero_g*zero_g'; % not going to be zero when mw on

% add both Hamiltonians	
H = H0 + (H1+H1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Decay rates

C = {};

% NV original parameters
carrier_decay = 77e6;
cross_decay = 1.5e6;
ones_to_meta_decay = 30e6;
meta_to_zero_decay = 3.3e6;

%test params to alter the dynamics
%carrier_decay = 77e6;
%cross_decay = 0;
%ones_to_meta_decay = 500e6;
%meta_to_zero_decay = (3.3/2)*1e6; %if added ir laser

% zero_e to zero_g
C{end+1} = sqrt(carrier_decay) * zero_g*zero_e';

% one_e to zero_g
C{end+1} = sqrt(cross_decay) * zero_g*one_e';

% m_one_e to zero_g
C{end+1} = sqrt(cross_decay) * zero_g*m_one_e';

% meta to zero_g
C{end+1} = sqrt(meta_to_zero_decay) * zero_g*meta';

% zero_e to one_g
C{end+1} = sqrt(cross_decay) * one_g*zero_e';

% one_e to one_g
C{end+1} = sqrt(carrier_decay) * one_g*one_e';

% zero_e to m_one_g
C{end+1} = sqrt(cross_decay) * m_one_g*zero_e';

% m_one_e to m_one_g
C{end+1} = sqrt(carrier_decay) * m_one_g*m_one_e';

% one_e to meta
C{end+1} = sqrt(ones_to_meta_decay) * meta*one_e';

% m_one_e to meta
C{end+1} = sqrt(ones_to_meta_decay) * meta*m_one_e';

% zero_g to zero_e via laser
C{end+1} = sqrt(sat*carrier_decay) * zero_e*zero_g';

% one_g to one_e via laser
C{end+1} = sqrt(sat*carrier_decay) * one_e*one_g';

% m_one_g to m_one_e via laser
C{end+1} = sqrt(sat*carrier_decay) * m_one_e*m_one_g';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the Liouvillian
L = -1i * (spre(H) - spost(H)); 

for k = 1:length(C)
	L = L + make_LiouC(C{k});
end

return

% ------------------------------------------------------------------------
function Lc = make_LiouC(C0)
% Liouvillian contribution of the collapse operator C0
C0dC0 = C0'*C0;
Lc = spre(C0)*spost(C0')-0.5*spre(C0dC0)-0.5*spost(C0dC0);
% ------------------------------------------------------------------------
