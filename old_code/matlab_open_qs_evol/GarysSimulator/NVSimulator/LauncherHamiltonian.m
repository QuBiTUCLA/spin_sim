%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Launcher: Hamiltonian creation %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Gary Wolfowicz, March 2011 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LauncherHamiltonian(NVSpinNumber,numberOfCarbon,carbonDistRadius,B0)
%Create full Hamiltonian for NV
%NVSpinNumber = 1 or 3/4 for 1 collapse to 1/2 spin
%numberOfCarbon = number of carbon in the bath
%carbonDistRadius = carbons are chosen randomly within [(1) (2)]
%B0 = static field strength in Tesla
%Solution is saved in 'HnatData' file

%Checking
carbonDistRadius = max(1.615, carbonDistRadius); %1.615 is the smallest radius given by Gali
if(carbonDistRadius(1) > carbonDistRadius(2))
    cdr = carbonDistRadius(2);
    carbonDistRadius(2) = carbonDistRadius(1);
    carbonDistRadius(1) = cdr;
end

%% Constants
%B0 given in Tesla
mu0 = 4*pi*1e-7; %H.m-1
hbar = 1.0545715964207855e-34; %rad 
gammaN14 = 3.0766e6; %N14 gyromagnetic value in Hz/T (known value)
gammaE = -2.8024953e10; %Hz/T
gammaCarbon = 10.705e6; %carbon13 known value, Hz/T
ZFS_NV = 2.87e9; %Hz

%% Create Hamiltonian (IS IN UNIT OF FREQUENCY HZ)
Hnat = Hamiltonian();

%Create applied static field
BStatic = Field([0 0 B0]); % field in Tesla
BStatic.Name = 'BStatic';
Hnat.addField(BStatic);

%Load hyperfine (dipole) tensor data
load('NVCarbonHyperfine'); %Creates a 'carbonData' variable containing the Hyperfine
                           %values from ab initio calculation

%Create NV and numberOfCarbon carbon 13
SpinCreation;

%Create all interactions
InteractionCreation;

%Use the defined interactions to create a full, exact Hamiltonian
Hnat.createFullHamiltonian();

%Create Rotating frame
Hnat.createRotatingFrame();

%Save Hamiltonian
HnatData = Hnat;
save('HnatData','HnatData','carbonListPosition','NVCarbonTensor');
end