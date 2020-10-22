#
# Dynamics of radical pair subjected to Geomagnetic field
#
from qutip import *
from scipy import *
from pylab import *
import math
import csv

def run():

# problem parameters:
    pi=math.pi
    muB = 5.788*10**-5  # Bohr Magneton in eV/Tesla
    g=2 # g-factor for electron
    gama=muB*g # in eV/Tesla
    gama=gama*0.1519756*10**16 # in sec-1/Tesla (Converted using natural units, by taking hbar=1)
    #B0=47*10**-6 # Geomagnetic field in Frankfurt (Tesla)
    #B0=46*10**-6 # Geomagnetic field in Frankfurt (Tesla)
    Brf=150*10**-9 # Disturbing RF field of strength 150 nanoTesla
    w = 2*pi*1.316*10**6      # Frequency of externally applied RF field
    sqrt2=math.sqrt(2)

# Initial State
    up=basis(2,0)
    down=basis(2,1)
    singlet=(tensor(up,down)-tensor(down,up))/sqrt2 # Initial State of the radical pair
    trip0=(tensor(up,down)+tensor(down,up))/sqrt2 # Triplet state with zero spin
    tripu= tensor(up,up) # Tirplet state with spin +1
    tripd=tensor(down,down) # Triplet state with spin -1
 
    # Initial State
    e_dm=singlet*singlet.dag()
    nuc_dm=0.5*qeye(2)
    rho=tensor(nuc_dm,e_dm)
    rho0=Qobj(pad(rho.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0)))

    # For Local Coherence
    elec_dm=singlet*singlet.dag()
    elec_array=elec_dm.full()
    elec_array[1][1]=0   # Making diagonal elements of density matrix to zero
    elec_array[2][2]=0
    elec_offd=Qobj(elec_array)
    nuc_dm=0.5*qeye(2)
    rho_d=tensor(nuc_dm,elec_offd)
    rho0_offd=Qobj(pad(rho_d.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0)))

# Defining Hamiltonian operators
    S1x=tensor(qeye(2),sigmax(),qeye(2))
    S2x=tensor(qeye(2),qeye(2),sigmax())
    S1x_10=Qobj(pad(S1x.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0))) # pad() function increases the size of array of S1x from  
    #8x8 to 10x10 and put 0s at the newly created places
    S2x_10=Qobj(pad(S2x.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0)))
    Sx=S1x_10+S2x_10
    S1y=tensor(qeye(2),sigmay(),qeye(2))
    S2y=tensor(qeye(2),qeye(2),sigmay())
    S1y_10=Qobj(pad(S1y.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0)))
    S2y_10=Qobj(pad(S2y.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0)))
    Sy=S1y_10+S2y_10
    S1z=tensor(qeye(2),sigmaz(),qeye(2))
    S2z=tensor(qeye(2),qeye(2),sigmaz())
    S1z_10=Qobj(pad(S1z.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0)))
    S2z_10=Qobj(pad(S2z.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0)))
    Sz=S1z_10+S2z_10
    #ISx8=tensor(sigmax(),sigmax(),qeye(2))
    ISx8=tensor(sigmax(),qeye(2),sigmax())
    #ISy8=tensor(sigmay(),sigmay(),qeye(2))
    ISy8=tensor(sigmay(),qeye(2),sigmay())
    #ISz8=tensor(sigmaz(),sigmaz(),qeye(2))
    ISz8=tensor(sigmaz(),qeye(2),sigmaz())
    ISx=Qobj(pad(ISx8.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0)))
    ISy=Qobj(pad(ISy8.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0)))
    ISz=Qobj(pad(ISz8.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0)))

    # Geomagnetic field angle (in radian)
    theta_list=linspace(0.0, pi/2, 25)
    #theta_list=linspace(pi/2, pi, 1)
    noa=len(theta_list)

# Defining Zeeman Field Window
    #B0=[0,10**-6,20**-6,30**-6,40**-6,50**-6,60**-6,70**-6,80**-6,90**-6,100**-6]
    #B=logspace(log10(float(0.01)),log10(100),100)
    B=logspace(log10(float(0.01)),log10(1000),50)
    B0=[x*10**-6 for x in B]
    loB=len(B0)
    print('B values',loB)
    # Magnetic field made zero for calculating coherence
    #B0_coh=0


# Hyperfine Coupling Constants (in eV)
    B_geo=47*10**-6
    print('Geomagnetic field = ', B_geo)
    #a=logspace(log10(10**-2),log10(10**2),200)
    #a=logspace(log10(5),log10(10**2),20)
    #a=10
    #ax=gama*B_geo*4
    #ay=gama*B_geo*3
    #az=0
    #a=logspace(log10(float(1.00)),log10(20),1)
    #az=[x*gama*B0 for x in a]
    #loa=len(az)
    #print('No of hyperfines = ', loa)

    #ax=0.100*5*10**-9
    #ay=0.100*5*10**-9
    #az=1.3*10**-8
    #ax=ax*0.1519756*10**16 # in per-sec
    #ay=ay*0.1519756*10**16 # in per-sec
    #az=az*0.1519756*10**16 # in per-sec
    ax=0.3455 # In Gauss
    ay=0.3455 # In Gauss
    az=8.5 # In Gauss

    ax=ax*10**-4 # in Tesla
    ay=ay*10**-4 # in Tesla
    az=az*10**-4 # in Tesla

    ax=ax*muB # in eV
    ay=ay*muB # in eV
    az=az*muB # in eV

    ax=ax*0.1519756*10**16 # in per-sec
    ay=ay*0.1519756*10**16 # in per-sec
    az=az*0.1519756*10**16 # in per-sec

    # Decay rates [These three rates are simulated]
    #k=5*10**4
    ks=2*10**4
    kt=2*10**4
    k=2*10**4
    #lok=len(k)
    print('ks = ', ks)
    print('kt = ', kt)

    # Time steps in Integration. It would take 6/k time to decay the population by almost 99%
    tlist=[]
    tlist = linspace(0.0, float(6)/k, 3000) # Defining time instants for simulation
    noe=len(tlist)
    dt=tlist[1]-tlist[0]
    print("dt = ", dt)

    # Defining Projection operators    
   
    P=[[0 for x in range(10)] for x in range(10)] # 10x10 array containing only 0s

    P1=[[0 for x in range(10)] for x in range(10)] # 10x10 array containing only 0s
    P1[8][1]=1/sqrt2
    P1[8][2]=-1/sqrt2
    P1=sqrt(ks)*Qobj(P1)
    
    P2=[[0 for x in range(10)] for x in range(10)] # 10x10 array containing only 0s
    P2[8][5]=1/sqrt2
    P2[8][6]=-1/sqrt2
    P2=sqrt(ks)*Qobj(P2)
   
    P3=[[0 for x in range(10)] for x in range(10)] # 10x10 array containing only 0s
    P3[9][1]=1/sqrt2
    P3[9][2]=1/sqrt2
    P3=sqrt(kt)*Qobj(P3)

    P4=[[0 for x in range(10)] for x in range(10)] # 10x10 array containing only 0s
    P4[9][5]=1/sqrt2
    P4[9][6]=1/sqrt2
    P4=sqrt(kt)*Qobj(P4)

    P5=[[0 for x in range(10)] for x in range(10)] # 10x10 array containing only 0s
    P5[9][0]=1
    P5=sqrt(kt)*Qobj(P5)

    P6=[[0 for x in range(10)] for x in range(10)] # 10x10 array containing only 0s
    P6[9][3]=1
    P6=sqrt(kt)*Qobj(P6)

    P7=[[0 for x in range(10)] for x in range(10)] # 10x10 array containing only 0s
    P7[9][4]=1
    P7=sqrt(kt)*Qobj(P7)
	
    P8=[[0 for x in range(10)] for x in range(10)] # 10x10 array containing only 0s
    P8[9][7]=1
    P8=sqrt(kt)*Qobj(P8)

    # Defining Collapse Operators
    c_op_list = [P1,P2,P3,P4,P5,P6,P7,P8]

    # Variables
    sensitivity=[]
    coherence=[]

    for lb in range(0,loB):  # Loop over various hyperfine constant values
        singlet_yield=[] # Storing Singlet yield for
        singlet_yield_coh=[] # Defining Singlet Yield variable
        triplet_yield=[] # Defining Triplet Yield variable
        Cohr=0
        for th in range(0,noa): # Loop over various angle values
            Bx=B0[lb]*sin(theta_list[th])
            By=0
            Bz=B0[lb]*cos(theta_list[th])
            #Brfx = Brf*sin(3*pi/2+theta_list[th])
            #Brfy = 0
            #Brfz = Brf*cos(3*pi/2+theta_list[th])

            H0 = 0.5*gama*Bx*Sx + 0.5*gama*By*Sy + 0.5*gama*Bz*Sz + 0.25*ax*ISx + 0.25*ay*ISy + 0.25*az*ISz
            #H0 = gama*Bx*Sx/2 + gama*By*Sy/2 + gama*Bz*Sz/2 + ax*ISx/4 + ay*ISy/4 + az*ISz/4
            #H0 = Bx*Sx/2 + By*Sy/2 + Bz*Sz/2 + ax*ISx/4 + ay*ISy/4 + az*ISz/4

            #H1 = gama*Brfx*Sx + gama*Brfy*Sy + gama*Brfz*Sz
            #args={'w':w}
            #H=[H0,[H1,'cos(w*t)']]
                # Running the Master Equation solver
            opts = Odeoptions(nsteps=5000)

            result = mesolve(H0, rho0, tlist, c_op_list, [],options=opts)
            #result_coh = mesolve(H0_coh, rho0_offd, tlist, c_op_list, [],options=opts)
            #print('theta = ', theta_list[th],"Done")
            #result = mesolve(H, rho0, tlist, c_op_list, [], args=args, options=None) # Hamiltonian with RF field also.
            #rh_ss=steadystate(H, c_op_list,maxiter=500, tol=1e-06, method='mesolve')

                # Output states of the differential equation solver
            rh=[]
            rh_coh=[]
            rh_static=[]
            rh=result.states
            singlet_yield.append(abs(rh[noe-1][8][0][8]))
            #singlet_yield_coh.append(abs(rh_coh[noe-1][8][0][8]))
            triplet_yield.append(abs(rh[noe-1][9][0][9]))

        print('B0 = ', B0[lb], 'Done!')
        Sens=max(singlet_yield)-min(singlet_yield)
        print('Sensitivity = ', Sens)
        sensitivity.append(Sens)
        #Cohr=max(singlet_yield_coh)
        #coherence.append(Cohr)

    # Saving data to a file
    ar_data=array([B,sensitivity])
    fl = open('data/FW_ax0.3455G_ay0.3455G_az8.5G_k2104.csv', 'w')
    writer = csv.writer(fl)
    writer.writerow(['Magnetic Field (microTesla)', 'Sensitivity']) #if needed
    for values in ar_data:
        writer.writerow(values)
    fl.close()

    #Plotting Sensitivity/Singlet Yield
    figure()
    #ylim([0.00,0.20])
    semilogx(B,sensitivity,'r')
    #plot(B,sensitivity,'r')
    xlabel("B ($\mu$T)",fontsize=18)
    ylabel('Sensitivity ($D_S$)', fontsize=18 )
    savefig('Figures/FW_ax0.3455G_ay0.3455G_az8.5G_k2104.eps')
    savefig('Figures/FW_ax0.3455G_ay0.3455G_az8.5G_k2104.pdf')

    figure()
    plot(B,sensitivity,'r')
    xlabel("B ($\mu$T)",fontsize=18)
    ylabel('Sensitivity ($D_S$)', fontsize=18 )
    savefig('Figures/FW_ax0.3455G_ay0.3455G_az8.5G_k2104_Linear.eps')
    savefig('Figures/FW_ax0.3455G_ay0.3455G_az8.5G_k2104_Linear.pdf')

    show()
if __name__ == '__main__':
    run()
