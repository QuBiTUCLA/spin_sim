from spin_sim_src import *
import pickle as pickle

if __name__ == "__main__":
    print("This is main.py\n")

    # run: python <filename>.py <unique_number> <tfinal> <timesteps> <B0_start> <B0_stop> <B_steps>
    # <theta_start> <theta_stop> <theta_steps>
    #
    # theta is in units of pi
    #
    # for example: python cluster_do_colorplot.py 001 0.5e-3 500 0 1000e-6 100 0 1 100

    # transform command line arguments into the parameters
    unique_number = np.float(sys.argv[1])
    tfinal = np.float(sys.argv[2])
    timesteps = np.int(sys.argv[3])

    B_start = np.float(sys.argv[4])
    B_stop = np.float(sys.argv[5])
    B_steps = np.int(sys.argv[6])

    theta_start = np.float(sys.argv[7]) * np.pi
    theta_stop = np.float(sys.argv[8]) * np.pi
    theta_steps = np.int(sys.argv[9])

    # B0 = 47e-6 #T, B field of Earth at Frankfurt
    # theta = np.pi/4.0 #inclination of B field wrt to axis of RP, [0,pi]

    phi = 0. * np.pi / 1.0  # azimuth of B field, , [0,2*pi]

    mu0h = 9.274e-24 / 6.6e-34  # mu0 / h, 1/(T*s)
    a = [0., 0., 0.]
    a[0] = mu0h * 0.345e-4  # 0.345G given by PRE ; ax is now in Hz
    a[1] = mu0h * 0.345e-4  # 0.345G given by PRE ; ay is now in Hz
    a[2] = mu0h * 9e-4  # 9G given by PRE ; az is now in Hz
    # "Biologically feasible regime of HF strengths, 0.1 to 10G, see ref 26 in PRE"

    k = [0., 0., 0., 0.]
    k[0] = 2e4  # Hz, rate of decay into |S>
    k[1] = 2e4  # Hz, rate of decay into |T->
    k[2] = 2e4  # Hz, rate of decay into |T0>
    k[3] = 2e4  # Hz, rate of decay into |T+>

    # time array to simulate
    tinit = 0.0
    # tfinal = 1.0e-3
    # notimepoints = 1000.0

    times = np.linspace(tinit, tfinal, timesteps)  # in seconds
    timedetail = (tfinal - tinit) / timesteps
    # time zero is creation of RP

    num_states = 12

    # preparing the initial density matrices for various scenarios
    rho0_vectors, rho0_description = prep_ini_state()

    # adding dimensions to rho0 to account for the shelving states
    rho0 = rho_add_dim(total_num_states=num_states, rho0_vector=rho0_vectors[1])

    # coupling to bath or decay operators
    c_op_list = []

    P1 = 1.0 / np.sqrt(2.0) * basis(num_states, 8) * (basis(num_states, 2).dag() - basis(num_states, 4).dag())
    P2 = 1.0 / np.sqrt(2.0) * basis(num_states, 8) * (basis(num_states, 3).dag() - basis(num_states, 5).dag())
    P3 = 1.0 / np.sqrt(2.0) * basis(num_states, 10) * (basis(num_states, 2).dag() + basis(num_states, 4).dag())
    P4 = 1.0 / np.sqrt(2.0) * basis(num_states, 10) * (basis(num_states, 3).dag() + basis(num_states, 5).dag())
    P5 = basis(num_states, 11) * basis(num_states, 0).dag()
    P6 = basis(num_states, 11) * basis(num_states, 1).dag()
    P7 = basis(num_states, 9) * basis(num_states, 6).dag()
    P8 = basis(num_states, 9) * basis(num_states, 7).dag()

    c_op_list.append(np.sqrt(k[0]) * P1)
    c_op_list.append(np.sqrt(k[0]) * P2)
    c_op_list.append(np.sqrt(k[1]) * P3)
    c_op_list.append(np.sqrt(k[1]) * P4)
    c_op_list.append(np.sqrt(k[2]) * P5)
    c_op_list.append(np.sqrt(k[2]) * P6)
    c_op_list.append(np.sqrt(k[3]) * P7)
    c_op_list.append(np.sqrt(k[3]) * P8)

    ex_op_list = []

    # define the arrays to calculate the B-field and the theta angles
    do_log_bfield_scan = False
    if do_log_bfield_scan:
        B0_arr = np.logspace(np.log10(B_start), np.log10(B_stop), B_steps)  # in T, 100 is usual
    else:
        B0_arr = np.linspace(B_start, B_stop, B_steps)  # in T, 100 is usual

    theta_arr = np.linspace(theta_start, theta_stop, theta_steps)  # rad, 100 is usual

    result_S = []
    result_S_formula = []
    # pops_tot = np.zeros([len(B0_arr), len(theta_arr), len(times), number_of_states])
    pops_tot = np.zeros([len(B0_arr), len(theta_arr), num_states])

    l0 = 0

    for B0 in B0_arr:

        l1 = 0

        print("B-field = " + str(B0))

        for theta in theta_arr:
            print("   Theta = " + str(theta))

            # building the Hamiltonian
            h = build_Hamiltonian(a, B0, theta, phi)

            # adding the shelving states as additional dimensions to the Hamiltonian
            H = hamiltonian_add_dim(total_num_states=num_states, hamiltonian=h)

            SYend, SYform, pops = time_evol(H, rho0, num_states, times, timesteps, timedetail, k, c_op_list, ex_op_list)
            result_S.append(SYend)
            result_S_formula.append(SYform)
            # pops_tot[l0,l1,:,:] = pops
            pops_tot[l0, l1, :] = pops[-1, :]  # only save the last time point

            l1 += 1

        l0 += 1

    pickle.dump([result_S, result_S_formula, pops_tot, B0_arr, theta_arr, phi, a[0], a[1], a[2], k[0], k[1], tfinal,
                 timedetail, rho0], open("Data_" + rho0_description + ".p", "wb"))

    print("Done.")
