#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# a script for Radical-Pair simulations
# author's email: reza8@ucla.edu
# --------------------------------------------
import traceback, math
import numpy as np
from qutip import *
import time

# log consts

LOGLEVEL = {0: "DEBUG", 1: "INFO", 2: "WARN", 3: "ERR", 4: "FATAL"}


# log function. I believe my log function is more friendly than logging package

def log(msg, l=1, end="\n", logfile=None):
    st = traceback.extract_stack()[-2]
    lstr = LOGLEVEL[l]
    now_str = "%s %03d" % (time.strftime("%y/%m/%d %H:%M:%S", time.localtime()), math.modf(time.time())[0] * 1000)
    if l < 3:
        tempstr = "%s<%s:%d,%s> %s%s" % (now_str, st.name, st.lineno, lstr, str(msg), end)
    else:
        tempstr = "%s<%s:%d,%s> %s:\n%s%s" % (now_str, st.name, st.lineno, lstr, str(msg), traceback.format_exc(limit=2), end)
    print(tempstr, end="")
    if l >= 1:
        if logfile is None:
            logfile = sys.argv[0].split(".")
            logfile[-1] = "log"
            logfile = ".".join(logfile)
        with open(logfile, "a") as f:
            f.write(tempstr)


def h_10st(a, b0, theta, phi):

    """The Hamiltonian for the theory developed in paper by
    Timmel et al.: https://doi.org/10.1080/00268979809483134"""

    # geomagnetic filed in micro-tesla
    b = b0 * np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])

    h_hyp = a[0] * (tensor(tensor(sigmax(), qeye(2)), sigmax())) + a[1] * (tensor(tensor(sigmay(), qeye(2)), sigmay())) + \
            a[2] * (tensor(tensor(sigmaz(), qeye(2)), sigmaz()))

    h_zeem = b[0] * (tensor(tensor(qeye(2), sigmax()), qeye(2)) + tensor(tensor(qeye(2), qeye(2)), sigmax())) + \
             b[1] * (tensor(tensor(qeye(2), sigmay()), qeye(2)) + tensor(tensor(qeye(2), qeye(2)), sigmay())) + \
             b[2] * (tensor(tensor(qeye(2), sigmaz()), qeye(2)) + tensor(tensor(qeye(2), qeye(2)), sigmaz()))

    h = h_hyp + h_zeem
    return h


### INSERT Brf
# Show: in tumbling, still sensitive to |B|

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
# from MakePdf import *
import pickle as pickle

import sys

# sys.path.append("/usr/bin") # necessary for the tex fonts
# sys.path.append("../Python modules/") # necessary for the tex fonts

fsizetit = 24
fsizepl = 20
sizex = 8
sizey = 6
fsizenb = 16
dpi_no = 80
lw = 2

number_of_states = 12

gyros = 28.0e9  # Hz/T

mu0h = 9.274e-24 / 6.6e-34  # mu0 / h, 1/(T*s)

#################################

# SPIN OPS - NO CHANGE NEEDED
#################################

S1_op = np.array(
    [tensor(sigmax(), qeye(2), qeye(2)), tensor(sigmay(), qeye(2), qeye(2)), tensor(sigmaz(), qeye(2), qeye(2))])
S2_op = np.array(
    [tensor(qeye(2), sigmax(), qeye(2)), tensor(qeye(2), sigmay(), qeye(2)), tensor(qeye(2), sigmaz(), qeye(2))])
nucleus = np.array(
    [tensor(qeye(2), qeye(2), sigmax()), tensor(qeye(2), qeye(2), sigmay()), tensor(qeye(2), qeye(2), sigmaz())])


# EVOLUTION FUNCTION - NO CHANGE NEEDED
#################################

def calculate_bird(rho0_vectors, B0, theta, phi, ax, ay, az, ks, kt):
    Bfield = gyros * B0 * np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])

    S2_mod_op = np.array([ax * tensor(qeye(2), sigmax(), qeye(2)), ay * tensor(qeye(2), sigmay(), qeye(2)),
                          az * tensor(qeye(2), sigmaz(), qeye(2))])

    H1 = np.dot(Bfield, 1.0 / 2.0 * (S1_op + S2_op))

    H2 = np.dot(1.0 / 2.0 * nucleus, 1.0 / 2.0 * S2_mod_op)

    # complete Hamiltonian
    H = H1 + H2

    # augmenting H with the 4 shelving states
    hlp = H.full()
    hlp2 = np.zeros([number_of_states, number_of_states], dtype=np.complex)
    hlp2[0:8, 0:8] = hlp
    H = Qobj(hlp2)

    # augmenting rho0 with the 4 shelving states
    hlp = rho0_vectors.full()
    hlp2 = np.zeros([12, 12], dtype=np.complex)
    hlp2[0:8, 0:8] = hlp
    rho0 = Qobj(hlp2)

    # coupling to bath or decay operators
    c_op_list = []

    hlp = number_of_states
    P1 = 1.0 / np.sqrt(2.0) * basis(hlp, 8) * (basis(hlp, 2).dag() - basis(hlp, 4).dag())
    P2 = 1.0 / np.sqrt(2.0) * basis(hlp, 8) * (basis(hlp, 3).dag() - basis(hlp, 5).dag())
    P3 = 1.0 / np.sqrt(2.0) * basis(hlp, 10) * (basis(hlp, 2).dag() + basis(hlp, 4).dag())
    P4 = 1.0 / np.sqrt(2.0) * basis(hlp, 10) * (basis(hlp, 3).dag() + basis(hlp, 5).dag())
    P5 = basis(hlp, 11) * basis(hlp, 0).dag()
    P6 = basis(hlp, 11) * basis(hlp, 1).dag()
    P7 = basis(hlp, 9) * basis(hlp, 6).dag()
    P8 = basis(hlp, 9) * basis(hlp, 7).dag()

    c_op_list.append(np.sqrt(ks) * P1)
    c_op_list.append(np.sqrt(ks) * P2)
    c_op_list.append(np.sqrt(kt) * P3)
    c_op_list.append(np.sqrt(kt) * P4)
    c_op_list.append(np.sqrt(kt) * P5)
    c_op_list.append(np.sqrt(kt) * P6)
    c_op_list.append(np.sqrt(kt) * P7)
    c_op_list.append(np.sqrt(kt) * P8)

    # full time evolution
    try:
        result = mesolve(H, rho0, times, c_op_list, [])
        states = result.states

        ## steady state -> does not work
        # rho_ss = steadystate(H, c_op_list)
        # print(rho_ss)

        # print(states.full())
        pops = np.zeros([len(states), number_of_states])  # time points, no states

        helper = 0
        for k in range(len(states)):

            for n in range(number_of_states):
                pops[k, n] = np.real(states[k][n, n])

            helper = helper + ks * np.exp(-ks * times[k]) * states[k].matrix_element(basis(hlp, 8),
                                                                                     basis(hlp, 8)) * timedetail

        return pops[-1, 8], 2 * np.real(helper), pops
        # Singlet yield using last point in time
        # Singlet yield = int_0^inf k exp(-kt) * <S|rho(t)|S> -> NOT YET WORKING
        #### not sure why I need this factor of 2 to make formulas match

    except:
        print('excepted')
        return np.nan, np.nan, np.nan * np.ones([int(notimepoints), int(number_of_states)])


#################################
# REAL START
#################################

def Do_colorplot(rho0_vectors, rho0_description, B_start, B_stop, B_steps, theta_start, theta_stop, theta_steps):
    do_log_bfield_scan = True
    # define the arrays to calculate the B-field and the theta angles
    if do_log_bfield_scan:
        B0_arr = np.logspace(np.log10(B_start), np.log10(B_stop), B_steps)  # in T, 100 is usual
    else:
        B0_arr = np.linspace(B_start, B_stop, B_steps)  # in T, 100 is usual

    theta_arr = np.linspace(theta_start, theta_stop, theta_steps)  # rad, 100 is usual

    result_S = []
    result_S_formula = []
    # pops_tot = np.zeros([len(B0_arr), len(theta_arr), len(times), number_of_states])
    pops_tot = np.zeros([len(B0_arr), len(theta_arr), number_of_states])

    l0 = 0

    for B0 in B0_arr:

        l1 = 0

        print("B-field = " + str(B0))

        for theta in theta_arr:
            print("   Theta = " + str(theta))

            SYend, SYform, pops = calculate_bird(rho0_vectors, B0, theta, phi, ax, ay, az, ks, kt)

            result_S.append(SYend)
            result_S_formula.append(SYform)
            # pops_tot[l0,l1,:,:] = pops
            pops_tot[l0, l1, :] = pops[-1, :]  # only save the last time point

            l1 += 1

        l0 += 1

    pickle.dump([result_S, result_S_formula, pops_tot, B0_arr, theta_arr, phi, ax, ay, az, ks, kt, tfinal, timedetail,
                 rho0_vectors], open("Data_" + rho0_description + ".p", "wb"))


###### START OF CODE

# run: python cluster_do_colorplot.py <unique_number> <tfinal> <notimepoints> <B0_start> <B0_stop> <B_steps> <theta_start> <theta_stop> <theta_steps>
#
# theta is in units of pi
#
# for example: python cluster_do_colorplot.py 001 0.5e-3 500 0 1000e-6 100 0 1 100

# transform command line arguments into the parameters
unique_number = np.float(sys.argv[1])
tfinal = np.float(sys.argv[2])
notimepoints = np.float(sys.argv[3])

B_start = np.float(sys.argv[4])
B_stop = np.float(sys.argv[5])
B_steps = np.float(sys.argv[6])

theta_start = np.float(sys.argv[7]) * np.pi
theta_stop = np.float(sys.argv[8]) * np.pi
theta_steps = np.float(sys.argv[9])

# B0 = 47e-6 #T, B field of Earth at Frankfurt
# theta = np.pi/4.0 #inclination of B field wrt to axis of RP, [0,pi]

phi = 0 * np.pi / 1.0  # azimuth of B field, , [0,2*pi]

ax = mu0h * 0.345e-4  # 0.345G given by PRE ; ax is now in Hz
ay = mu0h * 0.345e-4  # 0.345G given by PRE ; ay is now in Hz
az = mu0h * 9e-4  # 9G given by PRE ; az is now in Hz
# "Biologically feasible regime of HF strengths, 0.1 to 10G, see ref 26 in PRE"

ks = 2e4  # Hz, rate of decay into |S>
kt = 2e4  # Hz, rate of decay into |T->, |T0>, |T+>

# time array to simulate
tinit = 0.0
# tfinal = 1.0e-3
# notimepoints = 1000.0

times = np.linspace(tinit, tfinal, notimepoints)  # in seconds
timedetail = (tfinal - tinit) / notimepoints
# time zero is creation of RP

rho0_vectors = []
rho0_description = []

# 0 electrons depolarized
rho0_vectors.append(1.0 / 8.0 * tensor(qeye(4), qeye(2)))
rho0_description.append('Elecs. dep., nuc. dep.')

# 1 electrons singlet
psi0 = 1.0 / np.sqrt(4.0) * (tensor(basis(2, 0), basis(2, 1)) - tensor(basis(2, 1), basis(2, 0)))
rho0_vectors.append(tensor(tensor(psi0, psi0.dag()), qeye(2)))
rho0_description.append('S, nuc. dep.')

# 2 electrons T-, NORM
psi0 = 1.0 / np.sqrt(2.0) * (tensor(basis(2, 1), basis(2, 1)))
rho0_vectors.append(tensor(tensor(psi0, psi0.dag()), qeye(2)))
rho0_description.append('Tminus, nuc. dep.')

# 3 electrons T+, NORM
psi0 = 1.0 / np.sqrt(2.0) * (tensor(basis(2, 0), basis(2, 0)))
rho0_vectors.append(tensor(tensor(psi0, psi0.dag()), qeye(2)))
rho0_description.append('Tplus, nuc. dep.')

# 4 electrons T0, NORM
psi0 = 1.0 / np.sqrt(4.0) * (tensor(basis(2, 0), basis(2, 1)) + tensor(basis(2, 1), basis(2, 0)))
rho0_vectors.append(tensor(tensor(psi0, psi0.dag()), qeye(2)))
rho0_description.append('Tzero, nuc. dep.')

# calculate and then save
for k in np.arange(len(rho0_vectors)):
    Do_colorplot(rho0_vectors[k], rho0_description[k] + "_" + str(unique_number), B_start, B_stop, B_steps, theta_start,
                 theta_stop, theta_steps)

print("Done.")


