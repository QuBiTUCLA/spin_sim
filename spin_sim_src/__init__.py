#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# a script for Radical-Pair simulations
# author's email: reza8@ucla.edu
# --------------------------------------------
import traceback, math
import numpy as np
from qutip import *
import time
import sys


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


def build_Hamiltonian(a, b0, theta, phi):

    """The Hamiltonian for the theory developed in paper by
    Timmel et al.: https://doi.org/10.1080/00268979809483134"""

    # Spin operators
    # S1_op = np.array(
    #     [tensor(sigmax(), qeye(2), qeye(2)), tensor(sigmay(), qeye(2), qeye(2)), tensor(sigmaz(), qeye(2), qeye(2))])
    # S2_op = np.array(
    #     [tensor(qeye(2), sigmax(), qeye(2)), tensor(qeye(2), sigmay(), qeye(2)), tensor(qeye(2), sigmaz(), qeye(2))])
    # nucl = np.array(
    #     [tensor(qeye(2), qeye(2), sigmax()), tensor(qeye(2), qeye(2), sigmay()), tensor(qeye(2), qeye(2), sigmaz())])

    # geomagnetic field
    gyros = 28.0e9  # Hz/T
    b = gyros * b0 * np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])

    # S2_mod_op = np.array([a[0] * tensor(qeye(2), sigmax(), qeye(2)), a[1] * tensor(qeye(2), sigmay(), qeye(2)),
    #                       a[2] * tensor(qeye(2), sigmaz(), qeye(2))])
    #
    # h_zeem = np.dot(b, 1.0 / 2.0 * (S1_op + S2_op))
    #
    # h_hyp = np.dot(1.0 / 2.0 * nucl, 1.0 / 2.0 * S2_mod_op)

    h_hyp = a[0] * (tensor(tensor(sigmax(), qeye(2)), sigmax())) + \
            a[1] * (tensor(tensor(sigmay(), qeye(2)), sigmay())) + \
            a[2] * (tensor(tensor(sigmaz(), qeye(2)), sigmaz()))
    h_hyp = 1./4. * h_hyp
    h_zeem = b[0] * (tensor(tensor(qeye(2), sigmax()), qeye(2)) + tensor(tensor(qeye(2), qeye(2)), sigmax())) + \
             b[1] * (tensor(tensor(qeye(2), sigmay()), qeye(2)) + tensor(tensor(qeye(2), qeye(2)), sigmay())) + \
             b[2] * (tensor(tensor(qeye(2), sigmaz()), qeye(2)) + tensor(tensor(qeye(2), qeye(2)), sigmaz()))
    h_zeem = 1./2. * h_zeem
    h = h_hyp + h_zeem
    return h


def hamiltonian_add_dim(total_num_states, hamiltonian):

    """ Function for adding shelving states to the hamiltonian!"""

    # adding the shelving states as additional dimensions to the Hamiltonian
    temp_h = hamiltonian.full()
    temp_h2 = np.zeros([total_num_states, total_num_states], dtype=np.complex)
    temp_h2[0:8, 0:8] = temp_h
    ham = Qobj(temp_h2)

    return ham

def prep_ini_state():

    """Preparing the initials states and initial density matrices for various scenarios"""

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

    return rho0_vectors, rho0_description

def rho_add_dim(total_num_states, rho0_vector):

    """Function for adding extra dimensions to initial density matrices, accounting for shelving states"""

    tmp_st1 = rho0_vector.full()
    tmp_st2 = np.zeros([total_num_states, total_num_states], dtype=np.complex)
    tmp_st2[0:8, 0:8] = tmp_st1
    rho0 = Qobj(tmp_st2)

    return rho0

def time_evol(hamiltonian, rho0, num_states, times, timesteps,  timedetail, k, c_op_list, ex_op_list):

    """Time evolution, using meslove() method of QuTiP"""

    try:
        result = mesolve(hamiltonian, rho0, times, c_op_list, ex_op_list)
        states = result.states

        ## steady state -> does not work
        # rho_ss = steadystate(H, c_op_list)
        # print(rho_ss)

        # print(states.full())
        pops = np.zeros([len(states), num_states])  # time points, no states

        helper = 0
        for l in range(len(states)):

            for n in range(num_states):
                pops[l, n] = np.real(states[l][n, n])

            helper = helper + k[0] * np.exp(-k[0] * times[l]) * states[l].matrix_element(basis(num_states, 8),
                                                                                     basis(num_states, 8)) * timedetail

        return pops[-1, 8], 2 * np.real(helper), pops
        # Singlet yield using last point in time
        # Singlet yield = int_0^inf l exp(-k[1]) * <S|rho(t)|S> -> NOT YET WORKING
        #### not sure why I need this factor of 2 to make formulas match

    except:
        print('excepted')
        return np.nan, np.nan, np.nan * np.ones([int(timesteps), int(num_states)])



