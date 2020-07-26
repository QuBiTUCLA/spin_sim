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



