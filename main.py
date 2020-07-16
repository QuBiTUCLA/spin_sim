#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# a script for Radical-Pair simulations
# author's email: reza8@ucla.edu
#--------------------------------------------

#log consts
LOGLEVEL={0:"DEBUG",1:"INFO",2:"WARN",3:"ERR",4:"FATAL"}
#log function. I believe my log function is more friendly than logging package
import sys,traceback,math
def log(msg,l=1,end="\n",logfile=None):
    st=traceback.extract_stack()[-2]
    lstr=LOGLEVEL[l]
    now_str="%s %03d"%(time.strftime("%y/%m/%d %H:%M:%S",time.localtime()),math.modf(time.time())[0]*1000)
    if l<3:
        tempstr="%s<%s:%d,%s> %s%s"%(now_str,st.name,st.lineno,lstr,str(msg),end)
    else:
        tempstr="%s<%s:%d,%s> %s:\n%s%s"%(now_str,st.name,st.lineno,lstr,str(msg),traceback.format_exc(limit=2),end)
    print(tempstr,end="")
    if l>=1:
        if logfile==None:
            logfile=sys.argv[0].split(".")
            logfile[-1]="log"
            logfile=".".join(logfile)
        with open(logfile,"a") as f:
            f.write(tempstr)


import matplotlib.pyplot as plt
import numpy as np
from IPython.display import Image
from qutip import *

def H_10st(a, omega):
    """The Hamiltonian for the theory developed in paper by 
    Timmel et al.: https://doi.org/10.1080/00268979809483134"""
    
    h_hyp = a[0]*(tensor(tensor(sigmax(), qeye(2)), sigmax())) + a[1]*(tensor(tensor(sigmay(), qeye(2)), sigmay())) + \
        a[2]*(tensor(tensor(sigmaz(), qeye(2)), sigmaz()))

    h_zeem = omega[0]*(tensor(tensor(qeye(2), sigmax()), qeye(2)) + tensor(tensor(qeye(2), qeye(2)), sigmax())) + \
        omega[1]*(tensor(tensor(qeye(2), sigmay()), qeye(2)) + tensor(tensor(qeye(2), qeye(2)), sigmay())) + \
            omega[2]*(tensor(tensor(qeye(2), sigmaz()), qeye(2)) + tensor(tensor(qeye(2), qeye(2)), sigmaz()))
    
    H = h_hyp + h_zeem
    return H




if __name__=="__main__":
    print("This is main.py\n")

    #Parameters of Hamiltonian
    a = [1., 1., 1.]
    omega = [1., 1., 1.]
    k = [1., 1., 1., 1., 1., 1., 1., 1.]
    
    #Initial state
    psi_0 = 1/np.sqrt(2)*(basis(4, 2) - basis(4, 1))
    rho_0 = 0.5*tensor(qeye(2), tensor(psi_0,psi_0.dag()))
    trnsf = Qobj([[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],
    [0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0]])
    #rho_01 = trnsf*0.5*tensor(qeye(2), tensor(psi_0,psi_0.dag()))*trnsf.trans()
    
    #Collapse operators
    #p1 = Qobj([[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],
    #[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,1,-1,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0]])
    p2 = Qobj([[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,1,-1,0,0,0],[0,0,0,0,0,0,0,0,0,0]])
    p3 = Qobj([[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,1,1,0,0,0,0,0,0,0]])
    p4 = Qobj([[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,1,1,0,0,0]])
    p5 = Qobj([[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0]])
    p6 = Qobj([[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0]])
    p7 = Qobj([[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0]])
    p8 = Qobj([[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1,0,0]])

    #p=[p1,p2,p3,p4,p5,p6,p7,p8]
    p1 = Qobj([[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,1,-1,0,0]])
    p = [p1]
#     tem = 0.
#     for i in range(8):
#         tem = tem + lindblad_dissipator(p[i])
#         #tem = tem + k[i]*tensor(tensor(p[i], rho_0), p[i].dag()) #- 0.5*k[i]*(tensor(tensor(p[i].dag(),p[i]), rho_0) + tensor(tensor(rho_0, p[i].dag()), p[i]))
#     c_opt = tem   
    
    c_opt = p
    L = liouvillian(H_10st(a, omega), c_opt)
    tlist = np.linspace(0, 10, 100)
    result = mesolve(L, rho_0, tlist, [], [])
    print(result.states[-1])
    #print(tem)