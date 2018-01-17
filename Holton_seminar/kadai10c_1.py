#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
import sys

dt = 12.0/2000.0
tstep = 90
rho = 28.0
sigma = 10.0
beta = 8.0/3.0
spt = 0
state = np.array([0, 0, 0])
states = []

def L63core(state):
    x, y, z = state
    return np.array([sigma*(y-x), x*(rho-z)-y, x*y-beta*z])

def L63(i, j, k, p):
    state = np.array([i, j, k])
    states = []

    for n in range(tstep):
        s_tmp = state
        k1 = L63core(s_tmp)
        s_tmp = state + 0.5*dt*k1
        k2 = L63core(s_tmp)
        s_tmp = state + 0.5*dt*k2
        k3 = L63core(s_tmp)
        s_tmp = state + dt*k3
        k4 = L63core(s_tmp)
        
        state = state + (k1+2.0*k2+2.0*k3+k4)*(dt/6.0)
        states.append(state)

    states = np.array(states)
    
    #    fnt = tstep
    fnt = np.arange(tstep)
    plt.subplot(1, 3, p)
    #plt.plot(states[:,0], states[:,1], color='black', linewidth = 0.8)
    plt.plot(fnt, states[:,0], color='black', linewidth = 0.8)

L63(1, 2, 42, 1)
L63(1, 2, 9, 2)
L63(1, -1, 11, 3)
plt.show()
