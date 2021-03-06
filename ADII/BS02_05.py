#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## for distribution of f_c

import numpy as np
from matplotlib import pyplot as plt

## parameter
g = 9.80665
Ms0 = 3017.0
Mq0 = 2431.0
Msp = 0.0435
Mqp = 0.0507
E0 = 125.0
R0 = 150.0
a = 28.0
r = 0.2
cr = 0.0038
cq = 0.009
a1 = 2.32
a2 = 10.8
pi = np.pi
dt = 0.01
tmax = 300

## changable parameter
dsst = 0.0 + dt
T = 0.0

## varible
q = 0.0
fc = np.array([0.0]*tmax)
it = 0

while dsst <= 3.0:
    T = a*dsst / (cq + cr/(1+r))
    
    Ms = Ms0 + Msp*T
    Mq = Mq0 + Mqp*T
    C = Mq * (1+r) / ((Ms-Mq) - r*Mq)
    fc3 = 3*E0 / ((pi**2) * a * dsst * (1+C))
    fc[it] = np.cbrt(fc3)

    dsst = dsst + dt
    it = it + 1

tnum = np.arange(0, 3, dt)
#print(fig.shape)
plt.plot(tnum, fc)
plt.xlim(0.0, 3.0)
plt.ylim(0.0, 1.0)
plt.title("Convective area fraction")
plt.xlabel("$\Delta$SST [K]")
plt.ylabel("$f_c$")
plt.show()
