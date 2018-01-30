#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## for larger SST

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
a1 = 2.32
a2 = 10.8
pi = np.pi
dx = 0.001
xmax = 1000

## changable parameter
dsst = 2.0
T = 3100.0

## varible
q = 0.0
fc = 0.0

T = 3100.0
Ms = Ms0 + Msp*T
Mq = Mq0 + Mqp*T
C = Mq * (1+r) / ((Ms-Mq) - r*Mq)
fc3 = 3*E0 / ((pi**2) * a * dsst * (1+C))
fc = np.cbrt(fc3)

out = np.array([Mq]*xmax)
Rclr = R0 + cr*T
omg = -((g*Rclr) / Ms)
Mq = (-g*(E0-a*dsst)) / omg
ix = xmax - 1
out[ix] = Mq
x = 1.0 - dx

while x > fc:
    cpx = np.cos(pi*x)
    k1 = g*(E0 + a*dsst*cpx)/((1-x)*omg) + out[ix]/(1-x)
    
    out[ix-1] = out[ix] - dx * k1

    x = x - dx
    ix = ix - 1

xnum = np.arange(0, 1, dx)
#plt.plot(xnum[ix:xmax-1], out[ix:xmax-1])
#plt.plot(xnum[0:ix-1], out[0:ix-1])
plt.plot(xnum, out)
plt.xlim(0.0,1.0)
plt.title("distribution of Mq")
plt.xlabel("x")
plt.ylabel("Mq [ J/kg ]")
plt.show()
