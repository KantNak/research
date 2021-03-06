#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt

## parameter
g = 9.80665
Ms0 = 3017.0
Mq0 = 2431.0
Msp = 0.0435
Mqp = 0.0507
E0 = 125
R0 = 150
a = 28
r = 0.2
cq = 0.009
gms = 3.13
bs = 1.0
cr = 0.0038
a1 = 2.32
a2 = 10.8
pi = np.pi
dx = 0.01

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

Rclr = R0 + cr*T
omg = -g*Rclr / Ms
base = np.array([T]) ## value at fc
bases = []
bases.append(base)
x = fc
xnum = 1
print(fc)

while x < 1:
    #cp
    k1 = (omg*Mqp - g*cq)*base / (omg*Mqp*(1-x))
    k2 = g*(E0 + a*dsst*np.cos(pi*x))/((1-x)*omg*Mqp)
    k3 = Mq0 / (Mqp*(1-x))
    ksum = k1 + k2 + k3
    
    base = base + dx * ksum
    print(base)
    bases.append(base)

    #print(x)
    x = x + dx
    xnum = xnum + 1


#base = (g*(E0+a*dsst) - omg*Mq0) / (omg*Mqp)
bases.append(base)
xnum = xnum + 1
bases = np.array(bases)
xnum = np.arange(fc, x+dx, dx)
#print(fig.shape)
print(bases.shape)
print(xnum.shape)
plt.plot(xnum,bases)
plt.show()
