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
base = np.array([Mq]) ## value at fc
bases = []
bases.append(base)
x = fc
xnum = 1

print(omg)

while x < 1:
    #cp
    cpx = np.cos(pi*x)
    k1 = g*(E0 + a*dsst*cpx)/((1-x)*omg) + base/(1-x)
    
    base = base + dx * k1
    print(cpx)
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
