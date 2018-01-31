#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## for small SST

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
T = 0.0

## varible
q = 0.0
fc = 0.0

Mq = Mq0 + Mqp*T
out = np.array([Mq]*xmax)
xnum = np.arange(0, 1, dx)
#print(fig.shape)
plt.plot(xnum,out)
plt.xlim(0.0, 1.0)
plt.title("distribution of $M_q$")
plt.xlabel("x")
plt.ylabel("$M_q$ [ J/kg ]")
plt.show()
