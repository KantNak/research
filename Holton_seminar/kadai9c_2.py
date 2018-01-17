#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
import sys

# パラメータ設定
pi = np.arctan(1.0)*4.0
cp = 1000.
sorg = 1367.
rho = 1.2
h = 8300.
e = 0.6
sig = 5.67e-8
dt = 2592000
dm = 1.0

# 配列作成
T = np.array([0.]*2000)
T[0] = 289.58
Tn = np.array([0.]*15)
sn = np.array([sorg]*15)

## Runge-kutta のサブルーチン
j = np.array([0.]*2000)
k = 0

def RK(i, s):
    n = 0 + i*50
    while n <= 49 + i*50:
        # アルベドの違い
        if ( T[n] - 273.15 > -10 ):
            a = 0.30
        else:
            a = 0.62

        #print(a)
        # 計算中身
        k1 = (s*(1-a))/(4*cp*rho*h)-(e*sig*(T[n])**4)/(cp*rho*h)
        k2 = (s*(1-a))/(4*cp*rho*h)-(e*sig*(T[n]+(dt*dm*k1/2))**4)/(cp*rho*h)
        k3 = (s*(1-a))/(4*cp*rho*h)-(e*sig*(T[n]+(dt*dm*k2/2))**4)/(cp*rho*h)
        k4 = (s*(1-a))/(4*cp*rho*h)-(e*sig*(T[n]+(dt*dm*k3))**4)/(cp*rho*h)
        
        T[n+1] = T[n] + (dt*dm/6)*(k1+2*(k2+k3)+k4)
        
        j[n+1] = j[n] + 1
        n = n + 1

    ### for figure
    Tn[i] = T[n]
    sn[i] = s

i = 0
k = -1
while i <= 14:
    if ( i == 0 ):
        s = sorg
        RK(i, s)
    elif ( i <= 3 ):
        s = sorg*(1.0 + k*0.1)
        print(1.0 + k*0.1)
        RK(i, s)
        k = k - 1
    elif ( i <= 10 ):
        s = sorg * (1.0 + k*0.1)
        print(1.0 + k*0.1)
        RK(i, s)
        k = k + 1
    else:
        s = sorg * (1.0 + k*0.1)
        print(1.0 + k*0.1)
        RK(i, s)
        k = k - 1
    i = i + 1

plt.plot(sn[0:4], Tn[0:4], marker="o", color="royalblue")
plt.plot(sn[3:5], Tn[3:5], color="royalblue", linestyle="dashed")
plt.plot(sn[4:11],Tn[4:11], marker="o", color="y")
plt.plot(sn[10:12], Tn[10:12], color="y", linestyle="dashed")
plt.plot(sn[11:15],Tn[11:15], marker="o", color="g")
plt.xlabel("S [${W/m^2}$]")
plt.ylabel("T [K]")
plt.text(sn[0]-55, Tn[0]+3, "S=1367", color="royalblue")
plt.text(sn[4]+45, Tn[4]-1, "0.6S", color="y")
plt.text(sn[11]-75, Tn[11], "1.3S", color="g")
plt.text(sn[0]-400, Tn[0], "ICE-FREE")
plt.text(sn[0]+200, Tn[7], "ICE-COVERD")
plt.title("Hysteresis")
plt.show()
#plt.savefig("kadai9c_2.png")
