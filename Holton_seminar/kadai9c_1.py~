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
dt = 2592000.

# 配列作成
T = np.array([0.]*2000)
j = np.array([0.]*2000)
T[0] = 289.58

## Runge-kutta のサブルーチン

def RK(i, s, k, dm, clr):
    n = 0 + i*50
    while n <= 49 + i*50:
        # アルベドの違い
        if ( T[n] - 273.15 > -10 ):
            a = 0.30
        else:
            a = 0.62

        # 計算中身
        k1 = (s*(1-a))/(4*cp*rho*h)-(e*sig*(T[n])**4)/(cp*rho*h)
        k2 = (s*(1-a))/(4*cp*rho*h)-(e*sig*(T[n]+(dt*dm*k1/2))**4)/(cp*rho*h)
        k3 = (s*(1-a))/(4*cp*rho*h)-(e*sig*(T[n]+(dt*dm*k2/2))**4)/(cp*rho*h)
        k4 = (s*(1-a))/(4*cp*rho*h)-(e*sig*(T[n]+(dt*dm*k3))**4)/(cp*rho*h)
        
        T[n+1] = T[n] + (dt*dm/6)*(k1+2*(k2+k3)+k4)
        
        j[n+1] = j[n] + 1
        n = n + 1

    ### for figure    
    if ( i < 15):
        txt = 1.0 + k*0.1
        plt.plot(j[0+i*50:n], T[0+i*50:n], color=clr, alpha = 0.6)
        plt.text(5+50*i, T[n]+1, "%.1f"%txt + "S")

# dmの違いによるループのためのサブルーチン
def roop(sdm, sclr):
    i = 0
    k = -1
    dm = sdm
    clr = sclr
    while i <= 15:
        if ( i == 0 ):
            k = 0
            s = sorg
            RK(i, s, k, dm, clr)
            k = -1
        elif ( i <= 3 ):
            s = sorg*(1.0 + k*0.1)
            RK(i, s, k, dm, clr)
            k = k - 1
        elif ( i <= 10 ):
            s = sorg * (1.0 + k*0.1)
            RK(i, s, k, dm, clr)
            k = k + 1
        else:
            s = sorg * (1.0 + k*0.1)
            RK(i, s, k, dm, clr)
            k = k - 1
        i = i + 1

# ループそれぞれ
roop(1.0, "b")
roop(2.0, "r")
roop(3.0, "y")

plt.xlim(0, 750)
plt.xlabel("time step")
plt.ylabel("T[K]")
plt.text(400, 290, "dt=1m,2m", color="mediumvioletred")
plt.text(600, 220, "dt=3m", color="y")
#plt.show()
plt.savefig("kadai9c_1.png")
