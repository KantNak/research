#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
import random
import sys

dt = 12.0/2000.0
tstep = 82
rho = 28.0
sigma = 10.0
beta = 8.0/3.0
tem1 = 0
tem2 = 0
fnt = np.arange(tstep)
tstate1 = np.array([1, 2, 42])
tstate2 = np.array([1, 2, 9])
tstate3 = np.array([1, -1, 11])

## 描画位置設定
def subplot():
    fig = plt.figure()
    ax1 = fig.add_subplot(1,3,1)
    ax2 = fig.add_subplot(1,3,2)
    ax3 = fig.add_subplot(1,3,3)
    return (ax1, ax2, ax3)

## Lorenz 計算中身
def L63core(state):
    x, y, z = state
    return np.array([sigma*(y-x), x*(rho-z)-y, x*y-beta*z])

## 数値計算用ループ
def L63(cstate, em1, em2,  axp, axis):
    ## アンサンブルメンバー作成
    state = cstate
    em = np.array([0,em1,em2])
    state = state + em
    states = []
    
    for n in range(tstep):
        s_tmp = state
        k1 = L63core(s_tmp)  ## Lorenz 計算呼び出し
        s_tmp = state + 0.5*dt*k1
        k2 = L63core(s_tmp)
        s_tmp = state + 0.5*dt*k2
        k3 = L63core(s_tmp)
        s_tmp = state + dt*k3
        k4 = L63core(s_tmp)
        
        state = state + (k1+2.0*k2+2.0*k3+k4)*(dt/6.0)
        states.append(state)

    ## 配列の整理(計算結果)
    states = np.array(states)
    
    if (em1 == 0 and em2 == 0):  ## 真値？のプロット
        axp.plot(fnt, states[:,axis], color='green', linewidth = 2.5, alpha = 1)
    else:  ## アンサンブルメンバーのプロット
        axp.plot(fnt, states[:,axis], color='black', linewidth = 0.5, alpha = 0.6)
        return states

## アンサンブル用ループ
def L63ens(base_state, axp, axis):
    ens_states = 0
    for i in range(emn):
        cstate = base_state
        while True:  ## 円内ならループ抜け出す
            cem1 = random.random() - 0.5
            cem2 = random.random() - 0.5
            r2 = np.sqrt(cem1**2 + cem2**2)
            if (r2 < 0.5):
                break

        cem1 = random.random() - 0.5  ## アンサンブルメンバー用ランダム関数
        cem2 = random.random() - 0.5

        states = L63(cstate, cem1, cem2, axp, axis)
        ens_states = ens_states + states

    ens_mean = ens_states / emn  ## アンサンブル平均
    axp.plot(fnt, ens_mean[:,axis], color='red', linewidth = 2.5, alpha = 1)

## サブルーチンの実行
for j in range (3):
    axis = 2  ## x->0, y->1, z->2
    if (j == 0):  ## アンサンブルメンバー数
        emn = 10
    elif (j == 1):
        emn = 100
    elif (j == 2):
        emn = 1000
    
    ax1, ax2, ax3 = subplot()
    L63ens(tstate1, ax1, axis)
    L63(tstate1, tem1, tem2, ax1, axis)
    L63ens(tstate2, ax2, axis)
    L63(tstate2, tem1, tem2, ax2, axis)
    L63ens(tstate3, ax3, axis)
    L63(tstate3, tem1, tem2, ax3, axis)
    ax1.set_title("[1, 2, 42]")
    ax2.set_title("[1, 2, 9]")
    ax3.set_title("[1, -1, 11]")
    
plt.show()