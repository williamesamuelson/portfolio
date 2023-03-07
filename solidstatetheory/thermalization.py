# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 11:25:27 2022

@author: william
"""
import math as m
import numpy as np
import matplotlib.pyplot as plt


def euler_step(t, y, dt, rhs):
    """
    Takes one step in the direction of the derivative
    """
    return t + dt, y + dt*rhs(t, y)


def integrate(y0, t0, tf, dt, rhs):
    """
    Integrates the equation y'=rhs(t,y) from (t0,y0) until tf
    """
    tres = [t0]
    yres = [y0]
    while (tres[-1] < tf):
        dt = min(dt, tf-tres[-1])
        t, y = euler_step(tres[-1], yres[-1], dt, rhs)
        tres.append(t)
        yres.append(list(y))

    return tres, np.array(yres)


if __name__ == '__main__':
    tau = 2
    beta = 2
    deltaE = 1
    n_alpha = 1/(m.exp(deltaE*beta) - 1)

    def rhs_2levels(t, y):
        rhs1 = 1/tau*(1+n_alpha)*y[1]*(1-y[0]) - 1/tau*n_alpha*y[0]*(1-y[1])
        rhs2 = -rhs1
        return np.array([rhs1, rhs2])

    def rhs_3levels(t, y):
        rhs1 = 1/tau*(1+n_alpha)*y[1]*(1-y[0]) - 1/tau*n_alpha*y[0]*(1-y[1])
        rhs3 = 1/tau*n_alpha*y[1]*(1-y[2]) - 1/tau*(1+n_alpha)*y[2]*(1-y[1])
        rhs2 = -rhs1-rhs3
        return np.array([rhs1, rhs2, rhs3])

    dt = 0.01
    y0_2levels = [1/2, 1/2]
    y0_3levels = [0, 0, 1]
    t0 = 0
    tf = 17
    t2, y2 = integrate(y0_2levels, t0, tf, dt, rhs_2levels)
    t3, y3 = integrate(y0_3levels, t0, tf, dt, rhs_3levels)

    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    fig1, ax1 = plt.subplots(figsize=(5, 3), dpi=300)
    ax1.plot(t2, y2, t2, y2[:, 0] + y2[:, 1])
    ax1.set_title("Thermalization, 2-level system")
    ax1.legend([r'$f_1$', r'$f_2$', 'Total'])
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Occupation')
    plt.show()

    fig2, ax2 = plt.subplots(figsize=(5, 3), dpi=300)
    ax2.plot(t3, y3, t3, y3[:, 0]+y3[:, 1]+y3[:, 2])
    ax2.set_title("Thermalization, 3-level system")
    ax2.legend([r'$f_1$', r'$f_2$', r'$f_3$', 'Total'])
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Occupation')
    plt.show()
