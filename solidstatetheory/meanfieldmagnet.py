# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 14:59:07 2022

@author: william
"""
import numpy as np
import matplotlib.pyplot as plt
import math

#Constants

J = 1
g = 2
mu_b = 1
nu = 4
k_b = 1

#Plot params
    
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['lines.linewidth'] = 2

def fixed_point_iteration(x0, f, tol = 1e-5, max_iterations = 10000):
    """
    Solve x = f(x) with fixed point iteration using x0 as initial guess.
    The iteration stops when the iteration converges (tolerance = tol)
    or until max_iterations have been reached (no convergence).
    """
    x_old = x0
    x_new = f(x0)
    for i in range(max_iterations):     #loop until max steps reached
        if abs(x_new - x_old) < tol:    #break if convergence reached
            break
        #update values
        x_old = x_new
        x_new = f(x_new)
    else:
        raise Exception('Max steps (' + str(max_iterations) + ') reached')
    
    return x_new

def create_f(T, B, exact = True):
    """
    Create the rhs function in the the equation (m = f(m)) with a certain
    T and B. Can be with the exact equation (tanh) or with approximation.
    """
    def f_exact(x):
        return -1/2*math.tanh((g*mu_b*B - nu*J*x)/(2*k_b*T))
    def f_approx(x):
        argument = (g*mu_b*B - nu*J*x)/(2*k_b*T)
        return -1/2*(argument - argument**3/3)
    
    if exact:
        return f_exact
    else:
        return f_approx

def calc_m_vals(Bs, Ts, m0, use_exact = True):
    """
    Calculate m(B,T) for a range of values of B and T (Bs, Ts being vectors)).
    with initial guess m0. Exact/inexact rhs function (True/False)
    """
    m_vals = np.zeros(len(Ts))
    
    for i, T, B in zip(range(len(Ts)), Ts, Bs): #loop through B and T vectors
        f = create_f(T, B, use_exact)   #create rhs function for each B,T
        m_vals[i] = fixed_point_iteration(m0, f)
    
    return m_vals

def task1():
    Ts_inexact = np.linspace(0.51, 1.2)
    Ts_exact = np.linspace(0.2, 1.2)
    Bs = np.zeros(len(Ts_exact))
    m0 = 0.3
    m_vals_inexact = calc_m_vals(Bs, Ts_inexact, m0, False)
    m_vals_exact = calc_m_vals(Bs, Ts_exact, m0, True)
    
    fig1, ax1 = plt.subplots(figsize=(5,3), dpi = 300)
    ax1.plot(Ts_exact, m_vals_exact)
    ax1.plot(Ts_inexact, m_vals_inexact, '-.')
    ax1.set_title("Mean-field magnetization, B = 0")
    ax1.set_xlabel(r'$k_bT$')
    ax1.set_ylabel(r'$m$')
    ax1.legend(['Exact', 'With approximation'])
    
def task2a():
    Bs = np.linspace(-1, 1, 500)
    temps = [0.5, 0.9, 1.1, 1.5]
    m0 = 0.2
    fig1, ax1 = plt.subplots(figsize=(5,3), dpi = 300)
    leg = []
    for temp in temps:
        temp_vec = temp*np.ones(len(Bs))
        m_vals = calc_m_vals(Bs, temp_vec, m0)
        ax1.plot(Bs, m_vals)
        leg.append('T = ' + str(temp))

    ax1.set_title("Finite field magnetization")
    ax1.set_xlabel(r'$B$')
    ax1.set_ylabel(r'$m$')
    ax1.legend(leg)
    
def task2b():
    temps = np.linspace(1.1, 2)
    h = 0.001
    Bs = np.array([-h,h])
    ms = np.zeros((len(temps), 2))
    m0 = 0.3
    
    for i, temp in enumerate(temps):
        Ts = temp*np.ones(2)
        current_ms = calc_m_vals(Bs, Ts, m0)
        ms[i,:] = current_ms
        
    derivatives = (ms[:,1] - ms[:,0])/(2*h)
    
    
    fig1, ax1 = plt.subplots(figsize=(5,3), dpi = 300)
    ax1.plot(temps, derivatives)
    ax1.set_title("Magnetic susceptibility")
    ax1.set_xlabel(r'$T$')
    ax1.set_ylabel(r'$\chi$')
    
def task3a():
    Ts_beta = np.linspace(0.97, 0.99)
    Bs = np.zeros(len(Ts_beta))
    Tc = 1
    m0 = 0.3
    m_vals = calc_m_vals(Bs, Ts_beta, m0)
    x_beta = np.log((Tc-Ts_beta)/Tc)
    y = np.log(m_vals)
    fit_beta = np.polyfit(x_beta,y,1)
    
    fig1, ax1 = plt.subplots(figsize=(5,3), dpi = 300)
    ax1.plot(x_beta,y, '.')
    ax1.plot(x_beta, np.polyval(fit_beta,x_beta))
    ax1.set_xlabel(r'$\ln\frac{T_C-T}{T_C}$')
    ax1.set_ylabel(r'$\ln m$')
    ax1.legend(['Values from ex. 1', 'Linear fit'])
    
    
def task3b():
    temps = np.linspace(1.001, 1.005)
    h = 1e-7
    Bs = np.array([-h,h])
    ms = np.zeros((len(temps), 2))
    m0 = 0.3
    Tc = 1
    
    for i, temp in enumerate(temps):
        Ts = temp*np.ones(2)
        current_ms = calc_m_vals(Bs, Ts, m0)
        ms[i,:] = current_ms
        
    derivatives = (ms[:,1] - ms[:,0])/(2*h)
    
    x = np.log(temps-Tc)
    y = np.log(np.abs(derivatives))
    fit = np.polyfit(x,y,1)

if __name__ == '__main__':
    task3a()
    
    
    
        
        
    





