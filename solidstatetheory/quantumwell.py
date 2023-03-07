# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 21:22:50 2022

@author: willi
"""
import numpy as np
import matplotlib.pyplot as plt

def calc_eig(N, Nmiddle, Vlow, Vhigh, periodic):
    
    Nleft = (N - Nmiddle)/2

    
    V = [Vhigh if i < Nleft or i > Nleft + Nmiddle - 1 else Vlow for i in range(N)]
    T = -t * np.ones(N-1)
    H = np.diag(T,-1) + np.diag(V,0) + np.diag(T,1)

    if periodic:
        H[0,-1] = H[-1,0] = -t
    
    w, v = np.linalg.eig(H)
    sorted_indices = np.argsort(w)
    w = w[sorted_indices]
    v = v[:, sorted_indices]
    
    return w, v
    
    #plt.xlim([0,10])
if __name__=='__main__':
    Nmiddle = 40
    N = 240
    Vlow = 1
    Vhigh = 1
    t = 0.1
    energy_index = 0 #the wavefunction with the n'th lowest energy
    periodic = True #periodic/non-periodic boundary
    
    w, v = calc_eig(N, Nmiddle, Vlow, Vhigh, periodic)
    
    fig, ax = plt.subplots(figsize=(5,3),
                            dpi = 300)
    
    per = "periodic" if periodic else "non-periodic"
    
    title = 'Homogeneous potential with ' + per + ' BCs'
    #title = 'Inhomogeneous potential, ' + str(energy_index) + 'th eigenvector'
    
    #ax.set_ylim([0, 0.01])
    ax.set_xlabel('n', fontsize = 15)
    ax.set_ylabel(r'$|\psi(n)|^2$', fontsize = 15)
    ax.set_title(title, fontsize = 14)
    #plt.title(title, y=1.08, fontsize = 15)
    ax.tick_params(axis='both', labelsize=12)
    ax.plot(list(range(N)), np.abs(v[:, energy_index])**2)
    
    print(w)
