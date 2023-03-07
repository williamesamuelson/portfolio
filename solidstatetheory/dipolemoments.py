import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapz

def create_harm_osc(n, a):
    hermite_vec = np.zeros(n+1)
    hermite_vec[-1] = 1
    hermite = np.polynomial.hermite.herm2poly(hermite_vec)

    def H(x):
        hermite_val = np.polyval(hermite[::-1], x-a) 
        return 1/np.sqrt(2**n*math.factorial(n))*1/math.pi**(1/4)*math.exp(-(x-a)**2/2)*hermite_val
    return np.vectorize(H)

def create_integrand(n,m,a):
    psi = create_harm_osc(n,-a)
    phi = create_harm_osc(m, a)
    def integrand(x):
        return psi(x)*x*phi(x)
    return np.vectorize(integrand)

def calc_mu(a, limits):
    dim = 10
    x = np.linspace(-a-limits, a + limits, 1000)
    mu = np.zeros((dim,dim))
    for n in range(dim):
        for m in range(10):
             integrand = create_integrand(n,m,a/2)
             f = integrand(x)
             mu[n,m] = trapz(f,x)
    
    return mu

def bmatrix(a):
    if len(a.shape) > 2:
        raise ValueError('bmatrix can at most display two dimensions')
    lines = str(a).replace('[', '').replace(']', '').splitlines()
    rv = [r'\begin{bmatrix}']
    rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
    rv +=  [r'\end{bmatrix}']
    return '\n'.join(rv)

if __name__ == '__main__':
    a = 5
    limits = 10
    mu = calc_mu(a, limits)

    print(bmatrix(np.round(mu,2)))


