# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 14:17:49 2022

@author: willi
"""
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.titlesize'] = 14
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['lines.linewidth'] = 2

class Integrator():
    
    def __init__(self, problem):
        self.prob = problem
    
    def one_interval(self, x):
        raise NotImplementedError()
    
    def calc(self):
        res = 0
        for i in range(self.prob.steps - self.n):
            res += self.one_interval(i)
        return res
    
class Trapezoidal(Integrator):
    
    def __init__(self, problem):
        super().__init__(problem)
        self.n = 1
    
    def one_interval(self,i):
        return (self.prob.y[i] + self.prob.y[i+1])*self.prob.h/2
    
class Simpson(Integrator):
    def __init__(self, problem):
        super().__init__(problem)
        self.n = 2
    
    def one_interval(self, i):
        return (self.prob.y[i] + 4*self.prob.y[i+1] + self.prob.y[i+2])*self.prob.h/3
    
    
class Booles(Integrator):
    def __init__(self, problem):
        super().__init__(problem)
        self.n = 4
    
    def one_interval(self, i):
        return (7*self.prob.y[i] + 32*self.prob.y[i+1] + 12*self.prob.y[i+2] + 
                32*self.prob.y[i+3] + 7*self.prob.y[i+4])*self.prob.h*2/45    
    
    
class Problem():
    
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.h = abs(x[1]-x[0])
        self.steps = len(x)
        
        
def task2(a, b, omega, omega0s):
    result_a_trap = np.zeros(omega0s.shape)
    result_a_simp = np.zeros(omega0s.shape)
    result_a_boole = np.zeros(omega0s.shape)
    result_b_trap = np.zeros(omega0s.shape)
    result_b_simp = np.zeros(omega0s.shape)
    result_b_boole = np.zeros(omega0s.shape)
    
    #results = [result_a_trap, result_a_simp, result_a_boole, 
      #         result_b_trap, result_b_simp, result_b_boole]
    
    
    
    for i, omega0 in enumerate(omega0s):
        y_a = a/(omega - omega0)
        y_b = b/(omega - omega0)
        problem_a = Problem(omega, y_a)
        problem_b = Problem(omega, y_b)
        trap_a = Trapezoidal(problem_a)
        simpson_a = Simpson(problem_a)
        boole_a = Booles(problem_a)
        trap_b = Trapezoidal(problem_b)
        simpson_b = Simpson(problem_b)
        boole_b = Booles(problem_b)
        result_a_trap[i] = trap_a.calc()
        result_a_simp[i] = simpson_a.calc()
        result_a_boole[i] = boole_a.calc()
        result_b_trap[i] = trap_b.calc()
        result_b_simp[i] = simpson_b.calc()
        result_b_boole[i] = boole_b.calc()
        
    
    fig1, ax1 = plt.subplots(figsize=(5,3), dpi = 300)
    ax1.plot(omega0s, -1/np.pi*result_a_trap, omega0s, 1/np.pi*result_b_trap)
    ax1.set_title("Imaginary part of susceptibility calc. w/ Boole's rule")
    ax1.set_xlabel(r'$\omega$')
    ax1.set_ylabel(r'$Im \chi$')
    ax1.legend(['a', 'b'])
        
if __name__ == '__main__':
    file = 'kramerskronig.txt'
        
    with open(file) as f:
        content = f.readlines()
    lis = []

    for line in content:
        lis += line.split()
        
    omegas = np.array([float(om) for om in lis[0::3]])
    a = np.array([float(aa) for aa in lis[1::3]])
    b = np.array([float(bb) for bb in lis[2::3]])
    
    omega0s = omegas + 0.006
    
    im_a = np.zeros(omega0s.shape)
    im_b = np.zeros(omega0s.shape)
    
    for i, omega0 in enumerate(omega0s):
        y_a = a/(omegas - omega0)
        y_b = b/(omegas - omega0)
        problem_a = Problem(omegas, y_a)
        problem_b = Problem(omegas, y_b)
        bool_a = Booles(problem_a)
        bool_b = Booles(problem_b)
        im_a[i] = bool_a.calc()
        im_b[i] = bool_b.calc()
        
    im_a *= -1/np.pi
    im_b *= -1/np.pi
        
    re_a = np.zeros(omega0s.shape)
    re_b = np.zeros(omega0s.shape)
    
    for i, omega in enumerate(omegas):
        y_a = im_a/(omega0s-omega)
        y_b = im_b/(omega0s-omega)
        problem_a = Problem(omega0s, y_a)
        problem_b = Problem(omega0s, y_b)
        bool_a = Booles(problem_a)
        bool_b = Booles(problem_b)
        re_a[i] = bool_a.calc()
        re_b[i] = bool_b.calc()
        
    re_a *= 1/np.pi
    re_b *= 1/np.pi
        
    plt.plot(omegas, re_a, omegas, a)
    plt.show()
