# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 15:58:19 2022

@author: willi
"""
import numpy as np
from assimulo.explicit_ode import Explicit_ODE
from assimulo.problem import Explicit_Problem
from assimulo.ode import Explicit_ODE_Exception, ID_PY_OK, NORMAL
import math
import matplotlib.pyplot as plt

class abstractBDF(Explicit_ODE):
    
    def __init__(self, problem):
        Explicit_ODE.__init__(self, problem)
        
        #Solver options
        self.options["h"] = 0.01
        
        #Statistics
        self.statistics["nsteps"] = 0
        self.statistics["nfcns"] = 0
        
        def _set_h(self,h):
                self.options["h"] = float(h)

        def _get_h(self):
            return self.options["h"]
            
        h=property(_get_h,_set_h)
        
    def integrate(self, t, y, tf, opts):
        """
        _integrates (t,y) values until t > tf
        """
        h = self.options["h"]

        
        #Lists for storing the result
        tres = []
        yres = []
        
        for i in range(self.maxsteps):
            if i == 0:
                T, Y = [t], [y]
                
            if T[-1] >= tf:
                break
            self.statistics["nsteps"] += 1
                        
            h = min(h, abs(tf-T[-1]))
            
            if self.is_at_init(i):
                t_np1, y_np1 = self.initial_steps(T, Y, h, i)
            else:   
                t_np1, y_np1 = self.step_BDF(T, Y, h)
                
           
            self.shift_index(T, Y, t_np1, y_np1, i)

            
            tres.append(T[-1])
            yres.append(Y[-1].copy())
            
            h=min(h,np.abs(tf-T[-1]))
        else:
            raise Explicit_ODE_Exception('Final time not reached within maximum number of steps')
        
        return ID_PY_OK, tres, yres
    
    def print_statistics(self, verbose=NORMAL):
        self.log_message('Final Run Statistics            : {name} \n'.format(name=self.problem.name),        verbose)
        self.log_message(' Step-length                    : {stepsize} '.format(stepsize=self.options["h"]), verbose)
        self.log_message(' Number of Steps                : '+str(self.statistics["nsteps"]),          verbose)               
        self.log_message(' Number of Function Evaluations : '+str(self.statistics["nfcns"]),         verbose)
            
        self.log_message('\nSolver options:\n',                                    verbose)
        self.log_message(' Solver type       : Fixed step',                      verbose)
    
    def is_at_init(self, i): 
        raise NotImplementedError('asjd')

    def initial_steps(self, T, Y, h, i):
        raise NotImplementedError('asjd')
    
    def step_BDF(self, T, Y, h):
        raise NotImplementedError('asjd')
    
        
    def step_EE(self, t, y, h):
        """
        This calculates the next step in the integration with explicit Euler.
        """
        self.statistics["nfcns"] += 1
        
        f = self.problem.rhs
        return t + h, y + h*f(t, y)
    
    def shift_index(self, T, Y, t_np1, y_np1, i):
        if self.is_at_init(i):
            T.append(t_np1)
            Y.append(y_np1)
        else:
            T[0:-1] = T[1::]
            T[-1] = t_np1
            Y[0:-1] = Y[1::]
            Y[-1] = y_np1
        
            
    def _jacobian(self,G, y):
        n = len(y)
        J = np.zeros((n,n))
        delta = 10**-4
        Gy = G(y)
        
        for i in range(n):
            e_i = np.zeros(n)
            e_i[i] = 1.0
            J[:,i] = (G(y+delta*e_i)-Gy)/delta
        return J
        
    def newton(self, T, Y, alpha, h, f):
        """
        Solves the non-linear equation system using Newton iteration

        """
        
        t_np1 = T[-1] + h
        y_np1_i = Y[-1] #zero order predictor
    
        
        def G(y):
           #returns alpha[0]*y + alpha[1]*y_n+k-1 + alpha[2]*y_n+k-2 + ...  - h*f(t_np1,y)
           self.statistics["nfcns"] += 1
           res = alpha[0]*y - h*f(t_np1, y)
           for j in range(1, len(alpha)):
                 res += alpha[j] * Y[-j]
                
           return res
       

           #iterate until tolerance level is reached or maximum steps reached
           
        J = self._jacobian(G, y_np1_i)
        for i in range(self.maxit):
            #J = self._jacobian(G, y_np1_i)
            d = np.linalg.solve(J, G(y_np1_i))
            y_np1_ip1 = y_np1_i - d
            
            if np.linalg.norm(y_np1_ip1-y_np1_i) < self.tol:
                return t_np1,y_np1_ip1
            y_np1_i=y_np1_ip1
        else:
            raise Explicit_ODE_Exception('Corrector could not converge within % iterations'%i)
            

class BDF2(abstractBDF):
    tol=1.e-6    
    maxit=1000
    maxsteps=1000
    
    def is_at_init(self, i): 
        return i == 0

    
    def initial_steps(self, T, Y, h, i):
        t_n, y_n = T[-1], Y[-1]
        t_np1, y_np1 = self.step_EE(t_n, y_n, h)
        return t_np1, y_np1
        
    
    def step_BDF(self, T, Y, h):
        alpha=[3./2.,-2.,1./2]
        f=self.problem.rhs
        return self.newton(T,Y,alpha,h,f)
    
    def print_statistics(self, verbose=NORMAL):

        super().print_statistics(verbose=NORMAL)
        self.log_message(' Solver            : BDF2\n', verbose)

        
        

class BDF3(abstractBDF):
    tol=1.e-6    
    maxit=1000
    maxsteps=1000
    
    def is_at_init(self, i): 
        return i == 0 or i == 1

    
    def initial_steps(self, T, Y, h, i):
        t_n, y_n = T[-1], Y[-1]
        t_np1, y_np1 = self.step_EE(t_n, y_n, h)
        return t_np1, y_np1
        
    
    def step_BDF(self, T, Y, h):
        alpha = [11./6, -3, 3./2, -1./3]
        f=self.problem.rhs
        return self.newton(T,Y,alpha,h,f)

    
    def print_statistics(self, verbose=NORMAL):
        super().print_statistics(verbose=NORMAL)
        self.log_message(' Solver            : BDF3\n', verbose)
        
class BDF4(abstractBDF):
    tol=1.e-6    
    maxit=1000
    maxsteps=1000
    
    def is_at_init(self, i): 
        return i == 0 or i == 1 or i == 2

    
    def initial_steps(self, T, Y, h, i):
        t_n, y_n = T[-1], Y[-1]
        t_np1, y_np1 = self.step_EE(t_n, y_n, h)
        return t_np1, y_np1
        
    
    def step_BDF(self, T, Y, h):
        alpha = [25./12,-4,3,-4./3,1/4]
        f=self.problem.rhs
        return self.newton(T,Y,alpha,h,f)

    
    def print_statistics(self, verbose=NORMAL):
        super().print_statistics(verbose=NORMAL)
        self.log_message(' Solver            : BDF4\n', verbose)
            




        
        
        
