#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 15:08:44 2022

@author: antonfredriksson
"""
from assimulo.explicit_ode import Explicit_ODE
import numpy as np
from assimulo.ode import Explicit_ODE_Exception, ID_PY_OK
from scipy.sparse.linalg import spsolve

class SecondOrder(Explicit_ODE):
    """Base class for second order differential equation solvers"""
    
    def __init__(self, problem, delta_t):
        super().__init__(problem)
        
        #Solver options
        self.options["max_steps"] = 10000
        self.delta_t = delta_t
        
        
    def integrate(self, t, y0, tf, opts):
        u = y0[0:int(len(y0)/2)]
        udot = y0[int(len(y0)/2):]
        
        f, M, C, K = self.problem.f, self.problem.M, self.problem.C, self.problem._get_K(u)

        udotdot = spsolve(M, f(t) - C@udot - K@u)
        
        tres = []
        ures = []
        
        #lhs_matrix = self.get_lhs_matrix() #cannot do this if nonlinear K
        
        for i in range(self.options["max_steps"]):
            if t >= tf:
                break
            self.delta_t = min(self.delta_t, abs(tf-t))
            t = t + self.delta_t
            
            u_new = self.solve_u_n(t, u, udot, udotdot)
            udotdot_new = self.solve_udotdot_n(u, u_new, udot, udotdot, t)
            udot_new = self.solve_udot_n(u, u_new, udot, udotdot, udotdot_new)
            
            tres.append(t)
            ures.append(np.concatenate((u_new,udot_new)))
            
            u, udot, udotdot = u_new, udot_new, udotdot_new
        else:
            raise Explicit_ODE_Exception('Final time not reached within maximum number of steps')
        
        return ID_PY_OK, tres, ures
        
        
        
    def solve_u_n(self, t, u, udot, udotdot):
        raise NotImplementedError()
        
    def solve_udot_n(self, u, u_new, udot, udotdot, udotdot_new):
        raise NotImplementedError()
        
    def solve_udotdot_n(self, u, u_new, udot, udotdot, t):
        raise NotImplementedError()
    
    def get_fourth_term(self, u_old):
        raise NotImplementedError()
        
class ExplicitNewmark(SecondOrder):
    """Explicit Newmark method, i.e with beta = 0, gamma = 1/2 and no damping"""
    
    def __init__(self, problem, delta_t = 0.01):
        super().__init__(problem, delta_t)
  
    def solve_u_n(self, t, u, udot, udotdot):
        return u + udot*self.delta_t+udotdot*self.delta_t**2/2
    
    def solve_udot_n(self, u, u_new, udot, udotdot, udotdot_new):
        return udot + udotdot/2*self.delta_t + udotdot_new*self.delta_t/2
    
    def solve_udotdot_n(self, u, u_new, udot, udotdot, t):
        return np.linalg.solve(self.problem.M, self.problem.f(t) - self.problem._get_K(u_new)@u_new)
        
class ImplicitSecondOrder(SecondOrder):
    """Base class for the Implicit second order methods, e.g. HHT and Implicit Newmark"""
    
    def __init__(self, problem, beta, gamma, delta_t):
        super().__init__(problem, delta_t)
        self.beta = beta
        self.gamma = gamma
        
    def solve_u_n(self, tnew, u, udot, udotdot):
        f, M, C = self.problem.f, self.problem.M, self.problem.C
        beta, gamma, delta_t = self.beta, self.gamma, self.delta_t
        lhs_matrix = self.get_lhs_matrix(u)
        
        second_term = M @ (u/(beta*delta_t**2) + udot/(beta*delta_t) + \
                           (1/(2*beta) - 1)*udotdot)
        third_term = C @ (gamma*u/(beta*delta_t) - \
                          (1 - gamma/beta)*udot - \
                          (1 - gamma/(2*beta))*delta_t*udotdot)
        fourth_term = self.get_fourth_term(u)
            
        rhs = f(tnew) + second_term + third_term + fourth_term
            
        return spsolve(lhs_matrix, rhs)
        
    def solve_udot_n(self, u, u_new, udot, udotdot, udotdot_new):
        return self.gamma/self.beta*(u_new-u)/self.delta_t + (1-self.gamma/self.beta)*udot +  \
        (1-self.gamma/(2*self.beta))*self.delta_t*udotdot
        
    def solve_udotdot_n(self, u, u_new, udot, udotdot, t):
        return (u_new-u)/(self.beta*self.delta_t**2) - udot/(self.beta*self.delta_t) -\
        (1/(2*self.beta)-1)*udotdot
        
class ImplicitNewmark(ImplicitSecondOrder):
    """Implicit Newmark method"""
    
    def __init__(self, problem, beta, gamma, delta_t = 0.01):
        super().__init__(problem, beta, gamma, delta_t)
        

    def get_fourth_term(self, u_old):
        size = len(u_old)
        return np.zeros(size)
    
    def get_lhs_matrix(self, u):
        M, C, K = self.problem.M, self.problem.C, self.problem._get_K(u)
        beta, gamma, delta_t = self.beta, self.gamma, self.delta_t        
        return M/(beta*delta_t**2) + gamma*C/(beta*delta_t) + K
    
    
        
class HHT(ImplicitSecondOrder):
    """HHT alpha method"""
    
    def __init__(self, problem, beta, gamma, alpha, delta_t = 0.01):
        super().__init__(problem, beta, gamma, delta_t)
        self.alpha = alpha
        
    def get_fourth_term(self, u_old):
        return self.alpha*self.problem._get_K(u_old)@u_old
    
    def get_lhs_matrix(self, u):
        M, C, K = self.problem.M, self.problem.C, self.problem._get_K(u)
        beta, gamma, delta_t, alpha = self.beta, self.gamma, self.delta_t, self.alpha           
        return M/(beta*delta_t**2) + gamma*C/(beta*delta_t) + (1+alpha)*K
        

        
        