#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 13:51:06 2022

@author: antonfredriksson
"""
from assimulo.problem import Explicit_Problem
import numpy as np
from scipy.sparse.linalg import spsolve

class Explicit_Problem_2nd(Explicit_Problem):
    """Problem class for 2nd order differential equations"""
    
    def __init__(self, f, M, C, K, u0, udot0):
        super().__init__(self)
        self.f = f
        self.M = M
        self.C = C
        self.K = K
        self.y0 = np.concatenate((u0, udot0))
        
    def _get_K(self, u):
        if callable(self.K):
            return self.K(u)
        else:
            return self.K
    
    def rhs(self, t, y):
        f, M, C, K = self.f, self.M, self.C, self.K(y[0:int(len(y)/2)])
        M_inv_K = spsolve(M, K)
        M_inv_C = spsolve(M, C)
        shape = M.shape
        block_matrix = np.block([[np.zeros(shape), -np.eye(shape[0])], [M_inv_K, M_inv_C]])
        block_vector = np.concatenate((np.zeros(shape[0]), f(t)))
        
        return block_vector - block_matrix @ y
        
    