# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 15:20:55 2015

@author: nicola
"""

import numpy as np

class Quadrature:
    def __init__(self):
        self.weights = np.array([1./3.,1./3.,1./3.])
        
        self.points = np.array([[2./3.,1./3.,1./3.],
                                [1./3.,2./3.,1./3.],
                                [1./3.,1./3.,2./3.]])
        return