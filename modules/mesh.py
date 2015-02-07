# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 15:03:47 2015

@author: nicola
"""
import numpy as np

class Mesh:
    def __init__(self,file_to_be_read):
        self.filename = file_to_be_read
        return
        
    def read_data_from_input(self):
        self.topo = np.array([[]])
        self.x = np.array([])
        self.y = np.array([])
        return
        