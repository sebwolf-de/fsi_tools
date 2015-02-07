# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 14:59:14 2015

@author: nicola
"""
import numpy as np
from mesh import Mesh

vel_mesh = Mesh('vel_file')
vel_mesh.read_data() # now this is doing noting mik: implementation required
# you only need to locate the correct function and put int in the right place

prex_mesh = Mesh('prex_file')
prex_mesh.read_data() # now this is doing noting mik: implementation required

# here I artificially add therse two connettivity matrices
vel_mesh.topo = np.array([[0,1,2,3,4,5,6],[7,8,9,9,4,10,11]])
prex_mesh.topo = np.array([[0,1,2],[3,4,5,6]])

for v_row,p_row in zip(vel_mesh.topo,prex_mesh.topo):
    print v_row
    print p_row