# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 15:14:45 2015

@author: nicola
"""

from mesh import Mesh
#>>> a = np.array([[1,2,3,7],[4,5,6,8],[9,10,11,12]])
#>>> print a
#[[ 1  2  3  7]
# [ 4  5  6  8]
# [ 9 10 11 12]]
#>>> b = a[:,1]
#>>> print b
#[ 2  5 10]
#>>> c = a[:,[0,1]]
#>>> print c
#[[ 1  2]
# [ 4  5]
# [ 9 10]]
#>>> 
class Mapping:
    def __init__(self,Mesh):
                
        self.edge_topo = Mesh.topo[:,0:3]        # connettivity matrix corresponding only to the 
                                                 # edge nodes of the mesh 
                                                 # (deg2: exclude mid segment nodes)
        self.x = Mesh.x                                    # edge nodes
        self.y = Mesh.y                                  # edge nodes
        return