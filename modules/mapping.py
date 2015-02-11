# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 15:14:45 2015

@author: nicola
"""

from mesh import Mesh

class Mapping:
    def __init__(self,mesh):
        self.edge_topo =  mesh.topo[:0:3]        # connettivity matrix corresponding only to the 
                                                 # edge nodes of the mesh 
                                                 # (deg2: exclude mid segment nodes)
        self.x = mesh.x                                    # edge nodes
        self.y = mesh.y                                  # edge nodes
        return