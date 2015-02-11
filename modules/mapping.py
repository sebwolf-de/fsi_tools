# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 15:14:45 2015

@author: nicola
"""

from mesh import Mesh

class Mapping:
    def __init__(self,Mesh):
                
        self.edge_topo = Mesh.topo[:,0:3]        # connettivity matrix corresponding only to the 
                                                 # edge nodes of the mesh 
                                                 # (deg2: exclude mid segment nodes)
        self.x = Mesh.x                                    # edge nodes
        self.y = Mesh.y                                  # edge nodes
        return