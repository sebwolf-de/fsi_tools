# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 15:14:45 2015

@author: nicola
"""

class Mapping:
    def __init__(self,mesh):
        self.edge_topo # connettivity matrix corresponding only to the 
                       # edge nodes of the mesh 
                       # (deg2: exclude mid segment nodes)
        self.x # edge nodes
        self.y # edge nodes
        return