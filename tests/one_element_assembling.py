# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 15:15:32 2015

@author: nicola
"""

import numpy as np

from mesh import Mesh
from mapping import Mapping
from quadrature import Quadrature


if __name__== '__main__':
    
    mesh = Mesh('../gmsh_apps/one_element.msh') 
    mesh.read_data_from_input()
    
    mapping = Mapping(mesh)
    
    quad = Quadrature()
    
    print mesh.topo
    print mapping.edge_topo
    
    print mesh.x
    
    for row in mesh.topo:
        
        # this section shpuld go inside mapping

        edge_nodes = row[0:3]
        edge_x = mesh.x[edge_nodes]
        edge_y = mesh.y[edge_nodes]
        
        ref_simplex = np.vstack([edge_x,edge_y,np.ones((3,))])
        
        current_simplex = np.dot(ref_simplex,quad.points)
        
        x_qp = current_simplex[0,:]
        y_qp = current_simplex[1,:]
        
        #
        
        print '----------------'
        print x_qp
        print y_qp
        print '----------------'
        mesh.next()