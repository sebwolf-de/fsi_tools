"""
Created on Wed Jan 21 12:34:37 2015

@author: michele
"""

import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 

import lin_tri_mesh as lin_t3
import basis_func as shp
import assemble
import la_utils
import viewers
import esatta
#import nodes_split1
from mesh import Mesh
from mapping import Mapping
from quadrature import Quadrature
from reference_shape_function import ReferenceShapeFunction 

if __name__== '__main__':
    
    mesh = Mesh('../gmsh_apps/cerchio4_1.msh') 
    mesh.read_data_from_input()
    print mesh.topo
    print len(mesh.row)
    print mesh.nodes
       

    quad = Quadrature()
    mapping = Mapping(mesh,quad) 
    mapping.iso_map ()  
    xp = mapping.xp
    yp = mapping.yp
    det = mapping.det
    shape = ReferenceShapeFunction(2)
    value_list = shape.eval_basis_functions(xp,yp)
    load = esatta.load(xp,yp)
    f = np.zeros((len(mesh.nodes),1))
    f[1]= f[1]+np.sum(value_list[4]*load*quad.weights*det) 
    print f
##    x = cerchio.x
##    y = cerchio.y
##    nodes = cerchio.nodes
##    b_nodes = cerchio.b_nodes
##    connettivity = Mapping (cerchio)
##    topo_1 = connettivity.edge_topo
#    print topo
##    print topo_1
#
#    print topo.shape[1]
#    print quad.points
#    B = np.zeros(len(topo)+1,len(topo)+1)
#    for row in topo:
#        B[row,row] = B[row,row]+ A  
#    for i in np.arange(len(topo)+1):
#        for j in np.arange(len(topo)+1):
#            B = B 
#            
    
