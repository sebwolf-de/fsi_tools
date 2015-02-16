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
    
    cerchio = Mesh('../gmsh_apps/cerchio4_1.msh') 
    cerchio.read_data_from_input()
    topo = cerchio.topo

    quad = Quadrature()
    mapping = Mapping(cerchio,quad) 
    iso = mapping.iso_map ()  
    xq = iso.xp
    yq = iso.yp
    det = iso.det
#    x = cerchio.x
#    y = cerchio.y
#    nodes = cerchio.nodes
#    b_nodes = cerchio.b_nodes
#    connettivity = Mapping (cerchio)
#    topo_1 = connettivity.edge_topo
    print topo
#    print topo_1

    print topo.shape[1]
    print quad.points

