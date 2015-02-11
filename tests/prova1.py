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

if __name__== '__main__':
    
    cerchio = Mesh('../gmsh_apps/cerchio4_1.msh') 
    cerchio.read_data_from_input()
    topo = cerchio.topo
    x = cerchio.x
    y = cerchio.y
    nodes = cerchio.nodes
    b_nodes = cerchio.b_nodes
#    (topo, x, y, b_nodes, b_nodes_1, b_nodes_2, int_nodes) = nodes_split1.load_msh('../gmsh_apps/cerchio4_1.msh')
#    
    connettivity = Mapping (cerchio)
    topo_1 = connettivity.edge_topo
    print topo
#    print x
#    print y
#    print nodes
    print topo_1
#    print b_nodes_2
#    print int_nodes

#    print "esa_x",esa_x.shape
#    esa_y = esatta.der_sol_esatta_1(x,y)[1]
#    print "esa_y", esa_y.shape
#    (topo,x,y) = lin_t3.load_msh('../gmsh_apps/settore_circolare.msh')
#    A= assemble.gradu_gradv_p1(topo,x,y)
##    ndofs=A.shape[0]        
#    r = np.sqrt(x**2+y**2)    
#    bc_id0 = np.where( r > 0.999 )
##    print bc_id
#    bc_id1 = np.where ((r<=0.999) & (abs(x) < 1e-3 ) & (y <= 0 ) )
#    bc_id2 = np.where ( (r<=0.999) & (abs(y) < 1e-3 ) & ( x > 0 ) ) 
#      
#    bc_id3 = np.append (bc_id0 ,bc_id1)
#    bc_id = np.append (bc_id3 , bc_id2)
#    
##    bc_id = np.append(bc_id3,bc_id2)  
#    #print bc_id 
##    bc_id = np.array(bc_id3)
#    print bc_id    
##    i=1    
##    for k in bc_id:
##        for j in bc_id[i:]:         
##            if k==j: 
##               bc_id.remove(k)
##        i=i+1
##    print bc_id
#    
#    a= np.array([1,0,2,3])
#    A = sparse.crs_matrix((a.shape[0],a.shape[0]))
#    for i in a:
#        for j in a:
#            A[i,j]= (j-1)*(i+1) 
#            
#    print A