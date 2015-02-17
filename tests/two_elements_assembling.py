# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 15:15:32 2015

@author: nicola
"""

import numpy as np

from mesh import Mesh
from mapping import Mapping
from quadrature import Quadrature
from reference_shape_function import ReferenceShapeFunction 

if __name__== '__main__':
    
    mesh = Mesh('../gmsh_apps/two_elements.msh') 
    mesh.read_data_from_input()
#    print mesh.topo 
    
    quad = Quadrature()

    mapping = Mapping(mesh,quad)    
        
   
#    print mapping.edge_topo
#    
#    print mesh.x
    A = np.zeros((len(mesh.nodes),len(mesh.nodes)))
#    print B    
    for mesh.row in mesh.topo:
#        print mesh.row   
        # this section shpuld go inside mapping

#        edge_nodes = row[0:3]
#        edge_x = mesh.x[edge_nodes]
#        edge_y = mesh.y[edge_nodes]
#        
#        
#        ref_simplex = np.vstack([edge_x,edge_y,np.ones((3,))])
#        
#        current_simplex = np.dot(ref_simplex,quad.points)
#        
#        x_qp = current_simplex[0,:]
#        y_qp = current_simplex[1,:]

#===========================================================
#     Trasformazione da coordinate globali a locali  
        mapping.iso_map ()  
        xq = mapping.xp
        yq = mapping.yp
        det = mapping.det



#        xq = 1 - x_qp - y_qp        
#        yq = x_qp
#        
#        det = 1 + 0*xq  #(x1-x0)*(y2-y0)-(x2-x0)*(y1-y0)
#============================================================        
       
        print '----------------'
        print xq
        print yq
        print '----------------'
        
        shape = ReferenceShapeFunction(2)
        value_list = shape.eval_basis_functions(xq,yq)
        grad_list = shape.eval_basis_gradients(xq,yq)
        
        K = np.zeros((6,6)) 
        for i in [0,1,2,3,4,5]:
            for j in [0,1,2,3,4,5]:
                B = grad_list[j]*grad_list[i]
                C = np.sum(B, axis=1)
#                print C[:,0]
                A[i,j] = np.sum(C[:,0]*quad.weights*det)
        
#        print 'matrice di elemento \n' 
#        print A
#
#        B = np.zeros(len(topo)+1,len(topo)+1)
        h = 0  
        for i in mesh.row:
            k = 0            
            for j in mesh.row:
                B = grad_list[k]*grad_list[h]
                C = np.sum(B, axis=1)
                K[h,k] = np.sum(C[:,0]*quad.weights*det)                    
                A[i,j] = A[i,j] + K[h,k]
                k += 1
            h += 1
#        
    print A