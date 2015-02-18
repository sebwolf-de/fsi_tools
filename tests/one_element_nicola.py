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
    
    mesh = Mesh('../gmsh_apps/one_element.msh') 
    mesh.read_data_from_input()
    
    quad = Quadrature()

    mapping = Mapping()
    shape = ReferenceShapeFunction(2)
    
    mapping.reinit(s_eval,t_eval)
    
    for row in mesh.topo:
        x_map = mesh.x[row]
        y_map = mesh.y[row]
        
        (xq,yq) = mapping.evaluate(x_map,y_map)
        
        DF = mapping.get_gradient(x_map,y_map) # DF[nq,2,2]
        
        for i in np.arange(0,n_local_basis):
            grad_u_i = shape.get_gradient(i)# grad_u[nq,1,2]
            
            for j in np.arange(0,n_local_basis):
                grad_u_j = shape.get_gradient(j)# grad_u[nq,1,2]
                
                for df, gu_i, gu_j in zip(DF,grad_u_i,grad_u_j):
                    grad_corrente_i = np.dot(df,gu_i)
                    grad_corrente_j = np.dot(df,gu_j)
                    
                    local_matrix(i,j) += scalar_product(grad_corrente_i,grad_corrente_j) * w_det#non so come sia in numpy
        #
        # sono su un elemento
    
        print row
        
   
#    print mapping.edge_topo
#    
#    print mesh.x
    A = np.zeros((len(mesh.nodes),len(mesh.nodes)))
#    print B
#    for mesh.row in mesh.topo:
#
#        mapping.iso_map ()  
#        xq = mapping.xp
#        yq = mapping.yp
#        det = mapping.det
#       
#        print '----------------'
#        print xq
#        print yq
#        print '----------------'
#        
#        shape = ReferenceShapeFunction(2)
#        value_list = shape.eval_basis_functions(xq,yq)
#        grad_list = shape.eval_basis_gradients(xq,yq)
#        
#        K = np.zeros((6,6)) 
#        for i in [0,1,2,3,4,5]:
#            for j in [0,1,2,3,4,5]:
#                B = grad_list[j]*grad_list[i]
#                C = np.sum(B, axis=1)
##                print C[:,0]
#                A[i,j] = np.sum(C[:,0]*quad.weights*det)
#        
#        h = 0  
#        for i in mesh.row:
#            k = 0            
#            for j in mesh.row:
#                B = grad_list[k]*grad_list[h]
#                C = np.sum(B, axis=1)
#                K[h,k] = np.sum(C[:,0]*quad.weights*det)                    
#                A[i,j] = A[i,j] + K[h,k]
#                k += 1
#            h += 1
#        
    print A