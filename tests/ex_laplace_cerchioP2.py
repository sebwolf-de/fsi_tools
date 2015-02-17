# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 15:15:32 2015

@author: nicola
"""

import numpy as np
import scipy.linalg as sp_la
import matplotlib.pyplot as pl
import la_utils
import viewers
import esatta
import errore2D

from mesh import Mesh
from mapping import Mapping
from quadrature import Quadrature
from reference_shape_function import ReferenceShapeFunction 

if __name__== '__main__':
    
    mesh = Mesh('../gmsh_apps/cerchio4_1.msh') 
    mesh.read_data_from_input()
    print mesh.topo 
#    print mesh.b_nodes   
    quad = Quadrature()

    mapping = Mapping(mesh,quad)    
        
   
#    print mapping.edge_topo
#    
#    print mesh.x
    rhs = np.zeros((len(mesh.nodes),1))
    A = np.zeros((len(mesh.nodes),len(mesh.nodes)))
#    print A    
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
        xp = mapping.xp
        yp = mapping.yp
        det = mapping.det



#        xq = 1 - x_qp - y_qp        
#        yq = x_qp
#        
#        det = 1 + 0*xq  #(x1-x0)*(y2-y0)-(x2-x0)*(y1-y0)
#============================================================        
       
#        print '----------------'
#        print xp
#        print yp
#        print '----------------'
        
        shape = ReferenceShapeFunction(2)
        value_list = shape.eval_basis_functions(xp,yp)
        grad_list = shape.eval_basis_gradients(xp,yp)
        load = esatta.load_2(xp,yp)
#        print load
        K = np.zeros((len(mesh.row),len(mesh.row))) 
#        for i in [0,1,2,3,4,5]:
#            for j in [0,1,2,3,4,5]:
#                B = grad_list[j]*grad_list[i]
#                C = np.sum(B, axis=1)
##                print C[:,0]
#                A[i,j] = np.sum(C[:,0]*quad.weights*det)
        
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
            rhs[i] = rhs[i]+np.sum(value_list[h]*load*(1./6)*np.abs(det))  
            h += 1
#    print rhs     
#    print A
#    pl.spy(A)
#    pl.show
    for i in mesh.b_nodes:
        for j in mesh.nodes:
            if i!=j:            
                A[i,j] = 0 # la_utils.set_diag(A,mesh.b_nodes)
    rhs[mesh.b_nodes] = 0
    print A 
    print rhs
    pl.spy(A)
    pl.show

    sol = sp_la.solve(A,rhs)
#    print sol.shape
    sol1 = [] 
    for i in sol:
        sol1 = np.append(sol1,i)
#    print sol1.shape
#    print mesh.x.shape
#    print mesh.y.shape
    fig = pl.figure(2)
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(mesh.x, mesh.y, sol1)
    pl.show()
    
    fig = pl.figure(3)
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(mesh.x, mesh.y, esatta.sol_esatta_2(mesh.x,mesh.y))
    pl.show()

#    viewers.plot_sol_p1(mesh.x,mesh.y,sol,mesh.topo)

#    ( err_l2 , er_derx , er_dery ) = errore2D.err_l2(mesh.topo,x,y,sol)