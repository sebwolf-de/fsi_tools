# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 15:15:13 2015

@author: michele
"""

import numpy as np
import scipy.sparse.linalg as sp_la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 

import lin_tri_mesh as lin_t3
import basis_func as shp
import assemble
import la_utils
import viewers
import esatta
import errore2D
import basis_func as basis

if __name__== '__main__':
   mesh=[ "../gmsh_apps/settore_circolare.msh", "../gmsh_apps/settore_circolare1.msh", "../gmsh_apps/settore_circolare2.msh"]
   #mesh=[ "../gmsh_apps/settore_circolare.msh"]
   for m in mesh :
        (topo,x,y,nodes, b_nodes,int_nodes) = lin_t3.load_msh_1(m) #  numeration that starts from zero in all this elements 
        A= assemble.gradu_gradv_p1(topo,x,y)    # this matrix consider all nodes (b_nodes + int_nodes) and numeration starts from zero
                                                
        ndofs=A.shape[0]
        esattta = esatta.sol_esatta1 (x,y)
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_trisurf(x, y, esattta)
        plt.show()
    
        rhs = np.zeros((ndofs,1))       
        
        for row in topo:
            local_bc = np.zeros((3,1))
            a=x[row]
            b=y[row]
            surf_e = 1./2. * abs( a[0]*b[2]-a[0]*b[1]+a[1]*b[0]-a[1]*b[2]+a[2]*b[1]-a[2]*b[0] )
            tmpload = esatta.load1(a,b)
            tmpload = np.reshape ( tmpload, rhs[row].shape)
            local_rhs = 1./3. * tmpload * surf_e                         # three nodes integration formula
            rhs[row] = rhs[row] + local_rhs
       
        bc  = np.zeros((ndofs,1)) 
        
        bc_tmp = esatta.sol_esatta1 (x[b_nodes],y[b_nodes])
        k=0    
        for i in b_nodes:
            bc[i] = bc_tmp[k]
            k=k+1
        
        bc_rhs= A.dot(bc)
        rhs = rhs - bc_rhs
        
#        fig = plt.figure()
#        ax = fig.gca(projection='3d')
#        ax.plot_trisurf(x, y, bc)
#        plt.show()
                   
        A = la_utils.set_diag(A,b_nodes)
        rhs[b_nodes] = 0
        
        sol = sp_la.spsolve(A,rhs)
        sol = np.reshape(sol,bc.shape)
        sol = sol + bc
        solu=[]    
        for i in sol :
            solu= np.append(solu,i)
#            
#        fig = plt.figure()
#        ax = fig.gca(projection='3d')
#        ax.plot_trisurf(x, y, solu)
#        plt.show()
#        r = np.sqrt(x**2+y**2)
#        x_0=[]        
#        y_0=[]
#        solu_0=[]        
#        k=0        
#        for rho in r:
#            if rho > 1e-5:
#                x_0 = np.append(x_0,x[k])
#                y_0 = np.append(y_0,y[k])
#                solu_0 = np.append(solu_0,solu[k])
#                
#            k=k+1
                    
        ( err_l2 , er_derx , er_dery ) = errore2D.err_l2(topo,x,y,solu)
        print 'errore_L2' ,  err_l2 
        print 'errore_dx' ,  er_derx
        print 'errore_dy' ,  er_dery
        

       
