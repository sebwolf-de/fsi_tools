"""The Laplacian, the finite elements playground. 

.. moduleauthor:: Nicola Cavallini

"""

import numpy as np
import scipy.sparse.linalg as sp_la
import matplotlib.pyplot as pl

import lin_tri_mesh as lin_t3
import basis_func as shp
import assemble
import la_utils
import viewers


if __name__== '__main__':
    
    (topo,x,y) = lin_t3.load_msh('../gmsh_apps/unstr_square.msh')
    
    A = assemble.gradu_gradv_p1(topo,x,y)
    
    ndofs=A.shape[0]    
    
    rhs = np.zeros((ndofs,1))
    
    x_l = x[topo[0]]
    y_l = y[topo[0]]
    eval_points = np.zeros((0,2))
    (phi_dx,phi_dy,phi,omega) = shp.tri_p1(x_l,y_l,eval_points)
    
    for row in topo:
        a=x[row]
        b=y[row]
        surf_e = 1./2. * abs( a[0]*b[2] - a[0]*b[1] + a[1]*b[0] - a[1]*b[2] + a[2]*b[1] - a[2]*b[0] )
        local_rhs = 1./3. * np.ones((3,1)) * surf_e
        rhs[row] = rhs[row] + local_rhs
    
    bc_id = np.where( y == 0 )
    A = la_utils.set_diag(A,bc_id)
    rhs[bc_id] = 0
    
    bc_id = np.where( y == 1 )
    A = la_utils.set_diag(A,bc_id)
    rhs[bc_id] = 0
    
    bc_id = np.where( x == 1)
    A = la_utils.set_diag(A,bc_id)
    rhs[bc_id] = 0
    
    bc_id = np.where( x == 0 )
    A = la_utils.set_diag(A,bc_id)
    pl.spy(A)
    pl.show
        
    rhs[bc_id] = 0            
    sol = sp_la.spsolve(A,rhs)
    
    viewers.plot_sol_p1(x,y,sol,topo)
    
    