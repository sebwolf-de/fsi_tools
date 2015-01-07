import numpy as np
import scipy.sparse.linalg as sp_la
import matplotlib.pyplot as pl

import lin_tri_mesh as lin_t3
import basis_func as shp
import assemble
import la_utils
import viewers

if __name__== '__main__':
   
    (topo,x,y) = lin_t3.load_msh('../gmsh_apps/cerchio_unstr.msh')
    
    A= assemble.gradu_gradv_p1(topo,x,y)
    
    ndofs=A.shape[0]
    print (ndofs)    
    rhs = np.zeros((ndofs,1))


    x_l = x[topo[0]]
    y_l = y[topo[0]]
    print (x_l , y_l)
    eval_points = np.zeros((0,2))
    (phi_dx,phi_dy,phi,omega) = shp.tri_p1(x_l,y_l,eval_points)
    
    for row in topo:
        local_rhs = 1./3. * np.ones((3,1)) * omega
        rhs[row] = rhs[row] + local_rhs
    
    #print rhs	    
    r= np.sqrt(x**2+y**2)
    bc_id = np.where( r > 0.299)
    A = la_utils.set_diag(A,bc_id)
    pl.spy(A)
    pl.show()    
    rhs[bc_id] = 0
    sol = sp_la.spsolve(A,rhs)
    
    viewers.plot_sol_p1(x,y,sol,topo)
