# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 18:04:38 2015

@author: michele
"""

import numpy as np
import scipy.sparse.linalg as sp_la
import matplotlib.pyplot as pl

import lin_tri_mesh as lin_t3
import basis_func as basis
import assemble
import la_utils
import viewers
import errore2D 
import esatta

if __name__== '__main__':
    
    
#    a=np.array((1,2,3))
#    #print a
#    b=np.zeros((3,1))    
#    b[:,0]=a.transpose()
#    #print b    

    (topo,x,y) = lin_t3.load_msh('../gmsh_apps/unstr_square.msh')
    #print topo.shape
    #print x[topo[137,:]]
    aq = np.array([1./90.,1./90.,1./90.,16./225.,16./225.,16./225.,(49./120.)**2,(49./120.)**2,(49./120.)**2,81./320.]) # pesi quadratura
    xq = np.array([0.,1.,0.,0.5,0.5,0.,1./7.,5./7.,1./7.,1./3.]) # coord x nodi quadratura riferimento
    yq = np.array([0.,0.,1.,0.,0.5,0.5,5./7.,1./7.,1./7.,1./3.]) # coord y nodi quadratura riferimento

    A = assemble.gradu_gradv_p1(topo,x,y)
    
    ndofs=A.shape[0]    
    
    rhs = np.zeros((ndofs,1))
    
    for row in topo:
        a=x[row]
        b=y[row]
        surf_e = 1./2. * abs( a[0]*b[2] - a[0]*b[1] + a[1]*b[0] - a[1]*b[2] + a[2]*b[1] - a[2]*b[0] )
        local_rhs = 1./3. * np.ones((3,1)) * surf_e
        rhs[row] = rhs[row] + local_rhs    
            
    bc_id = np.where( y == 0)
    A = la_utils.set_diag(A,bc_id)
    rhs[bc_id] = 0
        
    bc_id = np.where( y == 1)
    A = la_utils.set_diag(A,bc_id)
    rhs[bc_id] = 0
        
    bc_id = np.where( x == 0)
    A = la_utils.set_diag(A,bc_id)
    rhs[bc_id] = 0

    bc_id = np.where( x == 1)
    A = la_utils.set_diag(A,bc_id)
    rhs[bc_id] = 0
        
    sol = sp_la.spsolve(A,rhs)

    for row in topo:
        A=x[row]
        B=y[row]
        uu=sol[row]
        b11=A[1]-A[0]
        b12=A[2]-A[0]
        b21=B[1]-B[0]
        b22=B[2]-B[0]
        det=b11*b22-b12*b21
        xp = A[1] + b11*xq +b12*yq
        yp = B[1] + b21*xq +b22*yq
        eval_p = np.array([xp,yp])
        #print eval_p    
        eval_p = np.reshape ( eval_p, (10,2))
     
        tmpesatta=esatta.sol_esatta(xp,yp)
        #print tmpesatta.shape
        tmpapprox=np.zeros((1,xp.shape[0]))        
        #print tmpapprox.shape
        tmpesatta = np.reshape(tmpesatta, (1,xp.shape[0]))
        #print tmpesatta.shape
        (phi_dx,phi_dy,phi,omega) = basis.tri_p1(A,B,eval_p) 
         
        tmpderx = esatta.dx_sol_esatta(xp,yp)
        print tmpderx.shape
        tmpderx = np.reshape ( tmpderx, (1,xp.shape[0]))
        print tmpderx.shape
        
        #print uu[0]*phi[:,0]
        #print phi.shape[0]   # (3,1)
        #print uu       # (3,1)
        #print tmpapprox  #(10,1)
        for i in np.array([0,1,2]):            
            tmpapprox = tmpapprox + uu[i]*phi[:,i]
            print tmpapprox.shape    

#               eval_points = np.zeros((0,2))
#        print eval_points.shape[0]
#        print a,b
        
#    r= np.sqrt(x**2+y**2)
#    #print r
#    bc_id = np.where( r > 0.999 )
#    #print bc_id
#    
#    esatta = (x[row]*(x[row]-1)*y[row]*(y[row]-1))
#    derx=1
#    dery=2
#    sol=3
#    pesi=errore2D.err_l2(topo,x,y,sol)
#    print pesi