"""
Created on Thu Jan 15 10:57:39 2015

@author: michele
"""
import numpy as np
import basis_func as basis
import esatta

def err_l2(topo,x,y,sol):
    aq = np.array([1./90.,1./90.,1./90.,16./225.,16./225.,16./225.,(49./120.)**2,(49./120.)**2,(49./120.)**2,81./320.]) # pesi quadratura
    xq = np.array([0.,1.,0.,0.5,0.5,0.,1./7.,5./7.,1./7.,1./3.]) # coord x nodi quadratura riferimento
    yq = np.array([0.,0.,1.,0.,0.5,0.5,5./7.,1./7.,1./7.,1./3.]) # coord y nodi quadratura riferimento
#    eval_p = np.array([xq,yq])
#    eval_p = np.reshape ( eval_p, (xq.shape[0],2))    
    er_l2 = 0.
    ex_l2 = 0.
    er_derx = 0.
    ex_derx = 0.
    er_dery = 0.
    ex_dery = 0.
    for row in topo:
        A=x[row]
        B=y[row]
        uu=sol[row]
        b11=A[1]-A[0]
        b12=A[2]-A[0]
        b21=B[1]-B[0]
        b22=B[2]-B[0]
        #det=b11*b22-b12*b21
        xp = A[1] + b11*xq +b12*yq
        yp = B[1] + b21*xq +b22*yq
        eval_p = np.array([xp,yp])
        eval_p = np.reshape ( eval_p, (xp.shape[0],2))
       #eval_p = eval_p.transpose
#       print eval_p.shape[0]
        
        tmpesatta = esatta.sol_esatta(xp,yp)
        tmpesatta = np.reshape ( tmpesatta, (1,xp.shape[0]))        
        tmpderx = esatta.dx_sol_esatta(xp,yp)
        tmpderx = np.reshape ( tmpderx, (1,xp.shape[0]))
        tmpdery = esatta.dy_sol_esatta(xp,yp)
        tmpdery = np.reshape (tmpdery , (1,xp.shape[0]))        
        tmpapprox = np.zeros((1,xp.shape[0]))
        tmphderx = np.zeros((1,xp.shape[0]))
        tmphdery = np.zeros((1,xp.shape[0]))
        
        (phi_dx,phi_dy,phi,omega) = basis.tri_p1(A,B,eval_p) 
        
        for i in np.array([0,1,2]):            
            tmpapprox = tmpapprox + uu[i]*phi[:,i]
            tmphderx = tmphderx + uu[i] * (phi_dx[:,i]*b22 - phi_dy[:,i]*b21)
            tmphdery = tmphdery + uu[i] * (phi_dy[:,i]*b11 - phi_dx[:,i]*b12)            
        
        tmp = np.sum ( aq * ( tmpesatta - tmpapprox )**2 )
        tmp2 = np.sum (aq * ( tmpesatta )**2)        
        tmpx=np.sum (aq * ( tmpderx - tmphderx )**2)
        tmp2x=np.sum ( aq * tmpderx**2 )
        tmpy=np.sum ( aq * ( tmpdery - tmphdery )**2 )
        tmp2y=np.sum ( aq * tmpdery**2 )
        er_derx = er_derx + tmpx
        ex_derx = ex_derx + tmp2x
        er_dery = er_dery + tmpy
        ex_dery = ex_dery + tmp2y
        
        er_l2 = er_l2 + tmp
        ex_l2 = ex_l2 + tmp2

    er_l2 = np.sqrt ( er_l2/ex_l2 )
    er_derx = np.sqrt(er_derx/ex_derx)
    er_dery = np.sqrt(er_dery/ex_dery)
    
    return er_l2 , er_derx , er_dery