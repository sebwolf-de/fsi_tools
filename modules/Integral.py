# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 14:25:36 2015

@author: michele
"""
import esatta
import shape_2D
def integral(A,B,f,eval_points):
    """
    INPUT:
    A = array with x-coordinates of triangle vertices
    B = array with y-coordinates of triangle vertices    
    f = function to integrate
    If function "f" is defined only on the reference triangle then A = [0,1,0] & A = [0,0,1]
    eval_points= number of elauations points
    
    OUTPUT:
    Integral of f on the triangle (A[0],B[0]), (A[1],B[1]), (A[2],B[2])
     """   
    
    if eval_points ==10:
        aq = np.array([1./90.,1./90.,1./90.,16./225.,16./225.,16./225.,(49./120.)**2,(49./120.)**2,(49./120.)**2,81./320.]) # pesi quadratura
        xq = np.array([0.,1.,0.,0.5,0.5,0.,1./7.,5./7.,1./7.,1./3.]) # coord x nodi quadratura riferimento
        yq = np.array([0.,0.,1.,0.,0.5,0.5,5./7.,1./7.,1./7.,1./3.]) # coord y nodi quadratura riferimento    
       
        b11=A[1]-A[0]
        b12=A[2]-A[0]
        b21=B[1]-B[0]        
        b22=B[2]-B[0]
        det=b11*b22-b12*b21        
        xp = A[0] + b11*xq +b12*yq
        yp = B[0] + b21*xq +b22*yq
        ff=f(xp,yp)
        Integral = np.sum (aq * ff * np.abs(det)) 
    return Integral
        