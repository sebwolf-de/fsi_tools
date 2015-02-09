# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 12:48:43 2015
@author: michele
"""

class ReferenceShapeFuncion:
    def __init__(self,in_degree):
        self.degree = in_degree
        return
#==============================================================================
    def eval_basis_functions(self,quad):

        x = quad[:0]
        y = quad[:1]

        if self.degree==1:
            value_list = [1.-x-y , x , y]
        if self.degree==2:      
            value_list = [2*(x**2+y**2)+4*x*y-3*(x+y)+1 , 2*(x**2)-x , 2*(y**2)-y , -4*(x**2)-4*x*y+4*x , 4*x*y , -4*(y**2)-4*x*y+4*y ]
        return value_list
#==============================================================================
    def eval_basis_gradients(self,quad):

        x = quad[:0]
        y = quad[:1]

        if self.degree==1:        
            grad_list = [[-1,-1],[1,0],[0,1]]      
        if self.degree==2:
            grad_list =[[4*x+4*y-3 , 4*y+4*x-3],[4*x-1 , 0],[0 , 4*y-1],[-8*x-4*y+4 , 4*x],[4*y , 4*x],[4*y , -8*y-4*x+4]]
             
        # for every quadrature point we create
        # a 2d numpy array. rows are the contravariant 
        # component, colums are the covariant.
        # grad \bf u = \partial_j u_i
        #
        # resituisci i gradienti
        return grad_list
#