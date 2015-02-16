# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 15:14:45 2015

@author: nicola
"""
import numpy as np
from mesh import Mesh
#from reference_shape_function import ReferenceShapeFunction
from quadrature import Quadrature
#>>> a = np.array([[1,2,3,7],[4,5,6,8],[9,10,11,12]])
#>>> print a
#[[ 1  2  3  7]
# [ 4  5  6  8]
# [ 9 10 11 12]]
#>>> b = a[:,1]
#>>> print b
#[ 2  5 10]
#>>> c = a[:,[0,1]]
#>>> print c
#[[ 1  2]
# [ 4  5]
# [ 9 10]]
#>>> 
class Mapping:
    def __init__(self, Mesh , Quadrature):
        
        self.topo = Mesh.topo        
        self.edge_topo = self.topo [:,0:3]        # connettivity matrix corresponding only to the 
                                                 # edge nodes of the mesh 
                                                 # (deg2: exclude mid segment nodes)
        self.x = Mesh.x                                    # edge nodes
        self.y = Mesh.y   
        self.row = Mesh.row
#        self.yp = np.array([])
#        self.det = 0                               
        return

    def iso_map (self):
        
        self.xp = np.array([])
        self.yp = np.array([])
        self.det = 0
               
        Qpt_bar = np.array([[2./3.,1./6.,1./6.],
                            [1./6.,2./3.,1./6.],
                            [1./6.,1./6.,2./3.]])                         #
        
        Qpoints = np.dot(np.array([[0,1,0],[0,0,1],[1,1,1]]),Qpt_bar)             
#        print Qpt_bar, Qpoints        
        xq = Qpoints[0,:] 
        yq = Qpoints[1,:]
#        print xq ,yq
#        self.xp = np.array([])
#        self.yp = np.array([])
#        self.det = 0
#        topo = Mesh.topo 
        
        if self.row in self.topo:
        
            a1=self.x[self.row][1]-self.x[self.row][0]
            a2=self.x[self.row][2]-self.x[self.row][0]

            b1=self.y[self.row][1]-self.y[self.row][0]
            b2=self.y[self.row][2]-self.y[self.row][0]
            
            if self.topo.shape[1]==6:
           
                a3=self.x[self.row][3]-self.x[self.row][0]
                a4=self.x[self.row][4]-self.x[self.row][0]
                a5=self.x[self.row][5]-self.x[self.row][0]
            
                b3=self.y[self.row][3]-self.y[self.row][0]
                b4=self.y[self.row][4]-self.y[self.row][0]
                b5=self.y[self.row][5]-self.y[self.row][0]
                
                A1=2*a1-4*a3
                A2=2*a2-4*a5
                A3=4*(a4-a3-a5)
                A4=4*a3-a1
                A5=4*a5-a2
                
                B1=2*b1-4*b3
                B2=2*b2-4*b5
                B3=4*(b4-b3-b5)
                B4=4*b3-b1
                B5=4*b5-b2
            
                        
                self.xp = self.x[self.row][0] + (A1*(xq**2) + A2*(yq**2) + A3*xq*yq + A4*xq + A5*yq)     
                self.yp = self.y[self.row][0] + (B1*(xq**2) + B2*(yq**2) + B3*xq*yq + B4*xq + B5*yq) 
                self.det = (A1*xq*2+A3*yq+A4)*(B2*yq*2+B3*xq+B5)-(A2*yq*2+A3*xq+A5)*(B1*xq*2+B3*yq+B4)
            else:
                self.xp = self.x[self.row][0] + a1*xq+a2*yq
                self.yp = self.y[self.row][0] + b1*xq+b2*yq
                self.det = a1*b2-a2*b1 
            
        return self.xp,self.yp,self.det