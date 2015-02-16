# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 16:28:37 2015

@author: michele
"""

for row in topo:
    a1=x[row][1]-x[row][0]
    a2=x[row][2]-x[row][0]
    a3=x[row][3]-x[row][0]
    a4=x[row][4]-x[row][0]
    a5=x[row][5]-x[row][0]

    b1=y[row][1]-y[row][0]
    b2=y[row][2]-y[row][0]
    b3=y[row][3]-y[row][0]
    b4=y[row][4]-y[row][0]
    b5=y[row][5]-y[row][0]
    
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
    
    x = x[row][0] + A1*xq^2 + A2*yq^2 + A3*xq*yq + A4*xq + A5*yq    
    y = y[row][0] + B1*xq^2 + B2*yq^2 + B3*xq*yq + B4*xq + B5*yq
    
    det = (A1*xq*2+A3*yq+A4)*(B2*yq*2+B3*xq+B5)-(A2*yq*2+A3*xq+A5)*(B1*xq*2+B3*yq+B4)
    
    
    

