# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 11:10:04 2015

@author: michele
"""
import numpy as np
def shape_2D(node,degree,der,x,y):
    """    
      INPUT
      node = 0,1,2,(3,4,5)    # shape function relative to "node"
      degree = 1,2  # degree of polynomials   
      der=0,1,2       
      x,y coordinates
    
      OUTPUT
      spape = value of the ith shape function
      der=0 --> function
      der=1 --> derivative in x
      der=2 --> derivative in y
    """
    if degree==1:
        print "degree value  ", degree
        if der==0:
            if node==0:
                val=1.-x-y
            if node==1:
                val=x
            if node==2:
                val=y
            else:
                print "node value  ", node ,"  not valid"
        if der==1:
            if node==0:
                val_x=-1
            if node==1:
                val_x=1
            if node==2:
                val_x=0
            else:
                print "node value  ", node ,"  not valid"
        if der==2:
            if node==0:
                val_x=-1
            if node==1:
                val_x=0
            if node==2:
                val_x=1
            else:
                print "node value  ", node ,"  not valid"
        else:
            print "derivative value  ", der ,"  not valid"
    
    if degree==2:
        print "degree value  ", degree
        if der==0:
            if node==0:
                val=2*(x**2+y**2)+4*x*y-3*(x+y)+1
            if node==1:
                val=2*(x**2)-x
            if node==2:
                val=2*(y**2)-y
            if node==3:
                val=-4*(x**2)-4*x*y+4*x 
            if node==4:
                val=4*x*y
            if node==5:
                val=-4*(y**2)-4*x*y+4*y                
            else:
                print "node value  ", node ,"  not valid"
        if der==1:
            if node==0:
                val_x=4*x+4*y-3
            if node==1:
                val_x=4*x-1 
            if node==2:
                val_x=0
            if node==3:
                val_x=-8*x-4*y+4 
            if node==4:
                val_x=4*y
            if node==5:
                val_x=4*y
            else:
                print "node value  ", node ,"  not valid"
        if der==2:
            if node==0:
                val_y=4*y+4*x-3
            if node==1:
                val_y=0
            if node==2:
                val_y=4*y-1
            if node==3:
                val_y=4*x 
            if node==4:
                val_y=4*x
            if node==5:
                val_y=-8*y-4*x+4 
            else:
                print "node value  ", node ,"  not valid"
    else:
            print "degree value  ", degree ,"  not valid"         

    return val, val_x, val_y
#
#        
#        error(['j=',j,' errore nella shape'])
#      end
#
#     elseif j==1
#      if i==1
#        val=-1.;
#      elseif i==2
#        val=1.;
#      elseif i==3
#        val=0.;
#      else
#        error(['j=',j,' errore nella shape'])
#      end
#
#      elseif j==2
#      if i==1
#        val=-1.;
#      elseif i==2
#        val=0.;
#      elseif i==3
#        val=1.;
#      else
#        error(['j=',j,' errore nella shape'])
#      end
#
#      else
#        print *,'j=',j,' errore nella shape'
#      end

