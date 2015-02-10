# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 11:05:15 2015

@author: michele
"""

import numpy as np


def load_msh(filename):
    f = open ( filename , 'r')
    x = np.array([])
    y = np.array([])
    topo = np.array([], dtype=int)
    nodes = np.array([], dtype=int)
    b_nodes = np.array([], dtype=int)
    b_nodes_1 = np.array([], dtype=int)
    b_nodes_2 = np.array([], dtype=int)
    int_nodes = np.array([], dtype=int)
    physical_group = np.array([], dtype=int)
#    for line in g:
#        if     

    for line in f:
        
        if line[0]=='$':
            print 'non fare un cippa'
        else:
            l = map(float,line.split())
            print l
            if len(l) == 4:
                x = np.append(x,l[1])
                y = np.append(y,l[2])
                nodes = np.append(nodes,int(l[0])-1)                
                #print 'ciao'
            if len(l) == 7:
                physical_group = np.append( physical_group  , l[3] )
                #print physical_group
                nod = l[5:7]
                for i in nod:
                    i= int(i-1)
                    if i not in b_nodes:  
                        b_nodes = np.append( b_nodes , i )
 
            if len(l) == 8:
                row = l[5:8]
                row = np.array(row, dtype=int)
                topo = np.append(topo,row)
    
    for  j in nodes:
            if j not in b_nodes:
                if j not in int_nodes:
                    int_nodes = np.append ( int_nodes , j )    
    
    j=0
    for k in b_nodes:
        if physical_group [j] == physical_group [0]:
             b_nodes_1 = np.append( b_nodes_1 , k )
        else:
            b_nodes_2 = np.append( b_nodes_2 , k )
        j=j+1
    
    topo = np.reshape(topo,(len(topo)/3,3))
    topo = topo-1
    r_id = 0 
    for row in topo:
        ck =      (x[row[1]]-x[row[0]])*(y[row[2]]-y[row[0]])
        ck = ck - (x[row[2]]-x[row[0]])*(y[row[1]]-y[row[0]])
        if ck < 0:
            topo[r_id,:] = np.array([[row[0],row[2],row[1]]])
        r_id+=1        
    print r_id    
    
    return topo, x, y, b_nodes, b_nodes_1, b_nodes_2, int_nodes   
    
    
    