# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 15:03:47 2015

@author: nicola
"""
import numpy as np

class Mesh:
    def __init__(self,file_to_be_read):
        self.filename = file_to_be_read
#        self.topo = np.array([[]])
#        self.x = np.array([])
#        self.y = np.array([])
#        self.nodes = np.array([])
#        self.b_nodes = np.array([])        
        return
        
    def read_data_from_input(self):
        f = open ( self.filename , 'r')
        self.topo = np.array([[]])
        self.x = np.array([])
        self.y = np.array([])
        self.nodes = np.array([])
        self.b_nodes = np.array([])
        for line in f: 
            if line[0]!='$':
#                print 'non fare un cippa'
#            else:
                l = map(float,line.split())

                if len(l) == 4:
                    self.x = np.append(self.x,l[1])
                    self.y = np.append(self.y,l[2])
                    self.nodes = np.append(self.nodes,int(l[0])-1)                

                if len(l) > 4:
                    if l[1]==1 or l[1]==8: 
                        nod = l[5:]
                        print l[1] ,nod
                        for i in nod:
                            i= int(i-1)
                            if i not in self.b_nodes:  
                                self.b_nodes = np.append( self.b_nodes , i )
                    
                    if l[1] == 2 or l[1]==9:
                        row = l[5:]
                        row = np.array(row, dtype=int)
                        self.topo = np.append(self.topo,row)
                        
        self.topo = np.reshape(self.topo,(len(self.topo)/len(row),len(row)))
        self.topo = self.topo-1
        r_id = 0 
        for row in self.topo:
            ck =      (self.x[row[1]]-self.x[row[0]])*(self.y[row[2]]-self.y[row[0]])
            ck = ck - (self.x[row[2]]-self.x[row[0]])*(self.y[row[1]]-self.y[row[0]])
            if ck < 0:
                self.topo[r_id,:] = np.array([[row[0],row[2],row[1]]])
            r_id+=1        
        print r_id    
        
        return 
        
        