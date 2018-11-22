# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 15:03:47 2015

@author: nicola
"""
import numpy as np
#from reference_shape_function import ReferenceShapeFunction
#from quadrature import Quadrature

class Mesh:
    def __init__(self,file_to_be_read):
        self.filename = file_to_be_read     
        self.el_id = 0
        return
    
#    def next(self):
#        self.el_id +=1
#        return
        
    def read_data_from_input(self):
        f = open ( self.filename , 'r')
        self.topo = np.array([[]], dtype=int)
        self.x = np.array([])
        self.y = np.array([])
        self.nodes = np.array([], dtype=int)
        self.b_nodes = np.array([], dtype=int)
        self.row = np.array([], dtype=int)         
        for line in f: 
            if line[0]!='$':
#                print 'non fare un cippa'
#            else:
                l = list(map(float,line.split()))

                if len(l) == 4:
                    self.x = np.append(self.x,l[1])
                    self.y = np.append(self.y,l[2])
                    self.nodes = np.append(self.nodes,int(l[0])-1)                

                if len(l) > 4:
                    if l[1]==1 or l[1]==8: 
                        nod = l[5:]
#                        print l[1] ,nod
                        for i in nod:
                            i= int(i-1)
                            if i not in self.b_nodes:  
                                self.b_nodes = np.append( self.b_nodes , i )
                    
                    if l[1] == 2 or l[1]==9:
                        row = l[5:]
                        self.row = np.array(row, dtype=int)-1
                        self.topo = np.append(self.topo,self.row)
                        

#                    print self.topo
#                        print self.row
#                        
#        self.topo = np.reshape(self.topo,len(self.topo)/6,6)
#        self.topo = self.topo-1        
        self.topo = np.reshape(self.topo,(len(self.topo)/len(self.row),len(self.row)))
        r_id = 0
        for self.row in self.topo:
            ck =      (self.x[self.row[1]]-self.x[self.row[0]])*(self.y[self.row[2]]-self.y[self.row[0]])
            ck = ck - (self.x[self.row[2]]-self.x[self.row[0]])*(self.y[self.row[1]]-self.y[self.row[0]])
            if ck < 0:
                self.topo[r_id,:] = np.array([[self.row[1],self.row[0],self.row[2],self.row[3],self.row[5],self.row[4]]])
            r_id+=1        
#        print r_id
        return


#        self.topo = self.topo-1
#        r_id = 0 
#
#        
#        for self.row in self.topo:
#            ck =      (self.x[self.row[1]]-self.x[self.row[0]])*(self.y[self.row[2]]-self.y[self.row[0]])
#            ck = ck - (self.x[self.row[2]]-self.x[self.row[0]])*(self.y[self.row[1]]-self.y[self.row[0]])
#            if ck < 0:
#                self.topo[r_id,:] = np.array([[self.row[0],self.row[2],self.row[1]]])
#            r_id+=1        
##        print r_id
#        return 
    
    

        
        