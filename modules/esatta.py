import numpy as np

def load(x,y):
      f = 2*(x*(1-x)+y*(1-y))
      return f

def sol_esatta(x,y):
	esa = x*(1-x)*y*(1-y)
	return esa

def dx_sol_esatta(x,y):
	dx_esa = (1-2*x)*y*(1-y)
	return dx_esa

def dy_sol_esatta(x,y):
	dy_esa = (1-2*y)*x*(1-x)
	return dy_esa

#===========================================================================

def load1(x,y):
      f = np.zeros(x.shape)
      return f


def sol_esatta1(x,y):
    esa = np.array([])
    r = np.sqrt(x**2+y**2)
#    print r
    k=0    
    for i in r:
        if  ( i > 1e-5 ):
            teta = np.arcsin(y[k]/i)
#            print 'teta'
#            print teta
            if  (x[k] >= 0):                                   
                esa_tmp = (i**2./3.)*np.cos (2*teta/3)  
                esa = np.append(esa,esa_tmp) 
#                print esa , '1'
            if  ( x[k] < 0 )&( y[k] > 0) :
                esa_tmp = (i**2./3.)*np.cos (2*(np.pi-teta)/3)
                esa = np.append(esa,esa_tmp)
#                print esa , '2'
            if  ( x[k] < 0 )&( y[k] < 0) :
                esa_tmp = (i**2./3.)*np.cos (2*(-teta-np.pi)/3)
                esa = np.append(esa,esa_tmp)

        else :
            esa = np.append(esa,0)
#            print esa , '3'
            #       # 
#           # b = np.where ( (x < 0) & (r>1e-5))
#            bc[j] = (i**2/3)*np.cos (-2*teta + np.pi/4)    #ok
#        #d = np.where ((x == 0) & ( y<0 ) )
        k=k+1
    return esa       
    