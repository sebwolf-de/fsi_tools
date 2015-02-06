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
    k=0    
    for i in r:
        if  ( i > 1e-5 ):
            teta = np.arcsin(y[k]/i)
            if  (x[k] >= 0):                                   
                esa_tmp = (i**(2./3.))*np.cos (2*teta/3)  
                esa = np.append(esa,esa_tmp) 
            if  ( x[k] < 0 )&( y[k] > 0) :
                esa_tmp = (i**(2./3.))*np.cos (2*(np.pi-teta)/3)
                esa = np.append(esa,esa_tmp)
            if  ( x[k] < 0 )&( y[k] < 0) :
                esa_tmp = (i**(2./3.))*np.cos (2*(-teta-np.pi)/3)
                esa = np.append(esa,esa_tmp)

        else:
            esa = np.append(esa,0)
        k=k+1
    return esa       

def der_sol_esatta_1(x,y):
    r_x = np.array([])
    r_y = np.array([])
#    teta = np.array([])
    teta_x = np.array([])
    teta_y = np.array([])
    esa_r = np.array([])
    esa_teta = np.array([])
    esa_x = np.array([])
    esa_y = np.array([])
    r = np.sqrt(x**2+y**2)
#    print r.shape
    k=0    
    for i in r:            
        if  ( i > 1e-5 ):
            a = x[k]/i
            b = y[k]/i 
            r_x = np.append(r_x,a)
            r_y = np.append(r_y,b)            
            teta = np.arcsin(y[k]/i)
#            print "teta",teta.shape
#            teta = np.append( teta, alfa )
            teta1 = np.array ([-x[k]*y[k]/((np.sqrt(i**2+y[k]**2))*(i**2))])
            teta2 = np.array ([(i**2-y[k]**2)/((np.sqrt(i**2+y[k]**2))*(i**2))])        
            teta_x = np.append(teta_x,teta1)
            teta_y = np.append(teta_y,teta2)
#            print "teta_x",teta_x.shape
#            print "teta_y",teta_y.shape
            if  (x[k] >= 0):                                   
                esa_tmp_r = (2./3.)*(i**(-1./3.))*np.cos (2*teta/3)
                esa_tmp_teta = -i*(2./3.)*(2./3.)*np.sin (2*teta/3)                
                                
                esa_r = np.append(esa_r,esa_tmp_r) 
                esa_teta = np.append(esa_teta,esa_tmp_teta)
#                print "esa_teta", esa_teta.shape
            if  ( x[k] < 0 )&( y[k] > 0) :
                esa_tmp_r = (2./3.)*(i**(-1./3.))*np.cos (2*(np.pi-teta)/3)
                esa_tmp_teta = -i*(2./3.)*(2./3.)*np.sin (2*(np.pi-teta)/3)                
                esa_r = np.append(esa_r,esa_tmp_r)
                esa_teta = np.append(esa_teta,esa_tmp_teta)
            if  ( x[k] < 0 )&( y[k] < 0) :
                esa_tmp_r = (2./3.)*(i**(-1./3.))*np.cos (2*(-teta-np.pi)/3)
                esa_tmp_teta = -i*(2./3.)*(2./3.)*np.sin (2*(-teta-np.pi)/3)                
                esa_r = np.append(esa_r,esa_tmp_r)
                esa_teta = np.append(esa_teta,esa_tmp_teta)
#            esa_x = esa_r*r_x + esa_teta*teta_x
#            esa_y = esa_r*r_y + esa_teta*teta_y
#            print "esa_x",esa_x.shape
#            print "esa_y",esa_y.shape 
        else:
            esa_r = np.append(esa_r,0)
            esa_teta = np.append(esa_teta,0)
            r_x =  np.append(r_x,0)
            r_y =  np.append(r_y,0)
            teta_x = np.append(teta_x,0)
            teta_y = np.append(teta_y,0)
        k=k+1
        
#    print "esa_teta", esa_teta.shape
#    print "esa_r" , esa_r.shape
#    print "r_x" , r_x.shape
#    print "r_y" , r_y.shape

    esa_x = esa_r*r_x + esa_teta*teta_x
    esa_y = esa_r*r_y + esa_teta*teta_y
#    

#    esa_x = esa_r*r_x + esa_teta*teta_x
#    esa_y = esa_r*r_y + esa_teta*teta_y
    return  esa_x, esa_y
#    r_x = x./r
#    r_y = np.array([y./r])    
#    teta_x = np.array ([-x*y/((np.sqrt(r**2+y**2))*(r**2))])
#    teta_y = np.array ([(r**2-y**2)/((np.sqrt(r**2+y**2))*(r**2))]) 
#    J = np.array([[x./r, y./r],[-x*y/((np.sqrt(r**2+y**2))*(r**2)),(r**2-y**2)/((np.sqrt(r**2+y**2))*(r**2))]])

    