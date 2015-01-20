import numpy as np


from func_lib import Parabola

if __name__== '__main__':
    
    x = np.linspace(0,3)
    
    par = Parabola(2,0,0)
    graph=par.plot(x,'r')    
        
    another_parabola = Parabola(3,0,1)
    graph=another_parabola.plot(x,'b')
    
    