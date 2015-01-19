import numpy as np
import matplotlib.pyplot as plt

from func_lib import Parabola

if __name__== '__main__':
    
    x = np.linspace(0,3)
    
    par = Parabola()
    
    par.a = 2
    par.b = 0
    par.c = 0

    y = par.evaluate(x)
    
    plt.plot(x,y,'b')
    
    another_parabola = Parabola()
    
    another_parabola.a = 3
    another_parabola.b = 0
    another_parabola.c = 0
    
    y = another_parabola.evaluate(x)
    
    plt.plot(x,y,'r')
    
    par.c = 1
    
    y = par.evaluate(x)
    
    plt.plot(x,y,'g')
    
    plt.show()