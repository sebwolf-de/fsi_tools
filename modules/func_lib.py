# -*- coding: utf-8 -*-
"""A Library of Functions.

.. moduleauthor:: Nicola Cavallini

"""

class Parabola:
    """
    This a first example on how to program a Parabola as a class.
    The classical formulation for a parabola is:
        
    .. math::
       y = a\cdot x^2 + b \cdot x + c
    """
    def evaluate(self,x):
        """
        Usually we want to evaluate the Function, then a good idea is the 
        evaluate function.
    .. code:: python
    
        def evaluate(self,x):        
            y = self.a* x**2+ self.b*x+self.c
            return y
            
    First of all notice that an the first argument of 
    a class function, is the class itself, the word 
    ``self``. This is a little peculiar but that's 
    the language, I think that is done to prevent people 
    coding super complicated classes.
    
    Adding ``self.`` in front of a variable makes it a 
    a public variable. Typing ``self.a``, ``self.b``, ``self.c``, 
    we defined three public variables corresponding to 
    the parabola coefficients. A public variable is a variable that 
    can be accessed and modified from outside the class.
        
        """
        
        y = self.a* x**2+ self.b*x+self.c
        return y