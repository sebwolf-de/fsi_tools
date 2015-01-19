``ex_classes.py`` Example of a Simple Class
===========================================

We first import the modules we need. Additionally, from 
the module ``func_lib.py`` we import our `Parabola`_ class.

.. _Parabola: ./func_lib.html


.. literalinclude:: ../../tests/ex_class.py
   :lines: 1-4


We move to the main section of the program, and build an 
abscissa using numpy.

.. literalinclude:: ../../tests/ex_class.py
   :lines: 6-8

We instantiate an object ``par`` of type ``Parabola``:

.. literalinclude:: ../../tests/ex_class.py
   :lines: 10

We access the public variables of the object ``par``,
and assign them a numerical value. 

.. literalinclude:: ../../tests/ex_class.py
   :lines: 12-15

I use the ``evaluate`` method to evaluate 
the parabola named ``par``. Then we plot its values.

.. literalinclude:: ../../tests/ex_class.py
   :lines: 16-18

I define ``another_parabola`` and I assign 
different values to the coefficients, I evaluate values 
for the new parabola and plot them.

.. literalinclude:: ../../tests/ex_class.py
   :lines: 20-28

Now I access the ``c`` variable of the 
``par`` parabola and modify it. I evaluate again
the parabola values and plot them.

.. literalinclude:: ../../tests/ex_class.py
   :lines: 30-34

The full code.

.. literalinclude:: ../../tests/ex_class.py


Exercises
---------

1. Add the function ``plot`` to the class 
``Parabola``.

2. Find a smart way to instantiate your class.
instead of:

.. code:: python

   par = Parabola()
   par.a = 2
   par.b = 0
   par.c = 0

You should get a more elegant:

.. code:: python

   par = Parabola(2,0,0)

This is done implementing a constructor. Practically
you need to add this function to your ``Parabola`` class:

.. code:: python

   def __init__(self,a,b,c,):
       self.a = a
       self.b = b
       self.c = c
       return

A constructor is a particular function that 
initializes a given class specifying some 
basic data.
