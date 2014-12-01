The Laplacian Example
=====================

The usual way to learn about finite elements 
is to fist solve the Laplace operator. This is 
simple enough to be understood, but it encapsulates most of the 
difficulties in coding finite elements.\

Our chosen environment is Python. One day I'll tell you why.\

It is now time to get our hands dirty, with some code. 
Unless you want to include all the code you need in a single file, 
resulting in a super confusing mess, you will need to ``include``
some modules that some nice guy prepared for you. This is 
exactly what you do at the beginning of your file. This can be 
standard libraries, like: 

.. literalinclude:: ../../tests/ex_laplace.py
   :lines: 7-8


Or some pieces of code included in this non-library:

.. literalinclude:: ../../tests/ex_laplace.py
   :lines: 10-14

In this lines we have respectively:

``lin_tri_mesh``:
  a module dedicated to mesh generation and reading.

``basis_func``:
  a module dedicated to basis functions.

``assemble``:
  a collection of assembling routines.

``la_utils``:
  small and useful linear algebra tricks.

``viewers``:
  functions useful for plotting.

Now we collect all of these functions into the ``__main__`` section of a script,
and we execute them one after another. This strange line tells the interpreter 
that this is the main section of the program.

.. literalinclude:: ../../tests/ex_laplace.py
   :lines: 16

Now we define the discretization parameters and we generate the mesh. 
The way I usually represent a mesh is with the ``x`` and ``y`` nodes coordinates, 
and with the ``topo``-ogical matrix (a.k.a. connectivity matrix). Each row correspond 
to an element. In each row the corresponding nodes are stored. (At some point mesh will become a class).

.. literalinclude:: ../../tests/ex_laplace.py
   :lines: 18-24

The "stiffness" matrix A, is assembled in the ``assemble`` module. We redirect to this
module documentation for details. And we actually recommend to check it, because people
in scientific computing usually want to control this section of their codes.

.. literalinclude:: ../../tests/ex_laplace.py
   :lines: 26

The ``rhs`` implementation is a little rough. I take advantage of having a structured mesh and I 
evaluate the element only once. To do this I use the ``tri_p1`` function then I loop over all 
the elements adding the contribution of a unitary function (our forcing) on each element.

.. literalinclude:: ../../tests/ex_laplace.py
   :lines: 30-37

Now we intend to set the boundary conditions. In this case we delete the A rows 
corresponding to fixed dofs, leaving only the diagonal entries to their original 
value. We locate the dofs corresponding to the domain border using the ``numpy`` 
function ``where``. 

.. literalinclude:: ../../tests/ex_laplace.py
   :lines: 39-41

Once we repeated this procedure for every side of the square 
we only need to solve the linear system. SciPy is the guy that helps us 
in this concern:

.. literalinclude:: ../../tests/ex_laplace.py
   :lines: 56

We only now need to plot the solution using one of the prepared functions in 
viewer:

.. literalinclude:: ../../tests/ex_laplace.py
   :lines: 58

This tutorial doesn't make much sense without a series of exercises 
that will help to get an acquittance with both Python utilities and the
buggy code I have wrote.

   1. Plot, using the ``matplotlib`` command ``spy`` the sparsity pattern for A.

   2. Play around with the mesh parameters, ultimately change from a structured mesh to an unstructured one.  ``gmsh`` is the tool to generate the mesh. Inside ``lin_tri_mesh`` there is a tool to read the mesh.

   3. Change the forcing to a manufactured solution.

   4. Evaluate the error and edit a convergence table. In the ``io_tools`` there should be something useful.

   5. Change the way to impose boundary conditions. Not only delete A rows, but also the columns (still keeping the diagonal entries unchanged). Then add the "lift" contribution to the right han side.

   6. Now you should have a symmetric matrix. Check it with ``spy`` command. It is time to change your solver to a conjugate gradient.

   7. This is difficult: try to implement null and ideal preconditioners for this linear system... Don't try this at home.

The plain code follows: \


.. literalinclude:: ../../tests/ex_laplace.py

