``ex_mesh.py`` Generating and Reading a Mesh
============================================

In this tutorial we generate and read a mesh. A standard tool 
for mesh generation is `gmsh`_. The web page contains a lot of documentation 
and tutorials. Also video tutorials... very clear.

.. _gmsh: http://geuz.org/gmsh/

In the directory ``fsi_tools/gmsh_apps/`` there are three ``.geo`` files. 
They are typical input files for ``gmsh``. I report here three different ways of meshing a 
square. In the first one we generate a point at the origin ``Point(1) = {0, 0.0, 0, cl1};``. The 
variable ``cl1`` is not useful at the moment. The point is ``Extrude {1, 0, 0} {...}`` along the 
fist axis for a length 1, obtaining a line. ``lyx`` specifies the number of layers you need to subdivide
the mesh. The line is extruded for a length 1, obtaining a square.

.. literalinclude:: ../../gmsh_apps/structerd.geo

Analogous is the the case for a symmetric mesh. The only difference is that we extrude the 
line in two opposite directions, obtaining the symmetry of the triangles. 

.. literalinclude:: ../../gmsh_apps/structerd_sym.geo

A completely unstructured mesh from the ``gmsh`` tutorial. ``ms`` is the mesh size.

.. literalinclude:: ../../gmsh_apps/unstr_square.geo

Calling the Mesh Generator
--------------------------

Invoke ``gmsh`` typing ``gmsh`` from command the terminal. ``gmsh`` has a nice 
graphical interface.

* Open the file: File->Open->Filename.
* Generate the mesh: Mesh->1D->2D.
* Save the mesh: File->Save Mesh. This will generate the ``.msh`` file, explore it!

Meshing in our Code
-------------------
In this tutorial we both generate the mesh using our hard coded commands,
and read ``.msh`` files. Both these approaches are pretty simple. They are implemented in
the ``lin_tri_mesh`` module, that we are importing giving it a shorter name: ``lin_t3``. 
This is the meaning of ``as`` after import. ``viewers`` module collects number of 
functions dedicated to visualization.

This special line separates the *header* section of the script from the *main*
section.

.. literalinclude:: ../../tests/ex_mesh_read.py
   :lines: 9


We first use our hard coded function, the classical and not so smart "by hand approach".

.. literalinclude:: ../../tests/ex_mesh_read.py
   :lines: 11-15

In the following lines we are smater and include the ``gmsh`` generated meshes. 
Filename is a string as we can see in the plain code.

.. code:: python

   (topo,x,y) = lin_t3.load_msh(filename)

We use the simple implemented viewer to plot the mesh. Check the 
`source!`_
Once you checked the documented source, go to real source, to see 
how thins look like. *This is important because it shows how to 
iterate through a mesh.*

.. _source!: ./viewers.html

.. code:: python

   viewers.tri_plot(x,y,topo)


Exercises
---------

This tutorial doesn't make much sense without a series of exercises 
that will help to get an acquittance with both Python utilities and the
buggy code I have wrote.


1. Define a regular mesh with different numbers of elements along the edges.
2. Define an unstructured mesh with different mesh sizes at the corners.
3. Mesh a circle.
4. Now try to iterate over a mesh. And print on the screen the dofs coords.

The plain code follows: \


.. literalinclude:: ../../tests/ex_mesh_read.py

