#! /usr/bin/env python

import os
import sys
import numpy as np

from pyevtk.hl import unstructuredGridToVTK
from pyevtk.vtk import VtkTriangle
x_u = np.array([0., 1., 0., 1., 1., 0.])
y_u = np.array([0., 0., 1., 0., 1., 1.])
topo_u = np.array([[0, 1, 2], [3, 4, 5]])
z_u = np.zeros(6)
n_triangles_u = 2
offset_u = (np.ones(n_triangles_u)*3).cumsum().astype(int)
topo_u = topo_u.reshape((n_triangles_u*3))
ctype_u = np.ones(n_triangles_u)*VtkTriangle.tid


unstructuredGridToVTK('test', x_u, y_u, z_u, connectivity = topo_u, offsets = offset_u, cell_types = ctype_u, cellData = None, pointData = {"ds": np.array([0, 0, 0, 3, 3, 3])})
