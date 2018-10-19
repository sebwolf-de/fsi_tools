#! /usr/bin/env python

import os
import sys
import numpy as np

from pyevtk.hl import unstructuredGridToVTK
from pyevtk.vtk import VtkTriangle

from parameters_handler import ParametersHandler

if len(sys.argv) > 1:
    ph = ParametersHandler(sys.argv[1])
else:
    ph = ParametersHandler('simulation_parameters_fsi.json')
ph.simulation_info()

time_list = np.where(ph.stampa)[0]

results_dir = ph.results_directory+'/'+ph.sim_prefix+'/binary_data/'

filename = results_dir+'mesh'
f = file(filename,"rb")
topo_p = np.load(f)
x_p = np.load(f)
y_p = np.load(f)
topo_u = np.load(f)
x_u = np.load(f)
y_u = np.load(f)
c2f = np.load(f)
f.close()

z_u = np.zeros(x_u.shape)
n_triangles = topo_u.shape[0]
offset = (np.ones(n_triangles)*3).cumsum().astype(int)
topo_u = topo_u.reshape((n_triangles*3))
ctype = np.ones(n_triangles)*VtkTriangle.tid

vtk_dir = ph.results_directory+'/'+ph.sim_prefix+'/vtk/'
if not os.path.exists(vtk_dir):
    os.makedirs(vtk_dir)

for cn_time in time_list:
    input_name = results_dir+'cn_time_'+str(cn_time).zfill(ph.time_index_digits)
    print input_name
    f = file(input_name,"rb")
    u = np.load(f)
    p = np.load(f)
    f.close()

    ndofs = u.shape[0]/2
    u_x = u[0:ndofs]
    u_y = u[ndofs:2*ndofs]

    unstructuredGridToVTK(vtk_dir+'cn_time_'+str(cn_time).zfill(ph.time_index_digits), x_u, y_u, z_u, connectivity = topo_u, offsets = offset, cell_types = ctype, cellData = None, pointData = {"vel": (u_x, u_y, z_u)})
