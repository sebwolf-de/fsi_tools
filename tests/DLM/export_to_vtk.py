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
    ph = ParametersHandler('simulation_parameters.json')
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
topo_s = np.load(f)
x_s = np.load(f)
y_s = np.load(f)
s_lgr = np.load(f)
f.close()

z_u = np.zeros(x_u.shape)
n_triangles_u = topo_u.shape[0]
offset_u = (np.ones(n_triangles_u)*3).cumsum().astype(int)
topo_u = topo_u.reshape((n_triangles_u*3))
ctype_u = np.ones(n_triangles_u)*VtkTriangle.tid

z_s = np.zeros(x_s.shape)
n_triangles_s = topo_s.shape[0]
offset_s = (np.ones(n_triangles_s)*3).cumsum().astype(int)
topo_s = topo_s.reshape((n_triangles_s*3))
ctype_s = np.ones(n_triangles_s)*VtkTriangle.tid

vtk_dir = ph.results_directory+'/'+ph.sim_prefix+'/vtk/'
if not os.path.exists(vtk_dir):
    os.makedirs(vtk_dir)

for cn_time in time_list:
    input_name = results_dir+'cn_time_'+str(cn_time).zfill(ph.time_index_digits)
    print input_name
    f = file(input_name,"rb")
    u = np.load(f)
    p = np.load(f)
    s_x = np.load(f) - x_s
    s_y = np.load(f) - y_s
    f.close()

    ndofs = u.shape[0]/2
    u_x = u[0:ndofs]
    u_y = u[ndofs:2*ndofs]

    unstructuredGridToVTK(vtk_dir+'fluid_cn_time_'+str(cn_time).zfill(ph.time_index_digits), x_u, y_u, z_u, connectivity = topo_u, offsets = offset_u, cell_types = ctype_u, cellData = None, pointData = {"vel": (u_x, u_y, z_u)})
    unstructuredGridToVTK(vtk_dir+'structure_cn_time_'+str(cn_time).zfill(ph.time_index_digits), x_s, y_s, z_s, connectivity = topo_s, offsets = offset_s, cell_types = ctype_s, cellData = None, pointData = {"ds": (s_x, s_y, z_s)})
