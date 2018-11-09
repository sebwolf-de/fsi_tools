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
f.close()

z_u = np.zeros(x_u.shape)
n_triangles_u = topo_u.shape[0]
offset_u = (np.ones(n_triangles_u)*3).cumsum().astype(int)
topo_u = topo_u.reshape((n_triangles_u*3))
ctype_u = np.ones(n_triangles_u)*VtkTriangle.tid

topo_p_old = topo_p
n_triangles_p = topo_p.shape[0]
x_p_mod = np.zeros((n_triangles_p*3))
y_p_mod = np.zeros((n_triangles_p*3))
topo_p_mod = np.zeros((topo_p.shape[0], 3))
k = 0
for row in topo_p:
     x_l = x_p[row[0:3]]
     y_l = y_p[row[0:3]]
     x_p_mod[3*k:3*k+3] = x_l
     y_p_mod[3*k:3*k+3] = y_l
     topo_p_mod[k,:] = [3*k, 3*k+1, 3*k+2]
     k = k+1
x_p = x_p_mod
y_p = y_p_mod
topo_p = topo_p_mod
z_p = np.zeros(x_p.shape)
offset_p = (np.ones(n_triangles_p)*3).cumsum().astype(int)
topo_p = topo_p[:,(0,1,2)].reshape((n_triangles_p*3))
ctype_p = np.ones(n_triangles_p)*VtkTriangle.tid

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

    p_vertex = np.zeros(x_p.shape)
    k = 0
    for row in topo_p_old:
        p_vertex[3*k:3*k+3] = np.reshape(p[row[0:3]] + p[row[3]], (3))
        k = k+1

    unstructuredGridToVTK(vtk_dir+'pressure_cn_time_'+str(cn_time).zfill(ph.time_index_digits), x_p, y_p, z_p, connectivity = topo_p, offsets = offset_p, cell_types = ctype_p, cellData = None, pointData = {"p": p_vertex})
    unstructuredGridToVTK(vtk_dir+'fluid_cn_time_'+str(cn_time).zfill(ph.time_index_digits), x_u, y_u, z_u, connectivity = topo_u, offsets = offset_u, cell_types = ctype_u, cellData = None, pointData = {"vel": (u_x, u_y, z_u)})
