#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os
from pyevtk.hl import unstructuredGridToVTK
from pyevtk.vtk import VtkTriangle
import scipy.sparse as sparse
import sys

sys.path.append('../../modules')
import basis_func
from parameters_handler import ParametersHandler

if len(sys.argv) > 1:
    ph = ParametersHandler(sys.argv[1])
else:
    ph = ParametersHandler('simulation_parameters.json')
ph.simulation_info()

time_list = np.where(ph.stampa)[0] + 1
time_list = np.append(np.array([0]), time_list)

results_dir = ph.results_directory+'/'+ph.sim_prefix+'/binary_data/'

filename = results_dir+'mesh'
f = open(filename,"rb")
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
t_lgr = np.load(f)
f.close()

topo_u_old = topo_u
z_u = np.zeros(x_u.shape)
n_triangles_u = topo_u.shape[0]
offset_u = (np.ones(n_triangles_u)*3).cumsum().astype(int)
topo_u = topo_u.reshape((n_triangles_u*3))
ctype_u = np.ones(n_triangles_u)*VtkTriangle.tid

mat_u_dx = sparse.lil_matrix((n_triangles_u, x_u.shape[0]))
mat_u_dy = sparse.lil_matrix((n_triangles_u, x_u.shape[0]))
k = 0
for row in topo_u_old:
    x_l = x_u[row]
    y_l = y_u[row]
    eval_points = np.array([[np.sum(x_l)/3, np.sum(y_l)/3]])
    (phi_dx,phi_dy,phi,omega) = basis_func.tri_p1(x_l,y_l,eval_points)
    if omega > 1e-8:
        mat_u_dx[k, row] = phi_dx
        mat_u_dy[k, row] = phi_dy
    k = k+1

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

topo_s_old = topo_s
z_s = np.zeros(x_s.shape)
n_triangles_s = topo_s.shape[0]
offset_s = (np.ones(n_triangles_s)*3).cumsum().astype(int)
topo_s = topo_s.reshape((n_triangles_s*3))
ctype_s = np.ones(n_triangles_s)*VtkTriangle.tid

mat_s_dx = sparse.lil_matrix((n_triangles_s, x_s.shape[0]))
mat_s_dy = sparse.lil_matrix((n_triangles_s, x_s.shape[0]))
k = 0
for row in topo_s_old:
    x_l = x_s[row]
    y_l = y_s[row]
    eval_points = np.array([[np.sum(x_l)/3, np.sum(y_l)/3]])
    (phi_dx,phi_dy,phi,omega) = basis_func.tri_p1(x_l,y_l,eval_points)
    if omega > 1e-8:
        mat_s_dx[k, row] = phi_dx
        mat_s_dy[k, row] = phi_dy
    k = k+1

vtk_dir = ph.results_directory+'/'+ph.sim_prefix+'/vtk/'
if not os.path.exists(vtk_dir):
    os.makedirs(vtk_dir)

for cn_time in time_list:
    input_name = results_dir+'cn_time_'+str(cn_time).zfill(ph.time_index_digits)
    print(input_name)
    f = open(input_name,"rb")
    u = np.load(f)
    p = np.load(f)
    s_x = np.load(f) - s_lgr
    s_y = np.load(f) - t_lgr
    f.close()

    ndofs = int(u.shape[0]/2)
    u_x = u[0:ndofs]
    u_y = u[ndofs:2*ndofs]

    p_vertex = np.zeros(x_p.shape)
    k = 0
    for row in topo_p_old:
         p_vertex[3*k:3*k+3] = p[row[0:3]] + p[row[3]]
         k = k+1

    u_dx = mat_u_dx.dot(u_x)
    v_dx = mat_u_dx.dot(u_y)
    u_dy = mat_u_dy.dot(u_x)
    v_dy = mat_u_dy.dot(u_y)

    s_dx = mat_s_dx.dot(s_x)
    t_dx = mat_s_dx.dot(s_y)
    s_dy = mat_s_dy.dot(s_x)
    t_dy = mat_s_dy.dot(s_y)

    norm_deform_grad = np.sqrt(s_dx**2+s_dy**2+t_dx**2+t_dy**2)

    unstructuredGridToVTK(vtk_dir+'pressure_cn_time_'+str(cn_time).zfill(ph.time_index_digits), x_p, y_p, z_p, connectivity = topo_p, offsets = offset_p, cell_types = ctype_p, cellData = None, pointData = {"p": p_vertex})
    unstructuredGridToVTK(vtk_dir+'fluid_cn_time_'+str(cn_time).zfill(ph.time_index_digits), x_u, y_u, z_u, connectivity = topo_u, offsets = offset_u, cell_types = ctype_u, cellData = {"u_dx": u_dx, "v_dx": v_dx, "u_dy": u_dy, "v_dy": v_dy}, pointData = {"vel": (u_x, u_y, z_u)})
    # unstructuredGridToVTK(vtk_dir+'fluid_cn_time_'+str(cn_time).zfill(ph.time_index_digits), x_u, y_u, z_u, connectivity = topo_u, offsets = offset_u, cell_types = ctype_u, cellData = None, pointData = {"vel": (u_x, u_y, z_u)})
    # unstructuredGridToVTK(vtk_dir+'structure_cn_time_'+str(cn_time).zfill(ph.time_index_digits), s_lgr, t_lgr, z_s, connectivity = topo_s, offsets = offset_s, cell_types = ctype_s, cellData = None, pointData = {"ds": (s_x, s_y, z_s)})
    unstructuredGridToVTK(vtk_dir+'structure_cn_time_'+str(cn_time).zfill(ph.time_index_digits), s_lgr, t_lgr, z_s, connectivity = topo_s, offsets = offset_s, cell_types = ctype_s, cellData = {"s_dx": s_dx, "t_dx": t_dx, "s_dy": s_dy, "t_dy": t_dy, "norm_deform_grad": norm_deform_grad}, pointData = {"ds": (s_x, s_y, z_s)})
