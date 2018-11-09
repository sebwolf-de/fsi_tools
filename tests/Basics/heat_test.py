#! /usr/bin/env python

import numpy as np
import os
from pyevtk.hl import unstructuredGridToVTK
from pyevtk.vtk import VtkTriangle
import scipy.sparse.linalg as sp_la

import assemble
import la_utils
import lin_tri_mesh
import viewers

### problem parameters
dt = 0.1
dx = 0.1
nx = 10
kappa = 0.01

### geometry setup
(topo_p, x_p, y_p, topo_u, x_u, y_u, c2f) = \
    lin_tri_mesh.load_t3_iso_t6_file('../mesh_collection/step.msh', '../mesh_collection/step_refined.msh')
viewers.tri_plot(x_u, y_u, topo_u)
viewers.tri_plot(x_p, y_p, topo_p)

### matrix assembly
mass = assemble.u_v_p1(topo, x, y)
stiffness = assemble.gradu_gradv_p1(topo, x, y)
BDF1 = 1/dt*mass + kappa*stiffness

### initial condition
u_n = np.zeros(x.shape)

### boundary conditions
bc_id = np.where(x < dx/10)
BDF1 = la_utils.set_diag(BDF1, bc_id)
bc_id = np.where(y < dx/10)
BDF1 = la_utils.set_diag(BDF1, bc_id)
bc_id = np.where(x > 1 - dx/10)
BDF1 = la_utils.set_diag(BDF1, bc_id)
bc_id = np.where(y > 1 - dx/10)
BDF1 = la_utils.set_diag(BDF1, bc_id)

f = np.zeros(u_n.shape)

### VTK setup
z = np.zeros(x.shape)
n_triangles = topo.shape[0]
offset = (np.ones(n_triangles)*3).cumsum().astype(int)
topo = topo.reshape((n_triangles*3))
ctype = np.ones(n_triangles)*VtkTriangle.tid

vtk_dir = 'heat_results'
if not os.path.exists(vtk_dir):
    os.makedirs(vtk_dir)


### time loop
for k in range(0, 250):
    ### assemble right hand side
    f = 1/dt * mass.dot(u_n)
    bc_id = np.where(x < dx/10)
    f[bc_id] = 100*np.sin(np.pi*y[bc_id])**2
    bc_id = np.where(y < dx/10)
    f[bc_id] = 0
    bc_id = np.where(x > 1 - dx/10)
    f[bc_id] = -100*np.sin(np.pi*y[bc_id])**2
    bc_id = np.where(y > 1 - dx/10)
    f[bc_id] = 0

    ### solve linear system
    u_n = sp_la.spsolve(BDF1, f)

    ### write output
    unstructuredGridToVTK(vtk_dir+'/k='+str(k).zfill(4), x, y, z, connectivity = topo, offsets = offset, cell_types = ctype, cellData = None, pointData = {"u" : u_n})
