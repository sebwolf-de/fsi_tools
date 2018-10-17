#! /usr/bin/env python

import time
import os
import numpy as np
from scipy import sparse
import scipy.sparse.linalg as sp_la
import matplotlib.pyplot as plt
import math as mth
import json

# nicola modules
import lin_tri_mesh as lin_t3
import basis_func as shp
import assemble
import la_utils
import viewers
import geom_utils as geom
from shapely.geometry import Polygon

from preconditioner import BlockPreconditioner
from parameters_handler import ParametersHandler

def write_mesh():
    filename = results_dir+'mesh'#'./mesh/'+sim_prefix
    f = file(filename,"wb")
    np.save(f,topo_p)
    np.save(f,x_p)
    np.save(f,y_p)
    np.save(f,topo_u)
    np.save(f,x_u)
    np.save(f,y_u)
    np.save(f,c2f)
    f.close()
    return

def assemble_blockwise_force():
    size = 2*ndofs_u+ndofs_p+1
    rhs = np.zeros((size))

    f_rhs = 1/ph.dt*M.dot(u_n)
    f_rhs_x = np.reshape(f_rhs[0:ndofs_u], (ndofs_u, 1))
    f_rhs_y = np.reshape(f_rhs[ndofs_u:2*ndofs_u], (ndofs_u, 1))

    bc_id = np.where(y_u > 1-delta_x/10)
    f_rhs_x[bc_id,:] = 1.
    f_rhs_y[bc_id,:] = 0

    bc_id = np.where(y_u < delta_x/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    bc_id = np.where(x_u > 1-delta_x/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    bc_id = np.where(x_u < delta_x/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    f_rhs_x = np.reshape(f_rhs_x,(ndofs_u))
    f_rhs_y = np.reshape(f_rhs_y,(ndofs_u))

    rhs[0:ndofs_u] = f_rhs_x
    rhs[ndofs_u:2*ndofs_u] = f_rhs_y

    return rhs

def assemble_blockwise_matrix():
    (K11, K12, K21, K22) = assemble.u_gradv_w_p1(topo_u, x_u, y_u, ux_n, uy_n)
    D11 = ph.rho_fluid/ph.dt*M11 + ph.nu*A11 + ph.rho_fluid*K11
    D22 = ph.rho_fluid/ph.dt*M11 + ph.nu*A11 + ph.rho_fluid*K11

    bc_id = np.where(y_u < delta_x/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D11, bc_id)
    K12 = la_utils.clear_rows(K12, bc_id)
    K21 = la_utils.clear_rows(K21, bc_id)

    bc_id = np.where(y_u > 1-delta_x/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D11, bc_id)
    K12 = la_utils.clear_rows(K12, bc_id)
    K21 = la_utils.clear_rows(K21, bc_id)

    bc_id = np.where(x_u < delta_x/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D11, bc_id)
    K12 = la_utils.clear_rows(K12, bc_id)
    K21 = la_utils.clear_rows(K21, bc_id)

    bc_id = np.where(x_u > 1-delta_x/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D11, bc_id)
    K12 = la_utils.clear_rows(K12, bc_id)
    K21 = la_utils.clear_rows(K21, bc_id)

    #### assembly of Navier-Stokes system
    mat = sparse.vstack([
        sparse.hstack([D11, ph.rho_fluid*K12, -BT1, sparse.csr_matrix((ndofs_u,1))]),
        sparse.hstack([ph.rho_fluid*K21, D22, -BT2, sparse.csr_matrix((ndofs_u,1))]),
        sparse.hstack([-BT1.transpose(), -BT2.transpose(), sparse.csr_matrix((ndofs_p,ndofs_p)), mean_p.transpose()]),
        sparse.hstack([sparse.csr_matrix((1,ndofs_u*2)), mean_p, sparse.csr_matrix((1,1))])
    ], "csr")

    #plt.spy(mat)
    #plt.show()
    #### assembly of Stokes system
    #mat = sparse.vstack([
    #    sparse.hstack([M11/ph.dt + A11, sparse.csr_matrix((ndofs_u, ndofs_u)), -BT1, sparse.csr_matrix((ndofs_u,1))]),
    #    sparse.hstack([sparse.csr_matrix((ndofs_u, ndofs_u)), M11/ph.dt + A11, -BT2, sparse.csr_matrix((ndofs_u,1))]),
    #    sparse.hstack([-BT1.transpose(), -BT2.transpose(), sparse.csr_matrix((ndofs_p,ndofs_p)), mean_p.transpose()]),
    #    sparse.hstack([sparse.csr_matrix((1,ndofs_u*2)), mean_p, sparse.csr_matrix((1,1))])
    #], "csr")
    #mat = apply_bc(mat)

    return mat

def apply_bc(A):
    bc_id = np.where(y_u < delta_x/10)
    A = la_utils.set_diag(A, bc_id)
    bc_id = bc_id + np.ones(len(bc_id[0])).astype(int)*ndofs_u
    A = la_utils.set_diag(A, bc_id)
    bc_id = np.where(x_u < delta_x/10)
    A = la_utils.set_diag(A, bc_id)
    bc_id = bc_id + np.ones(len(bc_id[0])).astype(int)*ndofs_u
    A = la_utils.set_diag(A, bc_id)
    bc_id = np.where(y_u > 1 - delta_x/10)
    A = la_utils.set_diag(A, bc_id)
    bc_id = bc_id + np.ones(len(bc_id[0])).astype(int)*ndofs_u
    A = la_utils.set_diag(A, bc_id)
    bc_id = np.where(x_u > 1 - delta_x/10)
    A = la_utils.set_diag(A, bc_id)
    bc_id = bc_id + np.ones(len(bc_id[0])).astype(int)*ndofs_u
    A = la_utils.set_diag(A, bc_id)

    return A

def write_output():
    filename = results_dir +'cn_time_'+str(cn_time).zfill(ph.time_index_digits)
    f = file(filename,"wb")
    np.save(f,u_n)
    np.save(f,p)
    f.close()
    print '--------------------------------------'
    print 'results saved to:'
    print filename
    print '--------------------------------------'
    return

def l2_norm(g):
    l2_g = np.dot(g.transpose(),g)
    l2_g = mth.sqrt(l2_g)
    return l2_g

np.set_printoptions(precision=4)
np.set_printoptions(suppress=True)

ph = ParametersHandler('simulation_parameters_ns.json')
ph.simulation_info()

nx_p = ph.n_delta_x
delta_x = 1./nx_p
ny_p = nx_p
delta_y = 1./ny_p
#(topo_p,x_p,y_p,
# topo_u,x_u,y_u,
# c2f) = lin_t3.mesh_t3_iso_t6(nx_p, ny_p,delta_x,delta_y)
#(topo_p,x_p,y_p) = lin_t3.mesh_t3_t0(nx_p,ny_p,delta_x,delta_y)

(topo_p, x_p, y_p, topo_u, x_u, y_u, c2f) = lin_t3.load_t3_iso_t6_file('mesh_collection/step.msh', 'mesh_collection/step_refined.msh')

print topo_p
if sum(ph.stampa) !=0:
    results_dir = ph.results_directory+'/'+ph.sim_prefix+'/binary_data/'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    ph.dump_to_json(ph.results_directory+'/'+ph.sim_prefix+'/simulation_parameters.json')
    write_mesh()

ndofs_u = max(x_u.shape)
ndofs_p = max(x_p.shape) + topo_p.shape[0]

u_n = np.zeros((2*ndofs_u,1))
ux_n = np.zeros((ndofs_u,1))
uy_n = np.zeros((ndofs_u,1))

M11 = assemble.u_v_p1(topo_u,x_u,y_u)
A11 = assemble.gradu_gradv_p1(topo_u,x_u,y_u)

(BT1,BT2) = assemble.divu_p_p1_iso_p2_p1p0(topo_p,x_p,y_p,
           topo_u,x_u,y_u,c2f)

bc_id = np.where(x_u < delta_x/10)
BT1 = la_utils.clear_rows(BT1,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)

bc_id = np.where(x_u > 1-delta_x/10)
BT1 = la_utils.clear_rows(BT1,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)

bc_id = np.where(y_u < delta_x/10)
BT1 = la_utils.clear_rows(BT1,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)

bc_id = np.where(y_u > 1-delta_x/10)
BT1 = la_utils.clear_rows(BT1,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)

BT = sparse.vstack([BT1,BT2])
B = BT.transpose()

M = sparse.vstack([
    sparse.hstack( [M11, sparse.csr_matrix((ndofs_u,ndofs_u))] ),
    sparse.hstack( [sparse.csr_matrix((ndofs_u,ndofs_u)), M11] )
    ])

mean_p = np.zeros((1,ndofs_p))
x_l = x_p[topo_p[0,0:3]]
y_l = y_p[topo_p[0,0:3]]
eval_p = np.zeros((0,2))
(phi_dx,phi_dy,phi,omega) = shp.tri_p1(x_l,y_l,eval_p)

for row in topo_p:
    mean_p[0,row] += omega * np.array([1./3.,1./3.,1./3.,1])

#### start time steppig procedure
for cn_time in range(0,len(ph.stampa)):
    step_t0 = time.time()

    mat = assemble_blockwise_matrix()
    force = assemble_blockwise_force()


    sol_t0 = time.time()
    sol = sp_la.spsolve(mat,force)
    sol_t1 = time.time()

    u = sol[0:2*ndofs_u]
    p = sol[2*ndofs_u:2*ndofs_u+ndofs_p]

    u_n = u
    ux_n1 = u[0      :  ndofs_u]
    uy_n1 = u[ndofs_u:2*ndofs_u]

    ux_n = np.reshape(ux_n1, ux_n.shape)
    uy_n = np.reshape(uy_n1, uy_n.shape)

    if ph.stampa[cn_time] == True:
        write_output()
    step_t1 = time.time()

    print '--------------------------------------'
    print 'cn_time   = ' + str(cn_time)
    print 't         = ' + str(cn_time*ph.dt)
    print 'l2 norm u = ' + str(l2_norm(u_n))
    print 'l2 norm p = ' + str(l2_norm(p))
    print 'step time = ' + str((step_t1-step_t0))
    print 'sol  time = ' + str((sol_t1-sol_t0))
    print '--------------------------------------'
