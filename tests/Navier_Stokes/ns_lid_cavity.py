#! /usr/bin/env python

import time
import os
import sys
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

def assemble_blockwise_force_BDF1():
    size = 2*ndofs_u+ndofs_p#+1
    rhs = np.zeros((size,1))

    f_rhs_x = 1/ph.dt*M11.dot(ux_n)
    f_rhs_y = 1/ph.dt*M22.dot(uy_n)

    #upper boundary
    bc_id = np.where(y_u > 1-delta_x/10)
    f_rhs_x[bc_id,:] = 1.
    f_rhs_y[bc_id,:] = 0

    #lower boundary
    bc_id = np.where(y_u < delta_x/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    #right boundary
    bc_id = np.where(x_u > 1-delta_x/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    #left boundary
    bc_id = np.where(x_u < delta_x/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    rhs[0:ndofs_u,:] = f_rhs_x
    rhs[ndofs_u:2*ndofs_u,:] = f_rhs_y

    return rhs

def assemble_blockwise_matrix_BDF1():
        (K11, K12, K21, K22) = assemble.u_gradv_w_p1(topo_u, x_u, y_u, ux_n, uy_n)
        D11 = 1/ph.dt*M11 + K11 + ph.nu*A11
        D22 = 1/ph.dt*M11 + K22 + ph.nu*A11

        #lower boundary
        bc_id = np.where(y_u < delta_x/10)
        D11 = la_utils.set_diag(D11, bc_id)
        D22 = la_utils.set_diag(D22, bc_id)
        K12 = la_utils.clear_rows(K12, bc_id)
        K21 = la_utils.clear_rows(K21, bc_id)

        #upper boundary
        bc_id = np.where(y_u > 1-delta_x/10)
        D11 = la_utils.set_diag(D11, bc_id)
        D22 = la_utils.set_diag(D22, bc_id)
        K12 = la_utils.clear_rows(K12, bc_id)
        K21 = la_utils.clear_rows(K21, bc_id)

        #left boundary
        bc_id = np.where(x_u < delta_x/10)
        D11 = la_utils.set_diag(D11, bc_id)
        D22 = la_utils.set_diag(D22, bc_id)
        K12 = la_utils.clear_rows(K12, bc_id)
        K21 = la_utils.clear_rows(K21, bc_id)

        #right boundary
        bc_id = np.where(x_u > 1-delta_x/10)
        D11 = la_utils.set_diag(D11, bc_id)
        D22 = la_utils.set_diag(D22, bc_id)
        K12 = la_utils.clear_rows(K12, bc_id)
        K21 = la_utils.clear_rows(K21, bc_id)

        #### assembly of Navier-Stokes system
        mat = sparse.vstack([
            sparse.hstack([D11, K12, -BT1]), #sparse.csr_matrix((ndofs_u,1))]),
            sparse.hstack([K21, D22, -BT2]), #sparse.csr_matrix((ndofs_u,1))]),
            sparse.hstack([-B, sparse.csr_matrix((ndofs_p,ndofs_p))]) #, mean_p.transpose()]),
            #sparse.hstack([sparse.csr_matrix((1,ndofs_u*2)), mean_p, sparse.csr_matrix((1,1))])
        ], "csr")
        return mat

def assemble_blockwise_force_BDF2():
    size = 2*ndofs_u+ndofs_p#+1
    rhs = np.zeros((size,1))

    f_rhs_x = 1/ph.dt*M11.dot(2*ux_n - 0.5*ux_n_old)
    f_rhs_y = 1/ph.dt*M22.dot(2*uy_n - 0.5*uy_n_old)

    #upper boundary
    bc_id = np.where(y_u > 1-delta_x/10)
    f_rhs_x[bc_id,:] = 1.
    f_rhs_y[bc_id,:] = 0

    #lower boundary
    bc_id = np.where(y_u < delta_x/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    #right boundary
    bc_id = np.where(x_u > 1-delta_x/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    #left boundary
    bc_id = np.where(x_u < delta_x/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    rhs[0:ndofs_u,:] = f_rhs_x
    rhs[ndofs_u:2*ndofs_u,:] = f_rhs_y

    return rhs

def assemble_blockwise_matrix_BDF2():
    (K11, K12, K21, K22) = assemble.u_gradv_w_p1(topo_u, x_u, y_u, ux_n, uy_n)
    (L11, L12, L21, L22) = assemble.u_gradv_w_p1(topo_u, x_u, y_u, ux_n_old, uy_n_old)
    D11 = 1.5/ph.dt*M11 + 2*K11 - L11 + ph.nu*A11
    D12 = 2*K12 - L12
    D21 = 2*K21 - L21
    D22 = 1.5/ph.dt*M11 + 2*K22 - L22 + ph.nu*A11

    #lower boundary
    bc_id = np.where(y_u < delta_x/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)
    D12 = la_utils.clear_rows(D12, bc_id)
    D21 = la_utils.clear_rows(D21, bc_id)

    #upper boundary
    bc_id = np.where(y_u > 1-delta_x/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)
    D12 = la_utils.clear_rows(D12, bc_id)
    D21 = la_utils.clear_rows(D21, bc_id)

    #left boundary
    bc_id = np.where(x_u < delta_x/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)
    D12 = la_utils.clear_rows(D12, bc_id)
    D21 = la_utils.clear_rows(D21, bc_id)

    #right boundary
    bc_id = np.where(x_u > 1-delta_x/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)
    D12 = la_utils.clear_rows(D12, bc_id)
    D21 = la_utils.clear_rows(D21, bc_id)

    #### assembly of Navier-Stokes system
    mat = sparse.vstack([
        sparse.hstack([D11, D12, -BT1]), #sparse.csr_matrix((ndofs_u,1))]),
        sparse.hstack([D21, D22, -BT2]), #sparse.csr_matrix((ndofs_u,1))]),
        sparse.hstack([-B, sparse.csr_matrix((ndofs_p,ndofs_p))]) #, mean_p.transpose()]),
        #sparse.hstack([sparse.csr_matrix((1,ndofs_u*2)), mean_p, sparse.csr_matrix((1,1))])
    ], "csr")
    return mat

# def apply_bc(A):
#     bc_id = np.where(y_u < delta_x/10)
#     A = la_utils.set_diag(A, bc_id)
#     bc_id = bc_id + np.ones(len(bc_id[0])).astype(int)*ndofs_u
#     A = la_utils.set_diag(A, bc_id)
#     bc_id = np.where(x_u < delta_x/10)
#     A = la_utils.set_diag(A, bc_id)
#     bc_id = bc_id + np.ones(len(bc_id[0])).astype(int)*ndofs_u
#     A = la_utils.set_diag(A, bc_id)
#     bc_id = np.where(y_u > 1 - delta_x/10)
#     A = la_utils.set_diag(A, bc_id)
#     bc_id = bc_id + np.ones(len(bc_id[0])).astype(int)*ndofs_u
#     A = la_utils.set_diag(A, bc_id)
#     bc_id = np.where(x_u > 1 - delta_x/10)
#     A = la_utils.set_diag(A, bc_id)
#     bc_id = bc_id + np.ones(len(bc_id[0])).astype(int)*ndofs_u
#     A = la_utils.set_diag(A, bc_id)
#
#     return A

def write_output():
    filename = results_dir +'cn_time_'+str(cn_time).zfill(ph.time_index_digits)
    f = file(filename,"wb")
    np.save(f,u_n)
    np.save(f,p_n)
    f.close()
    print '--------------------------------------'
    print 'results saved to:'
    print filename
    print '--------------------------------------'
    return

def l2_norm(M, g):
    #print g
    l2_g = M.dot(g)
    l2_g = np.dot(l2_g.transpose(),g)
    l2_g = mth.sqrt(l2_g)
    return l2_g

np.set_printoptions(precision=4)
np.set_printoptions(suppress=True)

if len(sys.argv) > 1:
    ph = ParametersHandler(sys.argv[1])
else:
    ph = ParametersHandler('simulation_parameters.json')
ph.simulation_info()

nx_p = ph.n_delta_x
delta_x = 1./nx_p
ny_p = nx_p
delta_y = 1./ny_p
(topo_p,x_p,y_p,
    topo_u,x_u,y_u,
    c2f) = lin_t3.mesh_t3_iso_t6(nx_p, ny_p,delta_x,delta_y)

(topo_p,x_p,y_p) = lin_t3.mesh_t3_t0(nx_p,ny_p,delta_x,delta_y)
tp = topo_p
#(topo_p, x_p, y_p, topo_u, x_u, y_u, c2f) = lin_t3.load_t3_iso_t6_file('../mesh_collection/step.msh', 'mesh_collection/step_refined.msh')
#(topo_p, x_p, y_p, topo_u, x_u, y_u, c2f) = lin_t3.load_t3_iso_t6_file('../mesh_collection/cavity_8.msh', 'mesh_collection/cavity_16.msh')

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
ux_n_old = np.zeros((ndofs_u,1))
uy_n_old = np.zeros((ndofs_u,1))

M11 = assemble.u_v_p1(topo_u,x_u,y_u)
A11 = assemble.gradu_gradv_p1(topo_u,x_u,y_u)

# D11 = 1/ph.dt*M11 + ph.nu*A11
# D22 = D11
M22 = M11

(BT1,BT2) = assemble.divu_p_p1_iso_p2_p1p0(topo_p,x_p,y_p,
           topo_u,x_u,y_u,c2f)

BT = sparse.vstack([BT1,BT2])
B = BT.transpose()

#left boundary
bc_id = np.where(x_u < delta_x/10)
#D11 = la_utils.set_diag(D11, bc_id)
#D22 = la_utils.set_diag(D22, bc_id)
BT1 = la_utils.clear_rows(BT1,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)

#right boundary
bc_id = np.where(x_u > 1-delta_x/10)
#D11 = la_utils.set_diag(D11, bc_id)
#D22 = la_utils.set_diag(D22, bc_id)
BT1 = la_utils.clear_rows(BT1,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)

#lower boundary
bc_id = np.where(y_u < delta_x/10)
#D11 = la_utils.set_diag(D11, bc_id)
#D22 = la_utils.set_diag(D22, bc_id)
BT1 = la_utils.clear_rows(BT1,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)

#upper boundary
bc_id = np.where(y_u > 1-delta_x/10)
#D11 = la_utils.set_diag(D11, bc_id)
#D22 = la_utils.set_diag(D22, bc_id)
BT1 = la_utils.clear_rows(BT1,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)


M = sparse.vstack([
     sparse.hstack( [M11, sparse.csr_matrix((ndofs_u,ndofs_u))] ),
     sparse.hstack( [sparse.csr_matrix((ndofs_u,ndofs_u)), M11] )
     ])

# mean_p = np.zeros((1,ndofs_p))
# x_l = x_p[topo_p[0,0:3]]
# y_l = y_p[topo_p[0,0:3]]
# eval_p = np.zeros((0,2))
# (phi_dx,phi_dy,phi,omega) = shp.tri_p1(x_l,y_l,eval_p)
#
# for row in topo_p:
#     mean_p[0,row] += omega * np.array([1./3.,1./3.,1./3.,1])
#
# mat = sparse.vstack([
#     sparse.hstack([D11, sparse.csr_matrix((ndofs_u, ndofs_u)), -BT1]), #sparse.csr_matrix((ndofs_u,1))]),
#     sparse.hstack([sparse.csr_matrix((ndofs_u, ndofs_u)), D22, -BT2]), #sparse.csr_matrix((ndofs_u,1))]),
#     sparse.hstack([-B, sparse.csr_matrix((ndofs_p,ndofs_p))])#,mean_p.transpose()]),
#     #sparse.hstack([sparse.csr_matrix((1,ndofs_u*2)), mean_p, sparse.csr_matrix((1,1))])
# ], "csr")

#plt.spy(mat)
#plt.show()

#### start time steppig procedure
for cn_time in range(0,len(ph.stampa)):
    step_t0 = time.time()


    mat = sparse.csr_matrix((2*ndofs_u + ndofs_p, 2*ndofs_u + ndofs_p))
    force = sparse.csr_matrix((2*ndofs_u + ndofs_p, 1))
    if ph.time_integration == 'BDF1':
        mat = assemble_blockwise_matrix_BDF1()
        force = assemble_blockwise_force_BDF1()
    elif ph.time_integration == 'BDF2':
        mat = assemble_blockwise_matrix_BDF2()
        force = assemble_blockwise_force_BDF2()

    ux_n_old = ux_n
    uy_n_old = uy_n

    sol_t0 = time.time()
    sol = sp_la.spsolve(mat,force)
    sol_t1 = time.time()

    u_n = sol[0:2*ndofs_u]
    p_n = sol[2*ndofs_u:2*ndofs_u+ndofs_p]

    ux_n = np.reshape(u_n[0      :  ndofs_u], ux_n.shape)
    uy_n = np.reshape(u_n[ndofs_u:2*ndofs_u], uy_n.shape)

    if ph.stampa[cn_time] == True:
        write_output()
    step_t1 = time.time()

    print '--------------------------------------'
    print 'cn_time   = ' + str(cn_time)
    print 't         = ' + str(cn_time*ph.dt)
    print 'l2 norm u = ' + str(l2_norm(M, u_n))
    print 'step time = ' + str((step_t1-step_t0))
    print 'sol  time = ' + str((sol_t1-sol_t0))
    print '--------------------------------------'
