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

def write_output():
    filename = results_dir +'cn_time_'+str(cn_time).zfill(ph.time_index_digits)
    f = file(filename,"wb")
    np.save(f,u_n)
    np.save(f,p_n)
    f.close()
    # print '--------------------------------------'
    # print 'results saved to:'
    # print filename
    # print '--------------------------------------'
    return

def write_analytical():
    filename = results_dir +'analytical'
    f = file(filename,"wb")
    np.save(f,analytical)
    np.save(f,p_n)
    f.close()
    # print '--------------------------------------'
    # print 'results saved to:'
    # print filename
    # print '--------------------------------------'
    return

def l2_norm(M, g):
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

#(topo_p,x_p,y_p) = lin_t3.mesh_t3_t0(nx_p,ny_p,delta_x,delta_y)

if sum(ph.stampa) !=0:
    results_dir = ph.results_directory+'/'+ph.sim_prefix+'/binary_data/'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    ph.dump_to_json(ph.results_directory+'/'+ph.sim_prefix+'/simulation_parameters.json')
    write_mesh()

ndofs_u = max(x_u.shape)
ndofs_p = max(x_p.shape)# + topo_p.shape[0]

print(ndofs_u)
print(ndofs_p)

M11 = assemble.u_v_p1(topo_u,x_u,y_u)
K11 = assemble.gradu_gradv_p1(topo_u,x_u,y_u)

(BT1,BT2) = assemble.divu_p_p1_iso_p2_p1(topo_p,x_p,y_p,
           topo_u,x_u,y_u,c2f)

mean_p = np.zeros((1,ndofs_p))

for row in topo_p:
    x_l = x_p[row]
    y_l = y_p[row]
    eval_p = np.zeros((0,2))
    (phi_dx,phi_dy,phi,omega) = shp.tri_p1(x_l,y_l,eval_p)
    mean_p[0,row] += omega * np.array([1./3.,1./3.,1./3.])

BT = sparse.vstack([BT1,BT2])
B = BT.transpose()

m_BDF1 =  1/ph.dt*M11   + K11
m_Theta = 1/ph.dt*M11   + 0.5*K11
m_BDF2 =  1.5/ph.dt*M11 + K11

#left bnd
bc_id = np.where(x_u < delta_x/10)
m_BDF1 = la_utils.set_diag(m_BDF1,bc_id)
m_BDF2 = la_utils.set_diag(m_BDF2,bc_id)
m_Theta = la_utils.set_diag(m_Theta,bc_id)
BT1 = la_utils.clear_rows(BT1, bc_id)
BT2 = la_utils.clear_rows(BT2, bc_id)

#right bnd
bc_id = np.where(x_u > 1-delta_x/10)
m_BDF1 = la_utils.set_diag(m_BDF1,bc_id)
m_BDF2 = la_utils.set_diag(m_BDF2,bc_id)
m_Theta = la_utils.set_diag(m_Theta,bc_id)
BT1 = la_utils.clear_rows(BT1, bc_id)
BT2 = la_utils.clear_rows(BT2, bc_id)

#lower bnd
bc_id = np.where(y_u < delta_x/10)
m_BDF1 = la_utils.set_diag(m_BDF1,bc_id)
m_BDF2 = la_utils.set_diag(m_BDF2,bc_id)
m_Theta = la_utils.set_diag(m_Theta,bc_id)
BT1 = la_utils.clear_rows(BT1, bc_id)
BT2 = la_utils.clear_rows(BT2, bc_id)

#upper bnd
bc_id = np.where(y_u > 1-delta_x/10)
m_BDF1 = la_utils.set_diag(m_BDF1,bc_id)
m_BDF2 = la_utils.set_diag(m_BDF2,bc_id)
m_Theta = la_utils.set_diag(m_Theta,bc_id)
BT1 = la_utils.clear_rows(BT1, bc_id)
BT2 = la_utils.clear_rows(BT2, bc_id)

BT = sparse.vstack([BT1, BT2])



K = sparse.vstack([
    sparse.hstack([K11, sparse.csr_matrix((ndofs_u, ndofs_u))]),
    sparse.hstack([sparse.csr_matrix((ndofs_u, ndofs_u)), K11])
], "csr")
M = sparse.vstack([
    sparse.hstack([M11, sparse.csr_matrix((ndofs_u, ndofs_u))]),
    sparse.hstack([sparse.csr_matrix((ndofs_u, ndofs_u)), M11])
], "csr")
mat_BDF1 = sparse.vstack([
    sparse.hstack([m_BDF1, sparse.csr_matrix((ndofs_u, ndofs_u)), -BT1, sparse.csr_matrix((ndofs_u, 1))]),
    sparse.hstack([sparse.csr_matrix((ndofs_u, ndofs_u)), m_BDF1, -BT2, sparse.csr_matrix((ndofs_u, 1))]),
    sparse.hstack([-B, sparse.csr_matrix((ndofs_p, ndofs_p)), mean_p.transpose()]),
    sparse.hstack([sparse.csr_matrix((1, 2*ndofs_u)), mean_p, sparse.csr_matrix((1,1))])
], "csr")
mat_BDF2 = sparse.vstack([
    sparse.hstack([m_BDF2, sparse.csr_matrix((ndofs_u, ndofs_u)), -BT1, sparse.csr_matrix((ndofs_u, 1))]),
    sparse.hstack([sparse.csr_matrix((ndofs_u, ndofs_u)), m_BDF2, -BT2, sparse.csr_matrix((ndofs_u, 1))]),
    sparse.hstack([-B, sparse.csr_matrix((ndofs_p, ndofs_p)), mean_p.transpose()]),
    sparse.hstack([sparse.csr_matrix((1, 2*ndofs_u)), mean_p, sparse.csr_matrix((1,1))])
], "csr")
mat_Theta = sparse.vstack([
    sparse.hstack([m_Theta, sparse.csr_matrix((ndofs_u, ndofs_u)), -0.5*BT1, sparse.csr_matrix((ndofs_u, 1))]),
    sparse.hstack([sparse.csr_matrix((ndofs_u, ndofs_u)), m_Theta, -0.5*BT2, sparse.csr_matrix((ndofs_u, 1))]),
    sparse.hstack([-B, sparse.csr_matrix((ndofs_p, ndofs_p)), mean_p.transpose()]),
    sparse.hstack([sparse.csr_matrix((1, 2*ndofs_u)), mean_p, sparse.csr_matrix((1,1))])
], "csr")

un_x = np.zeros((ndofs_u, 1))
un_y = np.zeros((ndofs_u, 1))
un_x_old = np.zeros((ndofs_u, 1))
un_y_old = np.zeros((ndofs_u, 1))
u_n = np.zeros((2*ndofs_u, 1))
p_n = np.zeros((ndofs_p, 1))


f_x = (2-12*x_u+12*x_u**2)*(2*y_u-6*y_u**2+4*y_u**3) + x_u**2*(1-x_u)**2*(-12+24*y_u) - 1
f_y = -(-12+24*x_u)*y_u**2*(1-y_u)**2-(2*x_u-6*x_u**2+4*x_u**3)*(2-12*y_u+12*y_u**2)
f_x = np.reshape(f_x, (ndofs_u, 1))
f_y = np.reshape(f_y, (ndofs_u, 1))
analytical_x = np.reshape(-x_u**2*(1-x_u)**2*(2*y_u-6*y_u**2+4*y_u**3), (ndofs_u, 1))
analytical_y = np.reshape((2*x_u-6*x_u**2+4*x_u**3)*y_u**2*(1-y_u)**2, (ndofs_u, 1))
analytical = np.reshape(np.append(analytical_x, analytical_y), (2*ndofs_u, 1))

err = np.zeros(len(ph.stampa))
for cn_time in range(0,len(ph.stampa)):
    if ph.time_integration == 'BDF1':
        rhs_u_x = 1/ph.dt*M11.dot(un_x) \
            + np.reshape(M11.dot(mth.cos(-(cn_time+1)*ph.dt)*analytical_x + (1-mth.sin(-(cn_time+1)*ph.dt))*f_x), (ndofs_u, 1))
        rhs_u_y = 1/ph.dt*M11.dot(un_y) \
            + np.reshape(M11.dot(mth.cos(-(cn_time+1)*ph.dt)*analytical_y + (1-mth.sin(-(cn_time+1)*ph.dt))*f_y), (ndofs_u, 1))
    elif ph.time_integration == 'BDF2':
        rhs_u_x = 1/ph.dt*M11.dot(2*un_x-0.5*un_x_old) \
            + np.reshape(M11.dot(mth.cos(-(cn_time+1)*ph.dt)*analytical_x + (1-mth.sin(-(cn_time+1)*ph.dt))*f_x), (ndofs_u, 1))
        rhs_u_y = 1/ph.dt*M11.dot(2*un_y-0.5*un_y_old) \
            + np.reshape(M11.dot(mth.cos(-(cn_time+1)*ph.dt)*analytical_y + (1-mth.sin(-(cn_time+1)*ph.dt))*f_y), (ndofs_u, 1))
    else:
        rhs_u_x = 1/ph.dt*M11.dot(un_x) - 0.5*K11.dot(un_x) + 0.5*BT1.dot(p_n) \
            - 0.5*np.reshape(M11.dot(mth.cos(-cn_time*ph.dt)*analytical_x + (1-mth.sin(-cn_time*ph.dt))*f_x), (ndofs_u, 1)) \
            - 0.5*np.reshape(M11.dot(mth.cos(-(cn_time+1)*ph.dt)*analytical_x + (1-mth.sin(-(cn_time+1)*ph.dt))*f_x), (ndofs_u, 1))
        rhs_u_y = 1/ph.dt*M11.dot(un_y) - 0.5*K11.dot(un_x) + 0.5*BT2.dot(p_n) \
            - 0.5*np.reshape(M11.dot(mth.cos(-cn_time*ph.dt)*analytical_y + (1-mth.sin(-cn_time*ph.dt))*f_y), (ndofs_u, 1)) \
            - 0.5*np.reshape(M11.dot(mth.cos(-(cn_time+1)*ph.dt)*analytical_y + (1-mth.sin(-(cn_time+1)*ph.dt))*f_y), (ndofs_u, 1))
    rhs_p = np.zeros((ndofs_p, 1))

    #left bnd
    bc_id = np.where(x_u < delta_x/10)
    rhs_u_x[bc_id,:] = 0
    rhs_u_y[bc_id,:] = 0

    #right bnd
    bc_id = np.where(x_u > 1-delta_x/10)
    rhs_u_x[bc_id,:] = 0
    rhs_u_y[bc_id,:] = 0

    #lower bnd
    bc_id = np.where(y_u < delta_x/10)
    rhs_u_x[bc_id,:] = 0
    rhs_u_y[bc_id,:] = 0

    #upper bnd
    bc_id = np.where(y_u > 1-delta_x/10)
    rhs_u_x[bc_id,:] = 0
    rhs_u_y[bc_id,:] = 0

    rhs = sparse.vstack([
        rhs_u_x,
        rhs_u_y,
        rhs_p,
        np.zeros((1,1))
    ])

    if ph.time_integration == 'BDF1':
        sol = sp_la.spsolve(mat_BDF1, rhs)
    elif ph.time_integration == 'Theta' or (ph.time_integration == 'BDF2' and cn_time == 0):
        sol = sp_la.spsolve(mat_Theta, rhs)
    elif ph.time_integration == 'BDF2':
        sol = sp_la.spsolve(mat_BDF2, rhs)

    un_x_old = un_x
    un_y_old = un_y
    un_x = np.reshape(sol[0        :ndofs_u],           (ndofs_u  , 1))
    un_y = np.reshape(sol[ndofs_u  :2*ndofs_u],         (ndofs_u  , 1))
    u_n =  np.reshape(sol[0        :2*ndofs_u],         (2*ndofs_u, 1))
    p_n =  np.reshape(sol[2*ndofs_u:2*ndofs_u+ndofs_p], (ndofs_p  , 1))

    write_output()

    err[cn_time] = l2_norm(M, u_n - (1-mth.sin(-cn_time*ph.dt))*analytical)
    print(err[cn_time])
print(err[len(err)-1])
print(np.mean(err))
