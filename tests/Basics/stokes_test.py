#! /bin/env python

import gc
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse.linalg as sp_la
from scipy import sparse
import time
import sys

sys.path.append('../../modules')
import assemble
import basis_func as shp
import la_utils
import lin_tri_mesh as lin_t3
import quadrature
import viewers

def analytical_u(t):
    analytical_x = 4*np.sin(t) * x_u**2 * (1-x_u)**2 * (2*y_u - 6*y_u**2 + 4*y_u**3)
    analytical_y = -4*np.sin(t) * (2*x_u - 6*x_u**2 + 4*x_u**3) * y_u**2 * (1-y_u)**2
    analytical = np.reshape(np.append(analytical_x, analytical_y), (2*ndofs_u, 1))
    return analytical

def analytical_u_1(t, x, y):
    analytical_x = 4*np.sin(t) * x**2 * (1-x)**2 * (2*y - 6*y**2 + 4*y**3)
    return analytical_x

def analytical_u_2(t, x, y):
    analytical_y = -4*np.sin(t) * (2*x - 6*x**2 + 4*x**3) * y**2 * (1-y)**2
    return analytical_y

def analytical_p(t):
    return np.reshape((x_p-0.5), (ndofs_p, 1))

def analytical(t):
    return sparse.vstack([analytical_u(t), analytical_p(t)])


def f(t):
    ## time derivative
    f_x =  4*np.cos(t) * x_u**2 * (1-x_u)**2 * (2*y_u - 6*y_u**2 + 4*y_u**3)
    f_y = -4*np.cos(t) * (2*x_u - 6*x_u**2 + 4*x_u**3) * y_u**2 * (1-y_u)**2
    ## convection term
    # f_x +=  16*np.sin(t)**2 * x_u**2 * (1-x_u)**2 * (2*y_u - 6*y_u**2 + 4*y_u**3)**2 * (2*x_u - 6*x_u**2 + 4*x_u**3)
    # f_y += -16*np.sin(t)**2 * x_u**2 * (1-x_u)**2 * (2*y_u - 6*y_u**2 + 4*y_u**3) * (2 - 12*x_u + 12*x_u**2) * y_u**2 * (1-y_u)**2
    # f_x += -16*np.sin(t)**2 * (2*x_u - 6*x_u**2 + 4*x_u**3) * y_u**2 * (1-y_u)**2 * x_u**2 * (1-x_u)**2 * (2 - 12*y_u + 12*y_u**2)
    # f_y +=  16*np.sin(t)**2 * (2*x_u - 6*x_u**2 + 4*x_u**3)**2 * y_u**2 * (1-y_u)**2 * (2*y_u - 6*y_u**2 + 4*y_u**3)
    ## diffusion term
    f_x += -4*np.sin(t) * ((2 - 12*x_u + 12*x_u**2) * (2*y_u - 6*y_u**2 + 4*y_u**3) + x_u**2 * (1-x_u)**2 * (-12 + 24*y_u))
    f_y += -4*np.sin(t) * (-(-12 + 24*x_u) * y_u**2 * (1-y_u)**2 - (2*x_u - 6*x_u**2 + 4*x_u**3) * (2 - 12*y_u + 12*y_u**2))
    ## pressure gradient
    f_x +=  -1
    f_y +=  0
    f_stacked = np.reshape(np.append(f_x, f_y), (2*ndofs_u, 1))
    return f_stacked

def analytical_f_1(t, x, y):
    f_x =  4*np.cos(t) * x**2 * (1-x)**2 * (2*y - 6*y**2 + 4*y**3)
    ## diffusion term
    f_x += -4*np.sin(t) * ((2 - 12*x + 12*x**2) * (2*y - 6*y**2 + 4*y**3) + x**2 * (1-x)**2 * (-12 + 24*y))
    ## pressure gradient
    f_x +=  -1
    return f_x

def analytical_f_2(t, x, y):
    ## time derivative
    f_y = -4*np.cos(t) * (2*x - 6*x**2 + 4*x**3) * y**2 * (1-y)**2
    ## diffusion term
    f_y += -4*np.sin(t) * (-(-12 + 24*x) * y**2 * (1-y)**2 - (2*x - 6*x**2 + 4*x**3) * (2 - 12*y + 12*y**2))
    ## pressure gradient
    f_y +=  0
    return f_y


def assemble_blockwise_force_BDF1(t, dt):
    size = 2*ndofs_u+ndofs_p+1
    rhs = np.zeros((size,1))

    g_1 = lambda x,y: analytical_f_1(t, x, y)
    g_2 = lambda x,y: analytical_f_2(t, x, y)
    g_1_rhs = quadrature.fem_rhs(g_1, x_u, y_u, topo_u)
    g_2_rhs = quadrature.fem_rhs(g_2, x_u, y_u, topo_u)
    f_rhs_x = 1./dt*M.dot(u_BDF1[0:ndofs_u]) + np.reshape(g_1_rhs, (ndofs_u, 1))
    f_rhs_y = 1./dt*M.dot(u_BDF1[ndofs_u:2*ndofs_u]) + np.reshape(g_2_rhs, (ndofs_u, 1))

    #upper boundary
    bc_id = np.where(y_u > 1-dx/10)
    f_rhs_x[bc_id,:] = 0.
    f_rhs_y[bc_id,:] = 0.

    #lower boundary
    bc_id = np.where(y_u < dx/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    #right boundary
    bc_id = np.where(x_u > 1-dx/10)
    f_rhs_x[bc_id,:] = 0.
    f_rhs_y[bc_id,:] = 0.

    #left boundary
    bc_id = np.where(x_u < dx/10)
    f_rhs_x[bc_id,:] = 0.
    f_rhs_y[bc_id,:] = 0.

    rhs[0:ndofs_u,:] = f_rhs_x
    rhs[ndofs_u:2*ndofs_u,:] = f_rhs_y

    return np.reshape(rhs, (size))

def assemble_blockwise_matrix_BDF1(dt):
    D11 = 1./dt*M + K
    D22 = 1./dt*M + K
    S12 = sparse.csr_matrix((ndofs_u, ndofs_u))
    S21 = sparse.csr_matrix((ndofs_u, ndofs_u))

    #lower boundary
    bc_id = np.where(y_u < dx/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)
    S12 = la_utils.clear_rows(S12, bc_id)
    S21 = la_utils.clear_rows(S21, bc_id)

    #upper boundary
    bc_id = np.where(y_u > 1-dx/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)
    S12 = la_utils.clear_rows(S12, bc_id)
    S21 = la_utils.clear_rows(S21, bc_id)

    #left boundary
    bc_id = np.where(x_u < dx/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)
    S12 = la_utils.clear_rows(S12, bc_id)
    S21 = la_utils.clear_rows(S21, bc_id)

    #right boundary
    bc_id = np.where(x_u > 1-dx/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)
    S12 = la_utils.clear_rows(S12, bc_id)
    S21 = la_utils.clear_rows(S21, bc_id)


    #### assembly of Navier-Stokes system
    mat = sparse.vstack([
        sparse.hstack([D11, S12, -BT1, sparse.csr_matrix((ndofs_u, 1))]),
        sparse.hstack([S21, D22, -BT2, sparse.csr_matrix((ndofs_u, 1))]),
        sparse.hstack([-B, sparse.csr_matrix((ndofs_p,ndofs_p)), mean_p.transpose()]),
        sparse.hstack([sparse.csr_matrix((1, 2*ndofs_u)), mean_p, sparse.csr_matrix((1,1))])
    ], "csr")
    return mat

def assemble_blockwise_force_BDF2(t, dt):
    size = 2*ndofs_u+ndofs_p+1
    rhs = np.zeros((size,1))

    g_1 = lambda x,y: analytical_f_1(t, x, y)
    g_2 = lambda x,y: analytical_f_2(t, x, y)
    g_1_rhs = quadrature.fem_rhs(g_1, x_u, y_u, topo_u)
    g_2_rhs = quadrature.fem_rhs(g_2, x_u, y_u, topo_u)
    f_rhs_x = 1./dt*M.dot(2*u_BDF2[0:ndofs_u] - 0.5*u_BDF2_old[0:ndofs_u]) + np.reshape(g_1_rhs, (ndofs_u, 1))
    f_rhs_y = 1./dt*M.dot(2*u_BDF2[ndofs_u:2*ndofs_u] - 0.5*u_BDF2_old[ndofs_u:2*ndofs_u]) + np.reshape(g_2_rhs, (ndofs_u, 1))

    #upper boundary
    bc_id = np.where(y_u > 1-dx/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    #lower boundary
    bc_id = np.where(y_u < dx/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    #right boundary
    bc_id = np.where(x_u > 1-dx/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    #left boundary
    bc_id = np.where(x_u < dx/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    rhs[0:ndofs_u,:] = f_rhs_x
    rhs[ndofs_u:2*ndofs_u,:] = f_rhs_y

    return np.reshape(rhs, (size))

def assemble_blockwise_matrix_BDF2(dt):
    D11 = 1.5/dt*M + K
    D22 = 1.5/dt*M + K
    S12 = sparse.csr_matrix((ndofs_u, ndofs_u))
    S21 = sparse.csr_matrix((ndofs_u, ndofs_u))

    #lower boundary
    bc_id = np.where(y_u < dx/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)
    S12 = la_utils.clear_rows(S12, bc_id)
    S21 = la_utils.clear_rows(S21, bc_id)

    #upper boundary
    bc_id = np.where(y_u > 1-dx/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)
    S12 = la_utils.clear_rows(S12, bc_id)
    S21 = la_utils.clear_rows(S21, bc_id)

    #left boundary
    bc_id = np.where(x_u < dx/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)
    S12 = la_utils.clear_rows(S12, bc_id)
    S21 = la_utils.clear_rows(S21, bc_id)

    #right boundary
    bc_id = np.where(x_u > 1-dx/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)
    S12 = la_utils.clear_rows(S12, bc_id)
    S21 = la_utils.clear_rows(S21, bc_id)

    #### assembly of Navier-Stokes system
    mat = sparse.vstack([
        sparse.hstack([D11, S12, -BT1, sparse.csr_matrix((ndofs_u, 1))]),
        sparse.hstack([S21, D22, -BT2, sparse.csr_matrix((ndofs_u, 1))]),
        sparse.hstack([-B, sparse.csr_matrix((ndofs_p,ndofs_p)), mean_p.transpose()]),
        sparse.hstack([sparse.csr_matrix((1, 2*ndofs_u)), mean_p, sparse.csr_matrix((1,1))])
    ], "csr")
    return mat

def assemble_blockwise_force_Theta(t,dt):
    size = 2*ndofs_u+ndofs_p+1
    rhs = np.zeros((size,1))

    g_1_now = lambda x,y: analytical_f_1(t, x, y)
    g_2_now = lambda x,y: analytical_f_2(t, x, y)
    g_1_rhs_now = quadrature.fem_rhs(g_1_now, x_u, y_u, topo_u)
    g_2_rhs_now = quadrature.fem_rhs(g_2_now, x_u, y_u, topo_u)

    g_1_prev = lambda x,y: analytical_f_1(t-dt, x, y)
    g_2_prev = lambda x,y: analytical_f_2(t-dt, x, y)
    g_1_rhs_prev = quadrature.fem_rhs(g_1_prev, x_u, y_u, topo_u)
    g_2_rhs_prev = quadrature.fem_rhs(g_2_prev, x_u, y_u, topo_u)
    f_rhs_x = 1./dt*M.dot(u_Theta[0:ndofs_u]) - 0.5*K.dot(u_Theta[0:ndofs_u])
    f_rhs_x +=  0.5*BT1.dot(u_Theta[2*ndofs_u:2*ndofs_u+ndofs_p])
    f_rhs_x += 0.5*np.reshape((g_1_rhs_now + g_1_rhs_prev), (ndofs_u, 1))

    f_rhs_y = 1./dt*M.dot(u_Theta[ndofs_u:2*ndofs_u]) - 0.5*K.dot(u_Theta[ndofs_u:2*ndofs_u])
    f_rhs_y += 0.5*BT2.dot(u_Theta[2*ndofs_u:2*ndofs_u+ndofs_p])
    f_rhs_y += 0.5*np.reshape((g_2_rhs_now + g_2_rhs_prev), (ndofs_u, 1))

    #upper boundary
    bc_id = np.where(y_u > 1-dx/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    #lower boundary
    bc_id = np.where(y_u < dx/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    #right boundary
    bc_id = np.where(x_u > 1-dx/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    #left boundary
    bc_id = np.where(x_u < dx/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    rhs[0:ndofs_u,:] = f_rhs_x
    rhs[ndofs_u:2*ndofs_u,:] = f_rhs_y

    return np.reshape(rhs, (size))

def assemble_blockwise_matrix_Theta(dt):
    D11 = 1./dt*M + 0.5*K
    D22 = 1./dt*M + 0.5*K
    S12 = sparse.csr_matrix((ndofs_u, ndofs_u))
    S21 = sparse.csr_matrix((ndofs_u, ndofs_u))

    #lower boundary
    bc_id = np.where(y_u < dx/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)
    S12 = la_utils.clear_rows(S12, bc_id)
    S21 = la_utils.clear_rows(S21, bc_id)

    #upper boundary
    bc_id = np.where(y_u > 1-dx/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)
    S12 = la_utils.clear_rows(S12, bc_id)
    S21 = la_utils.clear_rows(S21, bc_id)

    #left boundary
    bc_id = np.where(x_u < dx/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)
    S12 = la_utils.clear_rows(S12, bc_id)
    S21 = la_utils.clear_rows(S21, bc_id)

    #right boundary
    bc_id = np.where(x_u > 1-dx/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)
    S12 = la_utils.clear_rows(S12, bc_id)
    S21 = la_utils.clear_rows(S21, bc_id)

    S12 = 0.5*S12
    S21 = 0.5*S21

    #### assembly of Navier-Stokes system
    mat = sparse.vstack([
        sparse.hstack([D11, S12, -0.5*BT1, sparse.csr_matrix((ndofs_u, 1))]),
        sparse.hstack([S21, D22, -0.5*BT2, sparse.csr_matrix((ndofs_u, 1))]),
        sparse.hstack([-0.5*B, sparse.csr_matrix((ndofs_p,ndofs_p)), mean_p.transpose()]),
        sparse.hstack([sparse.csr_matrix((1, 2*ndofs_u)), mean_p, sparse.csr_matrix((1,1))])
    ], "csr")
    return mat

def l2_norm(M, g):
    l2_g = M.dot(g)
    l2_g = np.dot(l2_g.transpose(),g)
    l2_g = np.sqrt(l2_g)
    return l2_g

def apply_bc(g):
    bc_id = np.where(y_u > 1-dx/10)
    g[bc_id] = 0
    bc_id = bc_id + ndofs_u*np.ones(bc_id[0].shape).astype(int)
    g[bc_id] = 0

    bc_id = np.where(y_u < dx/10)
    g[bc_id] = 0
    bc_id = bc_id + ndofs_u*np.ones(bc_id[0].shape).astype(int)
    g[bc_id] = 0

    bc_id = np.where(x_u > 1-dx/10)
    g[bc_id] = 0
    bc_id = bc_id + ndofs_u*np.ones(bc_id[0].shape).astype(int)
    g[bc_id] = 0

    bc_id = np.where(x_u < dx/10)
    g[bc_id] = 0
    bc_id = bc_id + ndofs_u*np.ones(bc_id[0].shape).astype(int)
    g[bc_id] = 0

    return g

n = 140
dx = 1./n

T = 0.5*np.pi
Theta = 0.5
TOL = 1e-8
max_iter = 10

n_runs = 4

t0 = time.time()
(topo_p,x_p,y_p,topo_u,x_u,y_u,c2f) = lin_t3.mesh_t3_iso_t6(n,n,dx,dx)
t1 = time.time()
print('Mesh generation finished')
print('dofs u = ' + str(2*x_u.shape[0]))
print('dofs p = ' + str(x_p.shape[0]))
print('t mesh = ' + str(t1-t0))
ndofs_u = x_u.shape[0]
ndofs_p = x_p.shape[0]

t0 = time.time()
K = assemble.gradu_gradv_p1(topo_u,x_u,y_u)
M = assemble.u_v_p1(topo_u,x_u,y_u)
(BT1,BT2) = assemble.divu_p_p1_iso_p2_p1(topo_p,x_p,y_p,topo_u,x_u,y_u,c2f)
BT = sparse.vstack([BT1,BT2])
B = BT.transpose()

bc_id = np.where(y_u > 1-dx/10)
BT1 = la_utils.clear_rows(BT1,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)

bc_id = np.where(y_u < dx/10)
BT1 = la_utils.clear_rows(BT1,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)

bc_id = np.where(x_u > 1-dx/10)
BT1 = la_utils.clear_rows(BT1,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)

bc_id = np.where(x_u < dx/10)
BT1 = la_utils.clear_rows(BT1,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)

M_2D = sparse.vstack([
    sparse.hstack([M, sparse.csr_matrix(M.shape)]),
    sparse.hstack([sparse.csr_matrix(M.shape), M])
], "csr")
K_2D = sparse.vstack([
    sparse.hstack([K, sparse.csr_matrix(M.shape)]),
    sparse.hstack([sparse.csr_matrix(M.shape), K])
], "csr")

mean_p = np.zeros((1,ndofs_p))
x_l = x_p[topo_p[0,0:3]]
y_l = y_p[topo_p[0,0:3]]
eval_p = np.zeros((0,2))
(phi_dx,phi_dy,phi,omega) = shp.tri_p1(x_l,y_l,eval_p)

for row in topo_p:
    mean_p[0,row] += omega * np.array([1./3.,1./3.,1./3.])

t1 = time.time()
print('Assembled mass, stiffness and pressure matrix')
print('t = ' + str(t1-t0))

err_BDF1 = np.zeros((n_runs))
err_BDF2 = np.zeros((n_runs))
err_Theta = np.zeros((n_runs))
### start loop over different time steps
for t_ind in range(0, n_runs):
    dt = 0.5*T*2**(-t_ind)

    u_0 = analytical(0)
    u_1 = analytical(dt)

    u_BDF1 = u_1.toarray()
    u_BDF2 = u_1.toarray()
    u_BDF2_old = u_0.toarray()
    u_Theta = u_1.toarray()

    N = int(np.round(T/dt+1))
    print('dt = ' + str(dt) + ', ' + str(N) + ' time steps to solve')

    print('solve BDF1')
    M_BDF1 = assemble_blockwise_matrix_BDF1(dt)
    #lu_BDF1 = sp_la.splu(M_BDF1.tocsc())
    spilu = sp_la.spilu(M_BDF1.tocsc(), fill_factor=300, drop_tol=1e-6)
    M_x = lambda x: spilu.solve(x)
    precond = sp_la.LinearOperator((2*ndofs_u+ndofs_p+1, 2*ndofs_u+ndofs_p+1), M_x)
    ### start time loop for dt
    for k in range(2,N):
        t0_BDF1 = time.time()
        rhs_BDF1 = assemble_blockwise_force_BDF1(k*dt, dt)
        #sol = lu_BDF1.solve(rhs_BDF1)
        sol = sp_la.bicgstab(M_BDF1, rhs_BDF1, M=precond, tol=1e-8)[0]
        u_BDF1 = np.reshape(sol, (2*ndofs_u + ndofs_p + 1, 1))
        # print 'error of BDF1 solution for t = ' + str(k*dt) + ': ' + str(np.linalg.norm(u_BDF1[0:2*ndofs_u]-analytical_u(k*dt)))
        t1_BDF1 = time.time()
    gc.collect()

    print('solve BDF2')
    M_BDF2 = assemble_blockwise_matrix_BDF2(dt)
    #lu_BDF2 = sp_la.splu(M_BDF2.tocsc())
    spilu = sp_la.spilu(M_BDF2.tocsc(), fill_factor=300, drop_tol=1e-6)
    M_x = lambda x: spilu.solve(x)
    precond = sp_la.LinearOperator((2*ndofs_u+ndofs_p+1, 2*ndofs_u+ndofs_p+1), M_x)
    ### start time loop for dt
    for k in range(2,N):
        t0_BDF2 = time.time()
        rhs_BDF2 = assemble_blockwise_force_BDF2(k*dt, dt)
        #sol = lu_BDF2.solve(rhs_BDF2)
        sol = sp_la.bicgstab(M_BDF2, rhs_BDF2, M=precond, tol=1e-8)[0]
        u_BDF2_old = u_BDF2
        u_BDF2 = np.reshape(sol, (2*ndofs_u + ndofs_p + 1, 1))
        # print 'error of BDF2 solution for t = ' + str(k*dt) + ': ' + str(np.linalg.norm(u_BDF2[0:2*ndofs_u]-analytical_u(k*dt)))
        t1_BDF2 = time.time()
    gc.collect()

    print('solve Theta')
    M_Theta = assemble_blockwise_matrix_Theta(dt)
    #lu_Theta = sp_la.splu(M_Theta.tocsc())
    spilu = sp_la.spilu(M_Theta.tocsc(), fill_factor=300, drop_tol=1e-6)
    M_x = lambda x: spilu.solve(x)
    precond = sp_la.LinearOperator((2*ndofs_u+ndofs_p+1, 2*ndofs_u+ndofs_p+1), M_x)
    ### start time loop for dt
    for k in range(2,N):
        t0_Theta = time.time()
        rhs_Theta = assemble_blockwise_force_Theta(k*dt, dt)
        #sol = lu_Theta.solve(rhs_Theta)
        sol = sp_la.bicgstab(M_Theta, rhs_Theta, M=precond, tol=1e-8)[0]
        u_Theta = np.reshape(sol, (2*ndofs_u + ndofs_p + 1, 1))
        # print 'error of Theta solution for t = ' + str(k*dt) + ': ' + str(np.linalg.norm(u_Theta[0:2*ndofs_u]-analytical_u(k*dt)))
        t1_Theta = time.time()
    gc.collect()

    ### End of time loops

    f_x = lambda x, y: analytical_u_1(T, x, y)
    f_y = lambda x, y: analytical_u_2(T, x, y)
    e_x = quadrature.l2error_on_mesh(u_BDF1[0:ndofs_u], f_x, x_u, y_u, topo_u, 6)
    e_y = quadrature.l2error_on_mesh(u_BDF1[ndofs_u:2*ndofs_u], f_y, x_u, y_u, topo_u, 6)
    err_BDF1[t_ind] = np.sqrt(e_x**2 + e_y**2)
    e_x = quadrature.l2error_on_mesh(u_BDF2[0:ndofs_u], f_x, x_u, y_u, topo_u, 6)
    e_y = quadrature.l2error_on_mesh(u_BDF2[ndofs_u:2*ndofs_u], f_y, x_u, y_u, topo_u, 6)
    err_BDF2[t_ind] = np.sqrt(e_x**2 + e_y**2)
    e_x = quadrature.l2error_on_mesh(u_Theta[0:ndofs_u], f_x, x_u, y_u, topo_u, 6)
    e_y = quadrature.l2error_on_mesh(u_Theta[ndofs_u:2*ndofs_u], f_y, x_u, y_u, topo_u, 6)
    err_Theta[t_ind] = np.sqrt(e_x**2 + e_y**2)


    print('error BDF1:  ' + str(err_BDF1[t_ind]))
    print('error BDF2:  ' + str(err_BDF2[t_ind]))
    print('error Theta: ' + str(err_Theta[t_ind]))
    if t_ind > 0:
        print('Error decay BDF1:  '+str(err_BDF1[t_ind-1] / err_BDF1[t_ind]))
        print('Error decay BDF2:  '+str(err_BDF2[t_ind-1] / err_BDF2[t_ind]))
        print('Error decay Theta: '+str(err_Theta[t_ind-1] / err_Theta[t_ind]))

    ### End of loop over timesteps
print()
print('------')
print('dx = ' + str(dx))
print('------')

u_x = lambda x, y: analytical_u_1(T, x, y)
u_y = lambda x, y: analytical_u_2(T, x, y)
norm_u_x = quadrature.l2error_on_mesh(np.zeros(u_BDF1.shape), u_x, x_u, y_u, topo_u, 6)
norm_u_y = quadrature.l2error_on_mesh(np.zeros(u_BDF1.shape), u_y, x_u, y_u, topo_u, 6)
norm_u = np.sqrt(norm_u_x**2  + norm_u_y**2)

print('abs. error BDF1:  ' + str(err_BDF1))
print('abs. error BDF2:  ' + str(err_BDF2))
print('abs. error Theta: ' + str(err_Theta))
print('rel. error BDF1:  ' + str(np.divide(err_BDF1, norm_u)))
print('rel. error BDF2:  ' + str(np.divide(err_BDF2, norm_u)))
print('rel. error Theta: ' + str(np.divide(err_Theta, norm_u)))
print('Error decay BDF1:  '+str(np.divide(err_BDF1[0:n_runs-1], err_BDF1[1:n_runs])))
print('Error decay BDF2:  '+str(np.divide(err_BDF2[0:n_runs-1], err_BDF2[1:n_runs])))
print('Error decay Theta: '+str(np.divide(err_Theta[0:n_runs-1], err_Theta[1:n_runs])))
