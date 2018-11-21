#! /bin/env python

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse.linalg as sp_la
from scipy import sparse
import time
import sys

import assemble
import basis_func as shp
import la_utils
import lin_tri_mesh as lin_t3
import viewers

def analytical_u(t):
    analytical_x =  (4*np.sin(t)+5) * x_u**2 * (1-x_u)**2 * (2*y_u - 6*y_u**2 + 4*y_u**3)
    analytical_y = -(4*np.sin(t)+5) * (2*x_u - 6*x_u**2 + 4*x_u**3) * y_u**2 * (1-y_u)**2
    analytical = np.reshape(np.append(analytical_x, analytical_y), (2*ndofs_u, 1))
    return analytical

def analytical_p(t):
    return np.reshape((x_p-0.5), (ndofs_p, 1))

def analytical(t):
    return sparse.vstack([analytical_u(t), analytical_p(t)])


def f(t):
    ## time derivative
    f_x =  (4*np.cos(t)) * x_u**2 * (1-x_u)**2 * (2*y_u - 6*y_u**2 + 4*y_u**3)
    f_y = -(4*np.cos(t)) * (2*x_u - 6*x_u**2 + 4*x_u**3) * y_u**2 * (1-y_u)**2
    ## convection term
    f_x +=  (4*np.sin(t)+5)**2 * x_u**2 * (1-x_u)**2 * (2*y_u - 6*y_u**2 + 4*y_u**3)**2 * (2*x_u - 6*x_u**2 + 4*x_u**3)
    f_y += -(4*np.sin(t)+5)**2 * t**4 * x_u**2 * (1-x_u)**2 * (2*y_u - 6*y_u**2 + 4*y_u**3) * (2 - 12*x_u + 12*x_u**2) * y_u**2 * (1-y_u)**2
    f_x += -(4*np.sin(t)+5)**2 * (2*x_u - 6*x_u**2 + 4*x_u**3) * y_u**2 * (1-y_u)**2 * x_u**2 * (1-x_u)**2 * (2 - 12*y_u + 12*y_u**2)
    f_y +=  (4*np.sin(t)+5)**2 * t**4 * (2*x_u - 6*x_u**2 + 4*x_u**3)**2 * y_u**2 * (1-y_u)**2 * (2*y_u - 6*y_u**2 + 4*y_u**3)
    ## diffusion term
    f_x += -(4*np.sin(t)+5) * ((2 - 12*x_u + 12*x_u**2) * (2*y_u - 6*y_u**2 + 4*y_u**3) + x_u**2 * (1-x_u)**2 * (-12 + 24*y_u))
    f_y += -(4*np.sin(t)+5) * (-(-12 + 24*x_u) * y_u**2 * (1-y_u)**2 - (2*x_u - 6*x_u**2 + 4*x_u**3) * (2 - 12*y_u + 12*y_u**2))
    ## pressure gradient
    f_x +=  -1
    f_y +=  0
    f_stacked = np.reshape(np.append(f_x, f_y), (2*ndofs_u, 1))
    return f_stacked

def assemble_blockwise_force_BDF1(t):
    size = 2*ndofs_u+ndofs_p+1
    rhs = np.zeros((size,1))

    g = f(t)
    f_rhs_x = 1./dt*M.dot(u_BDF1[0:ndofs_u]) + M.dot(g[0:ndofs_u])
    f_rhs_y = 1./dt*M.dot(u_BDF1[ndofs_u:2*ndofs_u]) + M.dot(g[ndofs_u:2*ndofs_u])

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

def assemble_blockwise_matrix_BDF1():
    S11 = assemble.u_gradv_w_p1(topo_u, x_u, y_u, ux_n1, uy_n1)
    D11 = 1./dt*M + K + S11
    D22 = 1./dt*M + K + S11
    # D11 = 1./dt*M + K
    # D22 = 1./dt*M + K
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

def assemble_blockwise_force_BDF2(t):
    size = 2*ndofs_u+ndofs_p+1
    rhs = np.zeros((size,1))

    g = f(t)
    f_rhs_x = M.dot(2*u_BDF2[0:ndofs_u] - 0.5*u_BDF2_old[0:ndofs_u]) + dt*M.dot(g[0:ndofs_u])
    f_rhs_y = M.dot(2*u_BDF2[ndofs_u:2*ndofs_u] - 0.5*u_BDF2_old[ndofs_u:2*ndofs_u]) + dt*M.dot(g[ndofs_u:2*ndofs_u])

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

def assemble_blockwise_matrix_BDF2():
    S11 = assemble.u_gradv_w_p1(topo_u, x_u, y_u, ux_n1, uy_n1)
    D11 = 1.5*M + dt*K + dt*S11
    D22 = 1.5*M + dt*K + dt*S11
    # D11 = 1.5/dt*M + K
    # D22 = 1.5/dt*M + K
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

def assemble_blockwise_force_Theta(t):
    S11 = assemble.u_gradv_w_p1(topo_u, x_u, y_u, ux_n1, uy_n1)
    size = 2*ndofs_u+ndofs_p+1
    rhs = np.zeros((size,1))

    g_now = f(t)
    g_prev = f(t-dt)
    f_rhs_x = 1./dt*M.dot(u_Theta[0:ndofs_u]) - 0.5*K.dot(u_Theta[0:ndofs_u])
    f_rhs_x += - 0.5*S11.dot(u_Theta[0:ndofs_u]) - 0.5*S11.dot(u_Theta[ndofs_u:2*ndofs_u])
    f_rhs_x +=  0.5*BT1.dot(u_Theta[2*ndofs_u:2*ndofs_u+ndofs_p])
    f_rhs_x += 0.5*M.dot(g_now[0:ndofs_u] + g_prev[0:ndofs_u])

    f_rhs_y = 1./dt*M.dot(u_Theta[ndofs_u:2*ndofs_u]) - 0.5*K.dot(u_Theta[ndofs_u:2*ndofs_u])
    f_rhs_y += - 0.5*S11.dot(u_Theta[0:ndofs_u]) - 0.5*S11.dot(u_Theta[ndofs_u:2*ndofs_u])
    f_rhs_y += 0.5*BT2.dot(u_Theta[2*ndofs_u:2*ndofs_u+ndofs_p])
    f_rhs_y += 0.5*M.dot(g_now[ndofs_u:2*ndofs_u] + g_prev[ndofs_u:2*ndofs_u])

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

def assemble_blockwise_matrix_Theta():
    S11 = assemble.u_gradv_w_p1(topo_u, x_u, y_u, ux_n1, uy_n1)
    D11 = 1./dt*M + 0.5*K + 0.5*S11
    D22 = 1./dt*M + 0.5*K + 0.5*S11
    # D11 = 1./dt*M + 0.5*K
    # D22 = 1./dt*M + 0.5*K
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
        sparse.hstack([-B, sparse.csr_matrix((ndofs_p,ndofs_p)), mean_p.transpose()]),
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

if len(sys.argv) > 1:
    n = int(sys.argv[1])
else:
    n = 10
print n
dx = 1./n

T = 5
Theta = 0.5
TOL = 1e-7
max_iter = 10

n_runs = 3

t0 = time.time()
(topo_p,x_p,y_p,topo_u,x_u,y_u,c2f) = lin_t3.mesh_t3_iso_t6(n,n,dx,dx)
t1 = time.time()
print 'Mesh generation finished'
print 'dofs u = ' + str(2*x_u.shape[0])
print 'dofs p = ' + str(x_p.shape[0])
print 't mesh = ' + str(t1-t0)
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
print 'Assembled mass, stiffness and pressure matrix'
print 't = ' + str(t1-t0)

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
    print 'dt = ' + str(dt) + ', ' + str(N) + ' time steps to solve'
    ### start time loop for dt
    for k in range(2,N):
        print 't = ' + str(k*dt)
        t0_BDF1 = time.time()
        ux_n1 = np.reshape(u_BDF1[0:ndofs_u], (ndofs_u, 1))
        uy_n1 = np.reshape(u_BDF1[ndofs_u:2*ndofs_u], (ndofs_u, 1))
        rhs_BDF1 = assemble_blockwise_force_BDF1(k*dt)
        M_BDF1 = assemble_blockwise_matrix_BDF1()
        ### start nonlinear solver for BDF1
        for nonlin_ind in range(max_iter):
            sol = sp_la.spsolve(M_BDF1, rhs_BDF1)
            ux_n1 = np.reshape(sol[0:ndofs_u], (ndofs_u, 1))
            uy_n1 = np.reshape(sol[ndofs_u:2*ndofs_u], (ndofs_u, 1))
            M_BDF1 = assemble_blockwise_matrix_BDF1()
            res = np.linalg.norm(M_BDF1.dot(sol) - rhs_BDF1)
            print 'BDF1, res = ' + str(res)
            if res < TOL:
                break
        u_BDF1 = np.reshape(sol, (2*ndofs_u + ndofs_p + 1, 1))
        print 'error of BDF1 solution for t = ' + str(k*dt) + ': ' + str(l2_norm(M_2D, u_BDF1[0:2*ndofs_u] - analytical_u(k*dt)))
        t1_BDF1 = time.time()

        t0_BDF2 = time.time()
        ux_n1 = np.reshape(u_BDF2[0:ndofs_u], (ndofs_u, 1))
        uy_n1 = np.reshape(u_BDF2[ndofs_u:2*ndofs_u], (ndofs_u, 1))
        rhs_BDF2 = assemble_blockwise_force_BDF2(k*dt)
        M_BDF2 = assemble_blockwise_matrix_BDF2()
        ### start nonlinear solver for BDF2
        for nonlin_ind in range(max_iter):
            sol = sp_la.spsolve(M_BDF2, rhs_BDF2)
            ux_n1 = np.reshape(sol[0:ndofs_u], (ndofs_u, 1))
            uy_n1 = np.reshape(sol[ndofs_u:2*ndofs_u], (ndofs_u, 1))
            M_BDF2 = assemble_blockwise_matrix_BDF2()
            res = np.linalg.norm(M_BDF2.dot(sol) - rhs_BDF2)
            print 'BDF2, res = ' + str(res)
            if res < TOL:
                break
        u_BDF2_old = u_BDF2
        u_BDF2 = np.reshape(sol, (2*ndofs_u + ndofs_p + 1, 1))
        print 'error of BDF2 solution for t = ' + str(k*dt) + ': ' + str(l2_norm(M_2D, u_BDF2[0:2*ndofs_u] - analytical_u(k*dt)))
        t1_BDF2 = time.time()

        t0_Theta = time.time()
        ux_n1 = np.reshape(u_Theta[0:ndofs_u], (ndofs_u, 1))
        uy_n1 = np.reshape(u_Theta[ndofs_u:2*ndofs_u], (ndofs_u, 1))
        rhs_Theta = assemble_blockwise_force_Theta(k*dt)
        M_Theta = assemble_blockwise_matrix_Theta()
        ### Start nonlinear solver for Theta
        for nonlin_ind in range(max_iter):
            sol = sp_la.spsolve(M_Theta, rhs_Theta)
            ux_n1 = np.reshape(sol[0:ndofs_u], (ndofs_u, 1))
            uy_n1 = np.reshape(sol[ndofs_u:2*ndofs_u], (ndofs_u, 1))
            M_Theta = assemble_blockwise_matrix_Theta()
            res = np.linalg.norm(M_Theta.dot(sol) - rhs_Theta)
            print 'Theta, res = ' + str(res)
            if res < TOL:
                break
        u_Theta = np.reshape(sol, (2*ndofs_u + ndofs_p + 1, 1))
        print 'error of Theta solution for t = ' + str(k*dt) + ': ' + str(l2_norm(M_2D, u_Theta[0:2*ndofs_u] - analytical_u(k*dt)))
        t1_Theta = time.time()

        # print l2_norm(M_2D, analytical_u(k*dt))
        # print l2_norm(M_2D, u_BDF1[0:2*ndofs_u])
        #
        # viewers.quiver_vel(x_u,y_u,u_BDF1[0:2*ndofs_u],2*n,2*n,'asdf')
        # plt.show()
        # viewers.quiver_vel(x_u,y_u,analytical_u(k*dt),2*n,2*n,'asdf')
        # plt.show()

        ### End of time loop

    err_BDF1[t_ind] = l2_norm(M_2D, u_BDF1[0:2*ndofs_u]-analytical_u(T))
    err_BDF2[t_ind] = l2_norm(M_2D, u_BDF2[0:2*ndofs_u]-analytical_u(T))
    err_Theta[t_ind] = l2_norm(M_2D, u_Theta[0:2*ndofs_u]-analytical_u(T))
    # err_BDF1[t_ind] = np.linalg.norm(u_BDF1[0:2*ndofs_u]-analytical_u(T))
    # err_BDF2[t_ind] = np.linalg.norm(u_BDF2[0:2*ndofs_u]-analytical_u(T))
    # err_Theta[t_ind] = np.linalg.norm(u_Theta[0:2*ndofs_u]-analytical_u(T))
    print 't BDF1 per step  = ' + str(t1_BDF1-t0_BDF1)
    print 't BDF2 per step  = ' + str(t1_BDF2-t0_BDF2)
    print 't Theta per step = ' + str(t1_Theta-t0_Theta)
    print 'error BDF1:  ' + str(err_BDF1[t_ind])
    print 'error BDF2:  ' + str(err_BDF2[t_ind])
    print 'error Theta: ' + str(err_Theta[t_ind])
    if t_ind > 0:
        print 'Error decay BDF1:  '+str(err_BDF1[t_ind-1] / err_BDF1[t_ind])
        print 'Error decay BDF2:  '+str(err_BDF2[t_ind-1] / err_BDF2[t_ind])
        print 'Error decay Theta: '+str(err_Theta[t_ind-1] / err_Theta[t_ind])

    ### End of loop over timesteps
print
print '------'
print 'dx = ' + str(dx)
print '------'

print 'error BDF1:  ' + str(err_BDF1)
print 'error BDF2:  ' + str(err_BDF2)
print 'error Theta: ' + str(err_Theta)
print 'Error decay BDF1:  '+str(np.divide(err_BDF1[0:n_runs-1], err_BDF1[1:n_runs]))
print 'Error decay BDF2:  '+str(np.divide(err_BDF2[0:n_runs-1], err_BDF2[1:n_runs]))
print 'Error decay Theta: '+str(np.divide(err_Theta[0:n_runs-1], err_Theta[1:n_runs]))
