#! /bin/env python

import numpy as np
import scipy.sparse.linalg as sp_la
from scipy import sparse
import time

import assemble
import la_utils
import lin_tri_mesh as lin_t3

def analytical_u(t):
    analytical_x = np.sin(t) * x_u**2 * (1-x_u)**2 * (2*y_u - 6*y_u**2 + 4*y_u**3)
    analytical_y = -np.sin(t) * (2*x_u - 6*x_u**2 + 4*x_u**3) * y_u**2 * (1-y_u)**2
    analytical = np.reshape(np.append(analytical_x, analytical_y), (2*ndofs_u, 1))
    return analytical

def analytical_p(t):
    return np.reshape(np.sin(t) * (x_p-0.5), (ndofs_p, 1))

def analytical(t):
    return sparse.vstack([analytical_u(t), analytical_p(t)])


def f(t):
    f_x = (2 - 12*x_u + 12*x_u**2) * (2*y_u - 6*y_u**2 + 4*y_u**3) + x_u**2 * (1-x_u)**2 * (-12 + 24*y_u) - 1
    f_y = -(-12 + 24*x_u) * y_u**2 * (1-y_u)**2 - (2*x_u - 6*x_u**2 + 4*x_u**3) * (2 - 12*y_u + 12*y_u**2)
    f_stacked = np.reshape(np.append(f_x, f_y), (2*ndofs_u, 1))
    return np.cos(t) * analytical_u(t) - np.sin(t) * f_stacked

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

n = 80
dx = 1./n

t0 = time.time()
(topo_p,x_p,y_p,topo_u,x_u,y_u,c2f) = lin_t3.mesh_t3_iso_t6(n,n,dx,dx)
t1 = time.time()
print 'Mesh generation finished'
print 'dofs u = ' + str(2*x_u.shape[0])
print 'dofs p = ' + str(x_p.shape[0])
print 't mesh = ' + str(t1-t0)
ndofs_u = x_u.shape[0]
ndofs_p = x_p.shape[0]

T = 8
Theta = 0.5

t0 = time.time()
K = assemble.gradu_gradv_p1(topo_u,x_u,y_u)
M = assemble.u_v_p1(topo_u,x_u,y_u)
(BT1,BT2) = assemble.divu_p_p1_iso_p2_p1(topo_p,x_p,y_p,topo_u,x_u,y_u,c2f)
BT = sparse.vstack([BT1,BT2])
B = BT.transpose()
M_2D = sparse.vstack([
    sparse.hstack([M, sparse.csr_matrix(M.shape)]),
    sparse.hstack([sparse.csr_matrix(M.shape), M])
], "csr")
K_2D = sparse.vstack([
    sparse.hstack([K, sparse.csr_matrix(M.shape)]),
    sparse.hstack([sparse.csr_matrix(M.shape), K])
], "csr")
t1 = time.time()
print 'Assembled mass, stiffness and pressure matrix'
print 't = ' + str(t1-t0)

err_BDF1 = np.zeros((5))
err_BDF2 = np.zeros((5))
err_Theta = np.zeros((5))
for t_ind in range(0, 5):
    dt = 2**(2-t_ind)

    u_0 = analytical(0)
    u_1 = analytical(dt)

    M_BDF1 = M + dt*K
    M_BDF2 = 1.5*M + dt*K
    M_Theta = M + dt*Theta*K

    bc_id = np.where(y_u > 1-dx/10)
    M_BDF1 = la_utils.set_diag(M_BDF1,bc_id)
    M_BDF2 = la_utils.set_diag(M_BDF2,bc_id)
    M_Theta = la_utils.set_diag(M_Theta,bc_id)
    BT1 = la_utils.clear_rows(BT1,bc_id)
    BT2 = la_utils.clear_rows(BT2,bc_id)

    bc_id = np.where(y_u < dx/10)
    M_BDF1 = la_utils.set_diag(M_BDF1,bc_id)
    M_BDF2 = la_utils.set_diag(M_BDF2,bc_id)
    M_Theta = la_utils.set_diag(M_Theta,bc_id)
    BT1 = la_utils.clear_rows(BT1,bc_id)
    BT2 = la_utils.clear_rows(BT2,bc_id)

    bc_id = np.where(x_u > 1-dx/10)
    M_BDF1 = la_utils.set_diag(M_BDF1,bc_id)
    M_BDF2 = la_utils.set_diag(M_BDF2,bc_id)
    M_Theta = la_utils.set_diag(M_Theta,bc_id)
    BT1 = la_utils.clear_rows(BT1,bc_id)
    BT2 = la_utils.clear_rows(BT2,bc_id)

    bc_id = np.where(x_u < dx/10)
    M_BDF1 = la_utils.set_diag(M_BDF1,bc_id)
    M_BDF2 = la_utils.set_diag(M_BDF2,bc_id)
    M_Theta = la_utils.set_diag(M_Theta,bc_id)
    BT1 = la_utils.clear_rows(BT1,bc_id)
    BT2 = la_utils.clear_rows(BT2,bc_id)

    M_BDF1 = sparse.vstack([
        sparse.hstack([M_BDF1, sparse.csr_matrix(M_BDF1.shape), -dt*BT1]),
        sparse.hstack([sparse.csr_matrix(M_BDF1.shape), M_BDF1, -dt*BT2]),
        sparse.hstack([B, sparse.csr_matrix((ndofs_p, ndofs_p))])
    ], "csr")

    M_BDF2 = sparse.vstack([
        sparse.hstack([M_BDF2, sparse.csr_matrix(M_BDF2.shape), -dt*BT1]),
        sparse.hstack([sparse.csr_matrix(M_BDF2.shape), M_BDF2, -dt*BT2]),
        sparse.hstack([B, sparse.csr_matrix((ndofs_p, ndofs_p))])
    ], "csr")

    M_Theta = sparse.vstack([
        sparse.hstack([M_Theta, sparse.csr_matrix(M_Theta.shape), -Theta*dt*BT1]),
        sparse.hstack([sparse.csr_matrix(M_Theta.shape), M_Theta, -Theta*dt*BT2]),
        sparse.hstack([B, sparse.csr_matrix((ndofs_p, ndofs_p))])
    ], "csr")

    f_now = f(0)
    f_old = f(dt)

    u_BDF1 = u_1.toarray()
    u_BDF2 = u_1.toarray()
    u_BDF2_old = u_0.toarray()
    u_Theta = u_1.toarray()

    N = int(np.round(T/dt+1))
    print 'dt = ' + str(dt) + ', ' + str(N) + ' time steps to solve'
    for k in range(2,N):
        t0_BDF1 = time.time()
        f_now = f(k*dt)
        rhs_BDF1 = M_2D.dot(dt*f_now + u_BDF1[0:2*ndofs_u])
        rhs_BDF1 = apply_bc(rhs_BDF1)
        sol = sp_la.spsolve(M_BDF1, np.append(rhs_BDF1, np.zeros((ndofs_p, 1))))
        u_BDF1 = np.reshape(sol, (2*ndofs_u + ndofs_p, 1))
        t1_BDF1 = time.time()

        t0_BDF2 = time.time()
        f_now = f(k*dt)
        rhs_BDF2 = M_2D.dot(dt*f_now + 2*u_BDF2[0:2*ndofs_u] - 0.5*u_BDF2_old[0:2*ndofs_u])
        rhs_BDF2 = apply_bc(rhs_BDF2)
        sol = sp_la.spsolve(M_BDF2, np.append(rhs_BDF2, np.zeros((ndofs_p, 1))))
        u_BDF2_old = u_BDF2
        u_BDF2 = np.reshape(sol, (2*ndofs_u + ndofs_p, 1))
        t1_BDF2 = time.time()

        t0_Theta = time.time()
        f_now = f(k*dt)
        f_old = f((k-1)*dt)
        rhs_Theta = M_2D.dot(dt*Theta*f_now + dt*(1-Theta)*f_old + u_Theta[0:2*ndofs_u])
        rhs_Theta = rhs_Theta - (1-Theta)*dt*K_2D.dot(u_Theta[0:2*ndofs_u])
        rhs_Theta = rhs_Theta + (1-Theta)*dt*BT.dot(u_Theta[2*ndofs_u:2*ndofs_u+ndofs_p])
        rhs_Theta = np.ravel(rhs_Theta)
        rhs_Theta = apply_bc(rhs_Theta)
        sol = sp_la.spsolve(M_Theta, np.append(rhs_Theta, np.zeros((ndofs_p, 1))))
        u_Theta = np.reshape(sol, (2*ndofs_u + ndofs_p, 1))
        t1_Theta = time.time()

        if k*dt == T:
            err_BDF1[t_ind] = np.linalg.norm(u_BDF1[0:2*ndofs_u]-analytical_u(k*dt))
            err_BDF2[t_ind] = np.linalg.norm(u_BDF2[0:2*ndofs_u]-analytical_u(k*dt))
            err_Theta[t_ind] = np.linalg.norm(u_Theta[0:2*ndofs_u]-analytical_u(k*dt))
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

print 'error BDF1:  ' + str(err_BDF1)
print 'error BDF2:  ' + str(err_BDF2)
print 'error Theta: ' + str(err_Theta)
print 'Error decay BDF1:  '+str(np.divide(err_BDF1[0:4], err_BDF1[1:5]))
print 'Error decay BDF2:  '+str(np.divide(err_BDF2[0:4], err_BDF2[1:5]))
print 'Error decay Theta: '+str(np.divide(err_Theta[0:4], err_Theta[1:5]))
