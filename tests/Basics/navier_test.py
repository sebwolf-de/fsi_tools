#! /bin/env python

#import matplotlib.pyplot as plt
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
#from preconditioner import BlockPreconditioner
import quadrature
#import viewers

def analytical_u(t, x, y):
    analytical_x =  (4*np.sin(t)+5) * x**2 * (1-x)**2 * (2*y - 6*y**2 + 4*y**3)
    analytical_y = -(4*np.sin(t)+5) * (2*x - 6*x**2 + 4*x**3) * y**2 * (1-y)**2
    analytical = np.reshape(np.append(analytical_x, analytical_y), (2*len(x), 1))
    return analytical

def analytical_p(t, x_p, y_p):
    p_1 = x_p - 0.5
    p_0 = np.zeros(topo_p.shape[0])
    return np.reshape(np.append(p_1, p_0), (ndofs_p, 1))

def analytical(t, x_u, y_u, x_p, y_p):
    return sparse.vstack([analytical_u(t, x_u, y_u), analytical_p(t, x_p, y_p),0])


def f(t):
    # ## time derivative
    # f_x =  (4*np.cos(t)) * x_u**2 * (1-x_u)**2 * (2*y_u - 6*y_u**2 + 4*y_u**3)
    # f_y = -(4*np.cos(t)) * (2*x_u - 6*x_u**2 + 4*x_u**3) * y_u**2 * (1-y_u)**2
    # ## convection term
    # f_x +=  (4*np.sin(t)+5)**2 * x_u**2 * (1-x_u)**2 * (2*y_u - 6*y_u**2 + 4*y_u**3)**2 * (2*x_u - 6*x_u**2 + 4*x_u**3)
    # f_y += -(4*np.sin(t)+5)**2 * x_u**2 * (1-x_u)**2 * (2*y_u - 6*y_u**2 + 4*y_u**3) * (2 - 12*x_u + 12*x_u**2) * y_u**2 * (1-y_u)**2
    # f_x += -(4*np.sin(t)+5)**2 * (2*x_u - 6*x_u**2 + 4*x_u**3) * y_u**2 * (1-y_u)**2 * x_u**2 * (1-x_u)**2 * (2 - 12*y_u + 12*y_u**2)
    # f_y +=  (4*np.sin(t)+5)**2 * (2*x_u - 6*x_u**2 + 4*x_u**3)**2 * y_u**2 * (1-y_u)**2 * (2*y_u - 6*y_u**2 + 4*y_u**3)
    # ## diffusion term
    # f_x += -(4*np.sin(t)+5) * ((2 - 12*x_u + 12*x_u**2) * (2*y_u - 6*y_u**2 + 4*y_u**3) + x_u**2 * (1-x_u)**2 * (-12 + 24*y_u))
    # f_y += -(4*np.sin(t)+5) * (-(-12 + 24*x_u) * y_u**2 * (1-y_u)**2 - (2*x_u - 6*x_u**2 + 4*x_u**3) * (2 - 12*y_u + 12*y_u**2))
    # ## pressure gradient
    # f_x +=  -1
    # f_y +=  0
    f_x = analytical_f_1(t, x_u, y_u)
    f_y = analytical_f_2s(t, x_u, y_u)
    f_stacked = np.reshape(np.append(f_x, f_y), (2*ndofs_u, 1))
    return f_stacked

def analytical_f_1(t, x, y):
    ndofs = x.shape[0]
    ## time derivative
    f_x =  (4*np.cos(t)) * x**2 * (1-x)**2 * (2*y - 6*y**2 + 4*y**3)
    ## convection term
    f_x +=  (4*np.sin(t)+5)**2 * x**2 * (1-x)**2 * (2*y - 6*y**2 + 4*y**3)**2 * (2*x - 6*x**2 + 4*x**3)
    f_x += -(4*np.sin(t)+5)**2 * (2*x - 6*x**2 + 4*x**3) * y**2 * (1-y)**2 * y**2 * (1-x)**2 * (2 - 12*y + 12*y**2)
    ## diffusion term
    f_x += -(4*np.sin(t)+5) * ((2 - 12*x + 12*x**2) * (2*y - 6*y**2 + 4*y**3) + x**2 * (1-x)**2 * (-12 + 24*y))
    ## pressure gradient
    f_x +=  -1
    return f_x

def analytical_f_2(t, x, y):
    dofs = x.shape[0]
    ## time derivative
    f_y = -(4*np.cos(t)) * (2*x - 6*x**2 + 4*x**3) * y**2 * (1-y)**2
    ## convection term
    f_y += -(4*np.sin(t)+5)**2 * x**2 * (1-x)**2 * (2*y - 6*y**2 + 4*y**3) * (2 - 12*x + 12*x**2) * y**2 * (1-y)**2
    f_y +=  (4*np.sin(t)+5)**2 * (2*x - 6*x**2 + 4*x**3)**2 * y**2 * (1-y)**2 * (2*y - 6*y**2 + 4*y**3)
    ## diffusion term
    f_y += -(4*np.sin(t)+5) * (-(-12 + 24*x) * y**2 * (1-y)**2 - (2*x - 6*x**2 + 4*x**3) * (2 - 12*y + 12*y**2))
    ## pressure gradient
    f_y +=  0
    return f_y

def assemble_blockwise_force_BDF1(t, dt, u_BDF1):
    g_1 = lambda x,y: analytical_f_1(t, x, y)
    g_2 = lambda x,y: analytical_f_2(t, x, y)
    g_1_rhs = quadrature.fem_rhs(g_1, x_u, y_u, topo_u)
    g_2_rhs = quadrature.fem_rhs(g_2, x_u, y_u, topo_u)

    f_rhs_x = M.dot(u_BDF1[0:ndofs_u]) + dt*np.reshape(g_1_rhs, (ndofs_u, 1)) #M.dot(g[0:ndofs_u])
    f_rhs_y = M.dot(u_BDF1[ndofs_u:2*ndofs_u]) + dt*np.reshape(g_2_rhs, (ndofs_u, 1)) #M.dot(g[ndofs_u:2*ndofs_u])

    return apply_rhs_bc(f_rhs_x, f_rhs_y)

def assemble_blockwise_matrix_BDF1(dt, ux_n1, uy_n1):
    S11 = assemble.u_gradv_w_p1(topo_u, x_u, y_u, ux_n1, uy_n1)
    D11 = M + dt*K + dt*S11
    D22 = M + dt*K + dt*S11
    S12 = sparse.csc_matrix((ndofs_u, ndofs_u))
    S21 = sparse.csc_matrix((ndofs_u, ndofs_u))

    (D11, D22) = apply_mat_bc(D11, D22)

    #### assembly of Navier-Stokes system
    mat = sparse.vstack([
        sparse.hstack([D11, S12, -dt*BT1, sparse.csc_matrix((ndofs_u, 1))]),
        sparse.hstack([S21, D22, -dt*BT2, sparse.csc_matrix((ndofs_u, 1))]),
        sparse.hstack([-B, sparse.csc_matrix((ndofs_p,ndofs_p)), mean_p.transpose()]),
        sparse.hstack([sparse.csc_matrix((1, 2*ndofs_u)), mean_p, sparse.csc_matrix((1,1))])
    ], "csc")
    return mat

def assemble_blockwise_force_BDF2(t, dt, u_BDF2, u_BDF2_old):
    g_1 = lambda x,y: analytical_f_1(t, x, y)
    g_2 = lambda x,y: analytical_f_2(t, x, y)
    g_1_rhs = quadrature.fem_rhs(g_1, x_u, y_u, topo_u)
    g_2_rhs = quadrature.fem_rhs(g_2, x_u, y_u, topo_u)

    f_rhs_x = M.dot(2*u_BDF2[0:ndofs_u] - 0.5*u_BDF2_old[0:ndofs_u]) + np.reshape(dt*g_1_rhs, (ndofs_u, 1)) #M.dot(g[0:ndofs_u])
    f_rhs_y = M.dot(2*u_BDF2[ndofs_u:2*ndofs_u] - 0.5*u_BDF2_old[ndofs_u:2*ndofs_u]) + np.reshape(dt*g_2_rhs, (ndofs_u, 1)) #M.dot(g[ndofs_u:2*ndofs_u])

    return apply_rhs_bc(f_rhs_x, f_rhs_y)

def assemble_blockwise_matrix_BDF2(dt, ux_n1, uy_n1):
    S11 = assemble.u_gradv_w_p1(topo_u, x_u, y_u, ux_n1, uy_n1)
    D11 = 1.5*M + dt*K + dt*S11
    D22 = 1.5*M + dt*K + dt*S11
    S12 = sparse.csc_matrix((ndofs_u, ndofs_u))
    S21 = sparse.csc_matrix((ndofs_u, ndofs_u))

    (D11, D22) = apply_mat_bc(D11, D22)

    #### assembly of Navier-Stokes system
    mat = sparse.vstack([
        sparse.hstack([D11, S12, -dt*BT1, sparse.csc_matrix((ndofs_u, 1))]),
        sparse.hstack([S21, D22, -dt*BT2, sparse.csc_matrix((ndofs_u, 1))]),
        sparse.hstack([-B, sparse.csc_matrix((ndofs_p,ndofs_p)), mean_p.transpose()]),
        sparse.hstack([sparse.csc_matrix((1, 2*ndofs_u)), mean_p, sparse.csc_matrix((1,1))])
    ], "csc")
    return mat

def assemble_blockwise_force_Theta(t, dt, u_Theta, ux_n1, uy_n1):
    S11 = assemble.u_gradv_w_p1(topo_u, x_u, y_u, ux_n1, uy_n1)

    g_1_now = lambda x,y: analytical_f_1(t, x, y)
    g_2_now = lambda x,y: analytical_f_2(t, x, y)
    g_1_prev = lambda x,y: analytical_f_1(t-dt, x, y)
    g_2_prev = lambda x,y: analytical_f_2(t-dt, x, y)
    g_1_rhs_now = quadrature.fem_rhs(g_1_now, x_u, y_u, topo_u)
    g_2_rhs_now = quadrature.fem_rhs(g_2_now, x_u, y_u, topo_u)
    g_1_rhs_prev = quadrature.fem_rhs(g_1_prev, x_u, y_u, topo_u)
    g_2_rhs_prev = quadrature.fem_rhs(g_2_prev, x_u, y_u, topo_u)

    f_rhs_x = M.dot(u_Theta[0:ndofs_u]) - dt*0.5*K.dot(u_Theta[0:ndofs_u])
    f_rhs_x += - dt*0.5*S11.dot(u_Theta[0:ndofs_u])
    f_rhs_x +=  dt*0.5*BT1.dot(u_Theta[2*ndofs_u:2*ndofs_u+ndofs_p])
    f_rhs_x += np.reshape(dt*0.5*(g_1_rhs_now + g_1_rhs_prev), (ndofs_u, 1)) #M.dot(g_now[0:ndofs_u] + g_prev[0:ndofs_u])

    f_rhs_y = M.dot(u_Theta[ndofs_u:2*ndofs_u]) - dt*0.5*K.dot(u_Theta[ndofs_u:2*ndofs_u])
    f_rhs_y += - dt*0.5*S11.dot(u_Theta[ndofs_u:2*ndofs_u])
    f_rhs_y += dt*0.5*BT2.dot(u_Theta[2*ndofs_u:2*ndofs_u+ndofs_p])
    f_rhs_y += np.reshape(dt*0.5*(g_2_rhs_now + g_2_rhs_prev), (ndofs_u, 1)) #M.dot(g_now[ndofs_u:2*ndofs_u] + g_prev[ndofs_u:2*ndofs_u])

    return apply_rhs_bc(f_rhs_x, f_rhs_y)

def assemble_blockwise_matrix_Theta(dt, ux_n1, uy_n1):
    S11 = assemble.u_gradv_w_p1(topo_u, x_u, y_u, ux_n1, uy_n1)

    D11 = M + dt*0.5*K + dt*0.5*S11
    D22 = M + dt*0.5*K + dt*0.5*S11
    S12 = sparse.csc_matrix((ndofs_u, ndofs_u))
    S21 = sparse.csc_matrix((ndofs_u, ndofs_u))

    (D11, D22) = apply_mat_bc(D11, D22)

    #### assembly of Navier-Stokes system
    mat = sparse.vstack([
        sparse.hstack([D11, S12, -dt*0.5*BT1, sparse.csc_matrix((ndofs_u, 1))]),
        sparse.hstack([S21, D22, -dt*0.5*BT2, sparse.csc_matrix((ndofs_u, 1))]),
        sparse.hstack([-B, sparse.csc_matrix((ndofs_p,ndofs_p)), mean_p.transpose()]),
        sparse.hstack([sparse.csc_matrix((1, 2*ndofs_u)), mean_p, sparse.csc_matrix((1,1))])
    ], "csc")
    return mat

def l2_norm(M, g):
    l2_g = M.dot(g)
    l2_g = np.dot(l2_g.transpose(),g)
    l2_g = np.sqrt(l2_g)
    return l2_g

def apply_rhs_bc(f_rhs_x, f_rhs_y):
    size = 2*ndofs_u+ndofs_p+1
    rhs = np.zeros((size,1))

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

def apply_mat_bc(D11, D22):
    #lower boundary
    bc_id = np.where(y_u < dx/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)

    #upper boundary
    bc_id = np.where(y_u > 1-dx/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)

    #left boundary
    bc_id = np.where(x_u < dx/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)

    #right boundary
    bc_id = np.where(x_u > 1-dx/10)
    D11 = la_utils.set_diag(D11, bc_id)
    D22 = la_utils.set_diag(D22, bc_id)

    return D11, D22

def l2_norm(M,g):
    l2_g = M.dot(g)
    l2_g = np.dot(l2_g.transpose(),g)
    l2_g = np.sqrt(l2_g)
    return l2_g

if len(sys.argv) > 1:
    n = int(sys.argv[1])
else:
    n = 10
print(n)
dx = 1./n

T = 5
Theta = 0.5
TOL = 1e-8
max_iter = 10

n_runs = 5

(topo_p,x_p,y_p,topo_u,x_u,y_u,c2f) = lin_t3.mesh_t3_iso_t6(n,n,dx,dx)
(topo_p,x_p,y_p) = lin_t3.mesh_t3_t0(n,n,dx,dx)
print('Mesh generation finished')

t0 = time.time()
K = assemble.gradu_gradv_p1(topo_u,x_u,y_u)
M = assemble.u_v_p1(topo_u,x_u,y_u)
(BT1,BT2) = assemble.divu_p_p1_iso_p2_p1p0(topo_p,x_p,y_p,topo_u,x_u,y_u,c2f)
BT = sparse.vstack([BT1,BT2])
B = BT.transpose()

ndofs_u = BT1.shape[0]
ndofs_p = BT1.shape[1]
print('dofs u = ' + str(2*ndofs_u))
print('dofs p = ' + str(ndofs_p))

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
    sparse.hstack([M, sparse.csc_matrix(M.shape)]),
    sparse.hstack([sparse.csc_matrix(M.shape), M])
], "csc")
K_2D = sparse.vstack([
    sparse.hstack([K, sparse.csc_matrix(M.shape)]),
    sparse.hstack([sparse.csc_matrix(M.shape), K])
], "csc")

mean_p = np.zeros((1,ndofs_p))
x_l = x_p[topo_p[0,0:3]]
y_l = y_p[topo_p[0,0:3]]
eval_p = np.zeros((0,2))
(phi_dx,phi_dy,phi,omega) = shp.tri_p1(x_l,y_l,eval_p)

for row in topo_p:
    mean_p[0,row] += omega * np.array([1./3.,1./3.,1./3., 1.])

nodes_p = x_p.shape[0]
cells_p = topo_p.shape[0]
MP = sparse.vstack([
    sparse.hstack([assemble.u_v_p1(topo_p[:,0:3],x_p,y_p), sparse.csr_matrix((nodes_p, cells_p))]),
    sparse.hstack([sparse.csr_matrix((cells_p, nodes_p)), omega*sparse.eye(cells_p)])
    ])

big_mass_matrix = sparse.vstack([
    sparse.hstack([M_2D, sparse.csr_matrix((2*ndofs_u, ndofs_p + 1))]),
    sparse.hstack([sparse.csr_matrix((ndofs_p, 2*ndofs_u)), MP, sparse.csr_matrix((ndofs_p, 1))]),
    sparse.hstack([sparse.csr_matrix((1, 2*ndofs_u+ndofs_p)), sparse.eye(1)])
])
t1 = time.time()

print('Assembled mass, stiffness and pressure matrix, t = ' + str(t1-t0))

err_BDF1 = np.zeros((n_runs))
err_BDF2 = np.zeros((n_runs))
err_Theta = np.zeros((n_runs))
err_BDF1_ref = np.zeros((n_runs))
err_BDF2_ref = np.zeros((n_runs))
err_Theta_ref = np.zeros((n_runs))
BDF1 = np.zeros((2*ndofs_u, n_runs))
BDF2 = np.zeros((2*ndofs_u, n_runs))
Theta = np.zeros((2*ndofs_u, n_runs))
ref = np.zeros((2*ndofs_u))

nonlin_conv_ind = np.zeros((3, n_runs))

### calculate reference solution
dt_ref = 0.01
# dt = dt_ref
# N = int(np.round(T/dt+1))
# print('Calculate reference solution')
# print('dt = ' + str(dt) + ', ' + str(N) + ' time steps to solve')
#
# u_0 = analytical(0, x_u, y_u, x_p, y_p)
# u_1 = analytical(dt, x_u, y_u, x_p, y_p)
#
# u_Theta = u_1.toarray()
#
# ### start time loop for reference solution
# for k in range(2,N):
#     print('t = ' + str(k*dt))
#
#     f_x = lambda x, y: analytical_u(k*dt, x, y)[0:len(x)]
#     f_y = lambda x, y: analytical_u(k*dt, x, y)[len(x):2*len(x)]
#
#     t0_Theta = time.time()
#     ux_n1 = np.reshape(u_Theta[0:ndofs_u], (ndofs_u, 1))
#     uy_n1 = np.reshape(u_Theta[ndofs_u:2*ndofs_u], (ndofs_u, 1))
#     assemble_t0 = time.time()
#     rhs_Theta = assemble_blockwise_force_Theta(k*dt)
#     M_Theta = assemble_blockwise_matrix_Theta()
#     assemble_t1 = time.time()
#     ### Start nonlinear solver for Theta
#     for nonlin_ind in range(max_iter):
#         precond_t0 = time.time()
#         spilu = sp_la.spilu(M_Theta, fill_factor=300, drop_tol=1e-6)
#         M_x = lambda x: spilu.solve(x)
#         precond = sp_la.LinearOperator((2*ndofs_u+ndofs_p+1, 2*ndofs_u+ndofs_p+1), M_x)
#         precond_t1 = time.time()
#         sol_t0 = time.time()
#         sol = sp_la.bicgstab(M_Theta, rhs_Theta, M=precond, tol=1e-8)[0]
#         sol_t1 = time.time()
#         ux_n1 = np.reshape(sol[0:ndofs_u], (ndofs_u, 1))
#         uy_n1 = np.reshape(sol[ndofs_u:2*ndofs_u], (ndofs_u, 1))
#         res_t0 = time.time()
#         M_Theta = assemble_blockwise_matrix_Theta()
#         res = l2_norm(big_mass_matrix, M_Theta.dot(sol) - rhs_Theta)
#         res_t1 = time.time()
#         print('reference, res = ' + str(res))
#         if res < TOL:
#             break
#         ### End of nonlinear solver
#     u_Theta = np.reshape(sol, (2*ndofs_u + ndofs_p + 1, 1))
#
#     ### End of time loop
#
# ref = u_Theta
ref = np.zeros((2*ndofs_u + ndofs_p + 1, 1))

### End of computation of reference solution



### start loop over different time steps
for t_ind in range(0, n_runs):
    dt = 0.5*T*2**(-t_ind)
    N = int(np.round(T/dt+1))
    print('dt = ' + str(dt) + ', ' + str(N) + ' time steps to solve')

    u_0 = analytical(0, x_u, y_u, x_p, y_p)
    u_1 = analytical(dt, x_u, y_u, x_p, y_p)

    u_BDF1 = u_1.toarray()
    u_BDF2 = u_1.toarray()
    u_BDF2_old = u_0.toarray()
    u_Theta = u_1.toarray()

    ### start time loop for dt
    for k in range(2,N):
        print('t = ' + str(k*dt))

        f_x = lambda x, y: analytical_u(k*dt, x, y)[0:len(x)]
        f_y = lambda x, y: analytical_u(k*dt, x, y)[len(x):2*len(x)]

        t0_BDF1 = time.time()
        ux_n1 = np.reshape(u_BDF1[0:ndofs_u], (ndofs_u, 1))
        uy_n1 = np.reshape(u_BDF1[ndofs_u:2*ndofs_u], (ndofs_u, 1))
        assemble_t0 = time.time()
        rhs_BDF1 = assemble_blockwise_force_BDF1(k*dt, dt, u_BDF1)
        M_BDF1 = assemble_blockwise_matrix_BDF1(dt, ux_n1, uy_n1).tocsc()
        assemble_t1 = time.time()
        # print('assembled linear system, t = ' + str(assemble_t1-assemble_t0))
        ### start nonlinear solver for BDF1
        for nonlin_ind in range(max_iter):
            precond_t0 = time.time()
            spilu = sp_la.spilu(M_BDF1, fill_factor=300, drop_tol=1e-6)
            M_x = lambda x: spilu.solve(x)
            precond = sp_la.LinearOperator((2*ndofs_u+ndofs_p+1, 2*ndofs_u+ndofs_p+1), M_x)
            precond_t1 = time.time()
            # print('calculated preconditioner, t = ' + str(precond_t1-precond_t0))
            solve_t0 = time.time()
            sol = sp_la.bicgstab(M_BDF1, rhs_BDF1, M=precond, tol=1e-8)[0]
            ux_n1 = np.reshape(sol[0:ndofs_u], (ndofs_u, 1))
            uy_n1 = np.reshape(sol[ndofs_u:2*ndofs_u], (ndofs_u, 1))
            solve_t1 = time.time()
            # print('solved linear system, t = ' + str(solve_t1 - solve_t0))
            res_t0 = time.time()
            M_BDF1 = assemble_blockwise_matrix_BDF1(dt, ux_n1, uy_n1)
            res = l2_norm(big_mass_matrix, M_BDF1.dot(sol) - rhs_BDF1)
            res_t1 = time.time()
            # print('calculated residual, t = ' + str(res_t1 - res_t0))
            print('BDF1, res = ' + str(res))
            if res < TOL:
                # store the maximum number of nonlinear iterations
                nonlin_conv_ind[0,t_ind] = np.max(nonlin_ind+1, nonlin_conv_ind[0,t_ind])
                break
            if(nonlin_ind == max_iter-1):
                # nonlinear iterator did not converge
                nonlin_conv_ind[0,t_ind] = -1

        u_BDF1 = np.reshape(sol, (2*ndofs_u + ndofs_p + 1, 1))

        e_x = quadrature.l2error_on_mesh(u_BDF1[0:ndofs_u], f_x, x_u, y_u, topo_u, 6)
        e_y = quadrature.l2error_on_mesh(u_BDF1[ndofs_u:2*ndofs_u], f_y, x_u, y_u, topo_u, 6)
        print('error of BDF1 solution for t = ' + str(k*dt) + ': ' + str(np.sqrt(e_x**2 + e_y**2)))

        t1_BDF1 = time.time()

        t0_BDF2 = time.time()
        ux_n1 = np.reshape(u_BDF2[0:ndofs_u], (ndofs_u, 1))
        uy_n1 = np.reshape(u_BDF2[ndofs_u:2*ndofs_u], (ndofs_u, 1))
        assemble_t0 = time.time()
        rhs_BDF2 = assemble_blockwise_force_BDF2(k*dt, dt, u_BDF2, u_BDF2_old)
        M_BDF2 = assemble_blockwise_matrix_BDF2(dt, ux_n1, uy_n1)
        assemble_t1 = time.time()
        # print('assembled linear system, t = ' + str(assemble_t1 - assemble_t0))
        ### start nonlinear solver for BDF2
        for nonlin_ind in range(max_iter):
            precond_t0 = time.time()
            spilu = sp_la.spilu(M_BDF2, fill_factor=300, drop_tol=1e-6)
            M_x = lambda x: spilu.solve(x)
            precond = sp_la.LinearOperator((2*ndofs_u+ndofs_p+1, 2*ndofs_u+ndofs_p+1), M_x)
            precond_t1 = time.time()
            # print('calculated preconditioner, t = ' + str(precond_t1 - precond_t0))
            sol_t0 = time.time()
            sol = sp_la.bicgstab(M_BDF2, rhs_BDF2, M=precond, tol=1e-8)[0]
            sol_t1 = time.time()
            # print('solved linear system, t = ' + str(sol_t1 - sol_t0))
            ux_n1 = np.reshape(sol[0:ndofs_u], (ndofs_u, 1))
            uy_n1 = np.reshape(sol[ndofs_u:2*ndofs_u], (ndofs_u, 1))
            res_t0 = time.time()
            M_BDF2 = assemble_blockwise_matrix_BDF2(dt, ux_n1, uy_n1)
            res = l2_norm(big_mass_matrix, M_BDF2.dot(sol) - rhs_BDF2)
            res_t1 = time.time()
            # print('calculated residual, t = ' + str(res_t1 - res_t0))
            print('BDF2, res = ' + str(res))
            if res < TOL:
                # store the maximum number of nonlinear iterations
                nonlin_conv_ind[1,t_ind] = np.max(nonlin_ind+1, nonlin_conv_ind[1,t_ind])
                break
            if(nonlin_ind == max_iter-1):
                # nonlinear iterator did not converge
                nonlin_conv_ind[1,t_ind] = -1

        u_BDF2_old = u_BDF2
        u_BDF2 = np.reshape(sol, (2*ndofs_u + ndofs_p + 1, 1))

        e_x = quadrature.l2error_on_mesh(u_BDF2[0:ndofs_u], f_x, x_u, y_u, topo_u, 6)
        e_y = quadrature.l2error_on_mesh(u_BDF2[ndofs_u:2*ndofs_u], f_y, x_u, y_u, topo_u, 6)
        print('error of BDF2 solution for t = ' + str(k*dt) + ': ' + str(np.sqrt(e_x**2 + e_y**2)))

        t1_BDF2 = time.time()

        t0_Theta = time.time()
        ux_n1 = np.reshape(u_Theta[0:ndofs_u], (ndofs_u, 1))
        uy_n1 = np.reshape(u_Theta[ndofs_u:2*ndofs_u], (ndofs_u, 1))
        assemble_t0 = time.time()
        rhs_Theta = assemble_blockwise_force_Theta(k*dt, dt, u_Theta, ux_n1, uy_n1)
        M_Theta = assemble_blockwise_matrix_Theta(dt, ux_n1, uy_n1)
        assemble_t1 = time.time()
        # print('assembled linear system, t = ' + str(assemble_t1 - assemble_t0))
        ### Start nonlinear solver for Theta
        for nonlin_ind in range(max_iter):
            precond_t0 = time.time()
            spilu = sp_la.spilu(M_Theta, fill_factor=300, drop_tol=1e-6)
            M_x = lambda x: spilu.solve(x)
            precond = sp_la.LinearOperator((2*ndofs_u+ndofs_p+1, 2*ndofs_u+ndofs_p+1), M_x)
            precond_t1 = time.time()
            # print('calculated preconditioner, t = ' + str(precond_t1 - precond_t0))
            sol_t0 = time.time()
            sol = sp_la.bicgstab(M_Theta, rhs_Theta, M=precond, tol=1e-8)[0]
            sol_t1 = time.time()
            # print('solved linear system, t  = ' + str(sol_t1 - sol_t0))
            ux_n1 = np.reshape(sol[0:ndofs_u], (ndofs_u, 1))
            uy_n1 = np.reshape(sol[ndofs_u:2*ndofs_u], (ndofs_u, 1))
            res_t0 = time.time()
            M_Theta = assemble_blockwise_matrix_Theta(dt, ux_n1, uy_n1)
            res = l2_norm(big_mass_matrix, M_Theta.dot(sol) - rhs_Theta)
            res_t1 = time.time()
            # print('calculated residual, t = ' + str(res_t1 - res_t0))
            print('Theta, res = ' + str(res))
            if res < TOL:
                # store the maximum number of nonlinear iterations
                nonlin_conv_ind[2,t_ind] = np.max(nonlin_ind+1, nonlin_conv_ind[2,t_ind])
                break
            if(nonlin_ind == max_iter-1):
                # nonlinear iterator did not converge
                nonlin_conv_ind[2,t_ind] = -1

        u_Theta = np.reshape(sol, (2*ndofs_u + ndofs_p + 1, 1))

        e_x = quadrature.l2error_on_mesh(u_Theta[0:ndofs_u], f_x, x_u, y_u, topo_u, 6)
        e_y = quadrature.l2error_on_mesh(u_Theta[ndofs_u:2*ndofs_u], f_y, x_u, y_u, topo_u, 6)
        print('error of Theta solution for t = ' + str(k*dt) + ': ' + str(np.sqrt(e_x**2 + e_y**2)))
        t1_Theta = time.time()

        ### End of time loop

    ###Compute a lot of errors
    BDF1[:,t_ind] = u_BDF1[0:2*ndofs_u].ravel()
    BDF2[:,t_ind] = u_BDF2[0:2*ndofs_u].ravel()
    Theta[:,t_ind] = u_Theta[0:2*ndofs_u].ravel()
    f_x = lambda x, y: analytical_u(T, x, y)[0:len(x)]
    f_y = lambda x, y: analytical_u(T, x, y)[len(x):2*len(x)]
    zero_fun = lambda x, y: np.zeros(x.shape)
    norm_u_x = quadrature.l2error_on_mesh(np.zeros(u_BDF1.shape), f_x, x_u, y_u, topo_u, 6)
    norm_u_y = quadrature.l2error_on_mesh(np.zeros(u_BDF1.shape), f_y, x_u, y_u, topo_u, 6)
    norm_u = np.sqrt(norm_u_x**2  + norm_u_y**2)
    e_x = quadrature.l2error_on_mesh(u_BDF1[0:ndofs_u], f_x, x_u, y_u, topo_u, 6)
    e_y = quadrature.l2error_on_mesh(u_BDF1[ndofs_u:2*ndofs_u], f_y, x_u, y_u, topo_u, 6)
    err_BDF1[t_ind] = np.sqrt(e_x**2 + e_y**2)
    e_x = quadrature.l2error_on_mesh(u_BDF1[0:ndofs_u] - ref[0:ndofs_u], zero_fun, x_u, y_u, topo_u, 6)
    e_y = quadrature.l2error_on_mesh(u_BDF1[ndofs_u:2*ndofs_u] - ref[ndofs_u:2*ndofs_u], zero_fun, x_u, y_u, topo_u, 6)
    err_BDF1_ref[t_ind] = np.sqrt(e_x**2 + e_y**2)
    e_x = quadrature.l2error_on_mesh(u_BDF2[0:ndofs_u], f_x, x_u, y_u, topo_u, 6)
    e_y = quadrature.l2error_on_mesh(u_BDF2[ndofs_u:2*ndofs_u], f_y, x_u, y_u, topo_u, 6)
    err_BDF2[t_ind] = np.sqrt(e_x**2 + e_y**2)
    e_x = quadrature.l2error_on_mesh(u_BDF2[0:ndofs_u] - ref[0:ndofs_u], zero_fun, x_u, y_u, topo_u, 6)
    e_y = quadrature.l2error_on_mesh(u_BDF2[ndofs_u:2*ndofs_u] - ref[ndofs_u:2*ndofs_u], zero_fun, x_u, y_u, topo_u, 6)
    err_BDF2_ref[t_ind] = np.sqrt(e_x**2 + e_y**2)
    e_x = quadrature.l2error_on_mesh(u_Theta[0:ndofs_u], f_x, x_u, y_u, topo_u, 6)
    e_y = quadrature.l2error_on_mesh(u_Theta[ndofs_u:2*ndofs_u], f_y, x_u, y_u, topo_u, 6)
    err_Theta[t_ind] = np.sqrt(e_x**2 + e_y**2)
    e_x = quadrature.l2error_on_mesh(u_Theta[0:ndofs_u] - ref[0:ndofs_u], zero_fun, x_u, y_u, topo_u, 6)
    e_y = quadrature.l2error_on_mesh(u_Theta[ndofs_u:2*ndofs_u] - ref[ndofs_u:2*ndofs_u], zero_fun, x_u, y_u, topo_u, 6)
    err_Theta_ref[t_ind] = np.sqrt(e_x**2 + e_y**2)

    # err_BDF1[t_ind] = np.linalg.norm(u_BDF1[0:2*ndofs_u]-analytical_u(T, x_u, y_u, x_p, y_p))
    # err_BDF2[t_ind] = np.linalg.norm(u_BDF2[0:2*ndofs_u]-analytical_u(T, x_u, y_u, x_p, y_p))
    # err_Theta[t_ind] = np.linalg.norm(u_Theta[0:2*ndofs_u]-analytical_u(T, x_u, y_u, x_p, y_p))
    # print('t BDF1 per step  = ' + str(t1_BDF1-t0_BDF1))
    # print('t BDF2 per step  = ' + str(t1_BDF2-t0_BDF2))
    # print('t Theta per step = ' + str(t1_Theta-t0_Theta))
    # print('error BDF1:  ' + str(err_BDF1[t_ind]))
    # print('error BDF2:  ' + str(err_BDF2[t_ind]))
    # print('error Theta: ' + str(err_Theta[t_ind]))
    # if t_ind > 0:
    #     print('Error decay BDF1:  '+str(err_BDF1[t_ind-1] / err_BDF1[t_ind]))
    #     print('Error decay BDF2:  '+str(err_BDF2[t_ind-1] / err_BDF2[t_ind]))
    #     print('Error decay Theta: '+str(err_Theta[t_ind-1] / err_Theta[t_ind]))

    ### End of loop over timesteps

print()
print('------')
print('dx = ' + str(dx))
print('dt = ' + str(0.5*T*np.power(2., -np.arange(0,n_runs))))
print('dt_ref = ' + str(dt_ref))
print('------')

print('error comparted to analytical solution')
print('abs. error BDF1:  ' + str(err_BDF1))
print('abs. error BDF2:  ' + str(err_BDF2))
print('abs. error Theta: ' + str(err_Theta))
print('rel. error BDF1:  ' + str(np.divide(err_BDF1, norm_u)))
print('rel. error BDF2:  ' + str(np.divide(err_BDF2, norm_u)))
print('rel. error Theta: ' + str(np.divide(err_Theta, norm_u)))
print('Error decay BDF1:  '+str(np.divide(err_BDF1[0:n_runs-1], err_BDF1[1:n_runs])))
print('Error decay BDF2:  '+str(np.divide(err_BDF2[0:n_runs-1], err_BDF2[1:n_runs])))
print('Error decay Theta: '+str(np.divide(err_Theta[0:n_runs-1], err_Theta[1:n_runs])))

print()

print('error compared to reference solution')
print('abs. error BDF1:  ' + str(err_BDF1_ref))
print('abs. error BDF2:  ' + str(err_BDF2_ref))
print('abs. error Theta: ' + str(err_Theta_ref))
print('rel. error BDF1:  ' + str(np.divide(err_BDF1_ref, norm_u)))
print('rel. error BDF2:  ' + str(np.divide(err_BDF2_ref, norm_u)))
print('rel. error Theta: ' + str(np.divide(err_Theta_ref, norm_u)))
print('Error decay BDF1:  '+str(np.divide(err_BDF1_ref[0:n_runs-1], err_BDF1_ref[1:n_runs])))
print('Error decay BDF2:  '+str(np.divide(err_BDF2_ref[0:n_runs-1], err_BDF2_ref[1:n_runs])))
print('Error decay Theta: '+str(np.divide(err_Theta_ref[0:n_runs-1], err_Theta_ref[1:n_runs])))

print()

rate_u_BDF1 = np.zeros(n_runs-2)
rate_u_BDF2 = np.zeros(n_runs-2)
rate_u_Theta = np.zeros(n_runs-2)
for k in range(0,n_runs-2):
    rate_u_BDF1[k] = np.log2(l2_norm(M_2D, BDF1[:,k] - BDF1[:,k+1]) / l2_norm(M_2D, BDF1[:,k+1] - BDF1[:,k+2]))
    rate_u_BDF2[k] = np.log2(l2_norm(M_2D, BDF2[:,k] - BDF2[:,k+1]) / l2_norm(M_2D, BDF2[:,k+1] - BDF2[:,k+2]))
    rate_u_Theta[k] = np.log2(l2_norm(M_2D, Theta[:,k] - Theta[:,k+1]) / l2_norm(M_2D, Theta[:,k+1] - Theta[:,k+2]))
print('Empirical rate BDF1: ' + str(rate_u_BDF1))
print('Empirical rate BDF2: ' + str(rate_u_BDF2))
print('Empirical rate Theta: ' + str(rate_u_Theta))

print('Max iterations of the fixpoint iterations:')
print(nonlin_conv_ind)

filename = 'dx='+str(dx)
f = open(filename,"wb")
np.save(f,dx)
np.save(f,0.5*T*np.power(2., -np.arange(0,n_runs)))
np.save(f,dt_ref)
np.save(f,BDF1)
np.save(f,BDF2)
np.save(f,Theta)
np.save(f,ref)
f.close()
