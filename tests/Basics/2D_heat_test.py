#! /bin/env python

import numpy as np
import scipy.sparse.linalg as sp_la
import sys
import time

sys.path.append('../../modules')
import assemble
import la_utils
import lin_tri_mesh as lin_t3

def analytical(t):
    return np.sin(0.5*np.pi*x) * np.sin(np.pi*y) * np.cos(t)

def f(t):
    return np.sin(0.5*np.pi*x) * np.sin(np.pi*y) * (-np.sin(t) + 1.25*np.pi**2*np.cos(t))

n = 200
dx = 1./n

t0 = time.time()
(topo,x,y) = lin_t3.mesh_t3(n,n,dx,dx)
t1 = time.time()
print('Mesh generation finished')
print('dofs   = ' + str(x.shape[0]))
print('t mesh = ' + str(t1-t0))

T = 8
Theta = 0.5

t0 = time.time()
K = assemble.gradu_gradv_p1(topo,x,y)
M = assemble.u_v_p1(topo,x,y)
t1 = time.time()
print('Assembled mass and stiffness matrix')
print('t = ' + str(t1-t0))


err_BDF1 = np.zeros((5))
err_BDF2 = np.zeros((5))
err_Theta = np.zeros((5))
for t_ind in range(0, 5):
    dt = 2**(-t_ind)

    u_0 = analytical(0)
    u_1 = analytical(dt)

    M_BDF1 = M + dt*K
    M_BDF2 = 1.5*M + dt*K
    M_Theta = M + dt*Theta*K

    bc_id = np.where(y > 1-dx/10)
    M_BDF1 = la_utils.set_diag(M_BDF1,bc_id)
    M_BDF2 = la_utils.set_diag(M_BDF2,bc_id)
    M_Theta = la_utils.set_diag(M_Theta,bc_id)

    bc_id = np.where(y < dx/10)
    M_BDF1 = la_utils.set_diag(M_BDF1,bc_id)
    M_BDF2 = la_utils.set_diag(M_BDF2,bc_id)
    M_Theta = la_utils.set_diag(M_Theta,bc_id)

    bc_id = np.where(x > 1-dx/10)
    M_BDF1 = la_utils.set_diag(M_BDF1,bc_id)
    M_BDF2 = la_utils.set_diag(M_BDF2,bc_id)
    M_Theta = la_utils.set_diag(M_Theta,bc_id)

    bc_id = np.where(x < dx/10)
    M_BDF1 = la_utils.set_diag(M_BDF1,bc_id)
    M_BDF2 = la_utils.set_diag(M_BDF2,bc_id)
    M_Theta = la_utils.set_diag(M_Theta,bc_id)

    f_now = f(0)
    f_old = f(dt)

    rhs_BDF1 = M.dot(dt*f_now + u_0)
    rhs_Theta = M.dot(dt*Theta*f_now + dt*(1-Theta)*f_old + u_0) - (1-Theta)*dt*K.dot(u_0)
    rhs_Theta = np.ravel(rhs_Theta)

    bc_id = np.where(y > 1-dx/10)
    rhs_BDF1[bc_id] = 0
    rhs_Theta[bc_id] = 0

    bc_id = np.where(y < dx/10)
    rhs_BDF1[bc_id] = 0
    rhs_Theta[bc_id] = 0

    bc_id = np.where(x > 1-dx/10)
    rhs_BDF1[bc_id] = np.sin(np.pi*y[bc_id]) * np.cos(dt)
    rhs_Theta[bc_id] = np.sin(np.pi*y[bc_id]) * np.cos(dt)

    bc_id = np.where(x < dx/10)
    rhs_BDF1[bc_id] = 0
    rhs_Theta[bc_id] = 0

    u_BDF1 = sp_la.spsolve(M_BDF1, rhs_BDF1)
    u_Theta = sp_la.spsolve(M_Theta, rhs_Theta)

    u_BDF2 = u_Theta
    u_BDF2_old = u_0

    N = int(np.round(T/dt+1))
    for k in range(2,N):
        f_now = f(k*dt)
        f_old = f((k-1)*dt)

        rhs_BDF1 = M.dot(dt*f_now + u_BDF1)
        rhs_BDF2 = M.dot(dt*f_now + 2*u_BDF2 - 0.5*u_BDF2_old)
        rhs_Theta = M.dot(dt*Theta*f_now + dt*(1-Theta)*f_old + u_Theta) - (1-Theta)*dt*K.dot(u_Theta)
        rhs_Theta = np.ravel(rhs_Theta)

        bc_id = np.where(y > 1-dx/10)
        rhs_BDF1[bc_id] = 0
        rhs_BDF2[bc_id] = 0
        rhs_Theta[bc_id] = 0

        bc_id = np.where(y < dx/10)
        rhs_BDF1[bc_id] = 0
        rhs_BDF2[bc_id] = 0
        rhs_Theta[bc_id] = 0

        bc_id = np.where(x > 1-dx/10)
        rhs_BDF1[bc_id] = np.sin(np.pi*y[bc_id]) * np.cos(k*dt)
        rhs_BDF2[bc_id] = np.sin(np.pi*y[bc_id]) * np.cos(k*dt)
        rhs_Theta[bc_id] = np.sin(np.pi*y[bc_id]) * np.cos(k*dt)

        bc_id = np.where(x < dx/10)
        rhs_BDF1[bc_id] = 0
        rhs_BDF2[bc_id] = 0
        rhs_Theta[bc_id] = 0

        t0_BDF1 = time.time()
        u_BDF1 = sp_la.spsolve(M_BDF1, rhs_BDF1)
        t1_BDF1 = time.time()
        t0_BDF2 = time.time()
        sol = sp_la.spsolve(M_BDF2, rhs_BDF2)
        u_BDF2_old = u_BDF2
        u_BDF2 = sol
        t1_BDF2 = time.time()
        t0_Theta = time.time()
        u_Theta = sp_la.spsolve(M_Theta, rhs_Theta)
        t1_Theta = time.time()

        # print 'u_n('+str(k*dt)+') = '+str(u_n)
        # print 'u('+str(k*dt)+') = '+str(analytical(k*dt))
        if k*dt == T:
            err_BDF1[t_ind] = np.linalg.norm(u_BDF1-analytical(k*dt))
            err_BDF2[t_ind] = np.linalg.norm(u_BDF2-analytical(k*dt))
            err_Theta[t_ind] = np.linalg.norm(u_Theta-analytical(k*dt))
            print('dt = ' + str(dt))
            print('t BDF1  = ' + str(t1_BDF1-t0_BDF1))
            print('t BDF2  = ' + str(t1_BDF2-t0_BDF2))
            print('t Theta = ' + str(t1_Theta-t0_Theta))
        # print 'error Theta = '+str(np.linalg.norm(u_Theta-analytical(k*dt)))

print('error BDF1:  ' + str(err_BDF1))
print('error BDF2:  ' + str(err_BDF2))
print('error Theta: ' + str(err_Theta))
print('Error decay BDF1:  '+str(np.divide(err_BDF1[0:4], err_BDF1[1:5])))
print('Error decay BDF2:  '+str(np.divide(err_BDF2[0:4], err_BDF2[1:5])))
print('Error decay Theta: '+str(np.divide(err_Theta[0:4], err_Theta[1:5])))
