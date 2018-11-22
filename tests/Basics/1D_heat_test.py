#! /bin/env python

import numpy as np

def analytical(t):
    return np.sin(0.5*np.pi*x)*np.cos(t)

def f(t):
    return np.sin(0.5*np.pi*x)*(-np.sin(t) + 0.25*np.pi**2*np.cos(t))

n = 400
dx = 1./n
x = np.arange(0,1+10e-8,dx)
T = 20
Theta = 0.5


K = np.zeros((n+1, n+1))
for k in range(1,n):
    K[k,k-1] = 1/dx**2
    K[k,k] = -2/dx**2
    K[k,k+1] = 1/dx**2

err_BDF1 = np.zeros((5))
err_BDF2 = np.zeros((5))
err_Theta = np.zeros((5))
for t_ind in range(0, 5):

    dt = 2**(-t_ind)
    print('dt = ' + str(dt))

    u_0 = analytical(0)
    u_1 = analytical(dt)

    M_BDF1 = np.eye(n+1) - dt*K
    M_BDF1[0,0] = 1
    M_BDF1[n,n] = 1
    M_BDF2 = 1.5*np.eye(n+1) - dt*K
    M_BDF2[0,0] = 1
    M_BDF2[n,n] = 1
    M_Theta_l = np.eye(n+1) - dt*Theta*K
    M_Theta_l[0,0] = 1
    M_Theta_l[n,n] = 1
    M_Theta_r = np.eye(n+1) + dt*(1-Theta)*K

    u_BDF1 = u_1
    u_BDF2 = u_1
    u_BDF2_old = u_0
    u_Theta = u_1

    N = int(np.round(T/dt))
    for k in range(2,N):
        f_now = f(k*dt)
        f_old = f((k-1)*dt)

        rhs_BDF1 = dt*f_now + u_BDF1
        rhs_BDF1[0] = 0
        rhs_BDF1[n] = np.cos(k*dt)
        u_BDF1 = np.linalg.solve(M_BDF1, rhs_BDF1)

        rhs_BDF2 = dt*f_now + 2*u_BDF2 - 0.5*u_BDF2_old
        rhs_BDF2[0] = 0
        rhs_BDF2[n] = np.cos(k*dt)
        sol = np.linalg.solve(M_BDF2, rhs_BDF2)
        u_BDF2_old = u_BDF2
        u_BDF2 = sol

        rhs_Theta = M_Theta_r.dot(u_Theta) + dt*Theta*f_now + dt*(1-Theta)*f_old
        rhs_Theta[0] = 0
        rhs_Theta[n] = np.cos(k*dt)
        u_Theta = np.linalg.solve(M_Theta_l, rhs_Theta)
        # print 'u_n('+str(k*dt)+') = '+str(u_n)
        # print 'u('+str(k*dt)+') = '+str(analytical(k*dt))
        if k*dt == 18:
            err_BDF1[t_ind] = np.linalg.norm(u_BDF1-analytical(k*dt))
            err_BDF2[t_ind] = np.linalg.norm(u_BDF2-analytical(k*dt))
            err_Theta[t_ind] = np.linalg.norm(u_Theta-analytical(k*dt))
        # print 'error Theta = '+str(np.linalg.norm(u_Theta-analytical(k*dt)))

    print('error BDF1:  ' + str(err_BDF1))
    print('error BDF2:  ' + str(err_BDF2))
    print('error Theta: ' + str(err_Theta))
    print('Error decay BDF1:  '+str(np.divide(err_BDF1[0:4], err_BDF1[1:5])))
    print('Error decay BDF2:  '+str(np.divide(err_BDF2[0:4], err_BDF2[1:5])))
    print('Error decay Theta: '+str(np.divide(err_Theta[0:4], err_Theta[1:5])))
