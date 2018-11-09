#! /bin/env python

import numpy as np

def analytical(t):
    return np.array([np.sin(t), np.cos(t)])

dt = 0.0001
T = 200
Theta = 0.5

u_0 = np.array([0,1])
u_1 = analytical(dt)


M_BDF1 = np.array([[1, -dt], [dt, 1]])
M_BDF2 = np.array([[1.5, -dt], [dt, 1.5]])
M_Theta_l = np.array([[1, -dt*Theta],[dt*Theta, 1]])
M_Theta_r = np.array([[1, dt*(1-Theta)],[-dt*(1-Theta), 1]])

u_BDF1 = np.linalg.solve(M_BDF1, u_0)
u_BDF2 = u_BDF1#np.linalg.solve(M_Theta_l, M_Theta_r.dot(u_0))
u_BDF2_old = u_0
u_Theta = np.linalg.solve(M_Theta_l, M_Theta_r.dot(u_0))

for k in range(2,int(np.round(T/dt+1))):
    u_BDF1 = np.linalg.solve(M_BDF1, u_BDF1)
    sol = np.linalg.solve(M_BDF2, 2*u_BDF2 - 0.5*u_BDF2_old)
    u_BDF2_old = u_BDF2
    u_BDF2 = sol
    u_Theta = np.linalg.solve(M_Theta_l, M_Theta_r.dot(u_Theta))
    #print 'u_n('+str(k*dt)+') = '+str(u_n)
    #print 'u('+str(k*dt)+') = '+str(analytical(k*dt))
print 'error BDF1 =  '+str(np.linalg.norm(u_BDF1-analytical(T)))
print 'error BDF2 =  '+str(np.linalg.norm(u_BDF2-analytical(T)))
print 'error Theta = '+str(np.linalg.norm(u_Theta-analytical(T)))
