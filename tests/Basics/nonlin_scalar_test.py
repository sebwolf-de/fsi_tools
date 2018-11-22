#! /bin/env python

import numpy as np
import matplotlib.pyplot as plot

def analytical(t):
    return 1. / (1.+t)

dt = 0.05
T = 2.
Theta = 0.5
eps = 10e-8

u_0 = analytical(0)
u_1 = analytical(dt)
K = int(np.round(T/dt+1))

u_BDF1 = np.zeros(K)
u_BDF1[0] = u_0
u_BDF1[1] = u_1

u_BDF2 = np.zeros(K)
u_BDF2[0] = u_0
u_BDF2[1] = u_1

u_Theta = np.zeros(K)
u_Theta[0] = u_0
u_Theta[1] = u_1

u_BDF2_lin = np.zeros(K)
u_BDF2_lin[0] = u_0
u_BDF2_lin[1] = u_1

u_Theta_lin = np.zeros(K)
u_Theta_lin[0] = u_0
u_Theta_lin[1] = u_1

for k in range(2,K):
    v = u_BDF1[k-1]
    #print 't = ' + str(k*dt)
    #print 'start fixpoint iteration:'
    for n in range(0,10):
        v = -dt*v**2 + u_BDF1[k-1]
        res = np.abs(v + dt*v**2 - u_BDF1[k-1])
        #print 'n = ' + str(n) + ', res = ' + str(res)
        if res < eps:
            #print 'fixpoint iteration converged in ' + str(n+1) + ' steps.'
            break
    u_BDF1[k] = v


    v = u_BDF2[k-1]
    #print 't = ' + str(k*dt)
    #print 'start fixpoint iteration:'
    for n in range(0,10):
        v = -2./3.*dt*v**2 + 4./3.*u_BDF2[k-1] - 1./3.*u_BDF2[k-2]
        res = np.abs(v + 2./3.*dt*v**2 - 4./3.*u_BDF2[k-1] + 1./3.*u_BDF2[k-2])
        #print 'n = ' + str(n) + ', res = ' + str(res)
        if res < eps:
            #print 'fixpoint iteration converged in ' + str(n+1) + ' steps.'
            break
    u_BDF2[k] = v

    v = u_Theta[k-1]
    # print 't = ' + str(k*dt)
    # print 'start fixpoint iteration:'
    for n in range(0,10):
        v = -Theta*dt*v**2 - (1-Theta)*dt*u_Theta[k-1]**2 + u_Theta[k-1]
        res = np.abs(v + Theta*dt*v**2 + (1-Theta)*dt*u_Theta[k-1]**2 - u_Theta[k-1])
        # print 'n = ' + str(n) + ', res = ' + str(res)
        if res < eps:
            # print 'fixpoint iteration converged in ' + str(n+1) + ' steps.'
            break
    u_Theta[k] = v

    u_Theta_lin[k] = (u_Theta_lin[k-1] - dt*(1-Theta)*u_Theta_lin[k-1]**2)/(1+dt*Theta*(2*u_Theta_lin[k-1]-u_Theta_lin[k-2]))

    u_BDF2_lin[k] = (4./3.*u_BDF2_lin[k-1] - 1./3.*u_BDF2_lin[k-2]) / (1 + 2.*dt/3.*(2*u_BDF2_lin[k-1]-u_BDF2_lin[k-2]))

print('error BDF1:  ' + str(np.abs(u_BDF1[K-1] - analytical(T))))
print('error BDF2:  ' + str(np.abs(u_BDF2[K-1] - analytical(T))))
print('error Theta: ' + str(np.abs(u_Theta[K-1] - analytical(T))))
print('error BDF2 lin:  ' + str(np.abs(u_BDF2_lin[K-1] - analytical(T))))
print('error Theta lin: ' + str(np.abs(u_Theta_lin[K-1] - analytical(T))))

#plot.plot(np.arange(0,K)*dt, u_BDF1, 'r')
#plot.plot(np.arange(0,K)*dt, u_lin, 'b')
#plot.plot(np.arange(0,K)*dt, analytical(np.arange(0,K)*dt), 'y')
#plot.show()
