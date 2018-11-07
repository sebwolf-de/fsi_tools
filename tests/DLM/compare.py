import sys
import numpy as np
from scipy import sparse
import math as mth
import assemble

def l2_norm(M,g):
    l2_g = M.dot(g)
    l2_g = np.dot(l2_g.transpose(),g)
    l2_g = mth.sqrt(l2_g)
    return l2_g


results_dir = 'results/Convergence_Analysis_'+sys.argv[1]+'/'
#results_dir = 'results/BDF_Convergence_Analysis_Cavity/'
print results_dir

#Get mass matrix, all time integrations were being executed on the same mesh
filename = results_dir+'mesh'
f = file(filename,"rb")
topo_p = np.load(f)
x_p = np.load(f)
y_p = np.load(f)
topo_u = np.load(f)
x_u = np.load(f)
y_u = np.load(f)
c2f = np.load(f)
topo_s = np.load(f)
xs_n = np.load(f)
ys_n = np.load(f)
s_lgr = np.load(f)
f.close()

ndofs_u = x_u.shape[0]
mass_1D = assemble.u_v_p1(topo_u,x_u,y_u)
mass_matrix_u = sparse.vstack([
    sparse.hstack( [mass_1D, sparse.csr_matrix((ndofs_u,ndofs_u)) ] ),
    sparse.hstack( [sparse.csr_matrix((ndofs_u,ndofs_u)), mass_1D] )
    ])

ndofs_s = xs_n.shape[0]
mass_1D = assemble.u_v_p1(topo_s,xs_n,ys_n)
mass_matrix_s = sparse.vstack([
    sparse.hstack( [mass_1D, sparse.csr_matrix((ndofs_s,ndofs_s)) ] ),
    sparse.hstack( [sparse.csr_matrix((ndofs_s,ndofs_s)), mass_1D] )
    ])

stiffness_1D = assemble.gradu_gradv_p1(topo_u,x_u,y_u)
stiffness_matrix = sparse.vstack([
    sparse.hstack( [stiffness_1D, sparse.csr_matrix((ndofs_u,ndofs_u)) ] ),
    sparse.hstack( [sparse.csr_matrix((ndofs_u,ndofs_u)), stiffness_1D] )
    ])


input_name = results_dir+'reference'
f = file(input_name,"rb")
u_reference = np.load(f)
p_reference = np.load(f)
xs_reference = np.load(f)
ys_reference = np.load(f)
f.close()

N = 7

err_u_BDF1 = np.zeros(N-1)
err_u_BDF2 = np.zeros(N-1)
err_u_Theta = np.zeros(N-1)
err_s_BDF1 = np.zeros(N-1)
err_s_BDF2 = np.zeros(N-1)
err_s_Theta = np.zeros(N-1)
time_BDF1 = np.zeros(N-1)
time_BDF2 = np.zeros(N-1)
time_Theta = np.zeros(N-1)


for k in range(1,N):
    input_name = results_dir+'BDF1_'+str(k)
    f = file(input_name,"rb")
    u_BDF1 = np.load(f)
    p_BDF1 = np.load(f)
    xs_BDF1 = np.load(f)
    ys_BDF1 = np.load(f)
    f.close()

    #err_u_BDF1[k-1] = mth.sqrt(l2_norm(mass_matrix_u, u_BDF1 - u_reference)**2
    #                        + l2_norm(stiffness_matrix, u_BDF1 - u_reference)**2)
    err_u_BDF1[k-1] = l2_norm(mass_matrix_u, u_BDF1 - u_reference)
    err_s_BDF1[k-1] = l2_norm(mass_matrix_s,
        np.append(xs_BDF1, ys_BDF1) - np.append(xs_reference, ys_reference))

    #input_name = results_dir+'BDF1_'+str(2**k)+'_time'
    #f = file(input_name,"rb")
    #step_time = np.load(f)
    #sol_time = np.load(f)
    #f.close()

    #time_BDF1[k-1] = step_time

    input_name = results_dir+'BDF2_'+str(k)
    f = file(input_name,"rb")
    u_BDF2 = np.load(f)
    p_BDF2 = np.load(f)
    xs_BDF2 = np.load(f)
    ys_BDF2 = np.load(f)
    f.close()

    #err_u_BDF2[k-1] = mth.sqrt(l2_norm(mass_matrix_u, u_BDF2 - u_reference)**2
    #                       + l2_norm(stiffness_matrix, u_BDF2 - u_reference)**2)
    err_u_BDF2[k-1] = l2_norm(mass_matrix_u, u_BDF2 - u_reference)
    err_s_BDF2[k-1] = l2_norm(mass_matrix_s,
        np.append(xs_BDF2, ys_BDF2) - np.append(xs_reference, ys_reference))

    #input_name = results_dir+'BDF2_'+str(2**k)+'_time'
    #f = file(input_name,"rb")
    #step_time = np.load(f)
    #sol_time = np.load(f)
    #f.close()

    #time_BDF2[k-1] = step_time

    input_name = results_dir+'Theta_'+str(k)
    f = file(input_name,"rb")
    u_Theta = np.load(f)
    p_Theta = np.load(f)
    xs_Theta = np.load(f)
    ys_Theta= np.load(f)
    f.close()

    #err_u_Theta[k-1] = mth.sqrt(l2_norm(mass_matrix_u, u_Theta - u_reference)**2
    #                       + l2_norm(stiffness_matrix, u_Theta - u_reference)**2)
    err_u_Theta[k-1] = l2_norm(mass_matrix_u, u_Theta - u_reference)
    err_s_Theta[k-1] = l2_norm(mass_matrix_s,
        np.append(xs_Theta, ys_Theta) - np.append(xs_reference, ys_reference))

    #input_name = results_dir+'Theta_'+str(2**k)+'_time'
    #f = file(input_name,"rb")
    #step_time = np.load(f)
    #sol_time = np.load(f)
    #f.close()

    #time_Theta[k-1] = step_time

print 'BDF1 Error u:        '+str(err_u_BDF1)
print 'Error decay BDF1 u:  '+str(np.divide(err_u_BDF1[0:N-2], err_u_BDF1[1:N-1]))
print 'BDF2 Error u:        '+str(err_u_BDF2)
print 'Error decay BDF2 u:  '+str(np.divide(err_u_BDF2[0:N-2], err_u_BDF2[1:N-1]))
print 'Theta Error u:       '+str(err_u_Theta)
print 'Error decay Theta u: '+str(np.divide(err_u_Theta[0:N-2], err_u_Theta[1:N-1]))


print 'BDF1 Error s:        '+str(err_s_BDF1)
print 'Error decay BDF1 s:  '+str(np.divide(err_s_BDF1[0:N-2], err_s_BDF1[1:N-1]))
print 'BDF2 Error s:        '+str(err_s_BDF2)
print 'Error decay BDF2 s:  '+str(np.divide(err_s_BDF2[0:N-2], err_s_BDF2[1:N-1]))
print 'Theta Error s:       '+str(err_s_Theta)
print 'Error decay Theta s: '+str(np.divide(err_s_Theta[0:N-2], err_s_Theta[1:N-1]))
#print 'Time BDF1: '+str(time_BDF1)
#print 'Time BDF2: '+str(time_BDF2)
#print 'Time Theta: '+str(time_Theta)
