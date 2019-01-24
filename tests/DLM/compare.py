import math as mth
import numpy as np
from scipy import sparse
import sys

sys.path.append('../../modules')
import assemble

def l2_norm(M,g):
    l2_g = M.dot(g)
    l2_g = np.dot(l2_g.transpose(),g)
    l2_g = mth.sqrt(l2_g)
    return l2_g

def l1_norm(M,g):
    l1_g = M.dot(np.abs(g))
    l1_g = np.dot(l1_g.transpose(),np.ones(g.shape))
    return l1_g


results_dir = 'results/Convergence_Analysis_'+sys.argv[1]+'/'
#results_dir = 'results/BDF_Convergence_Analysis_Cavity/'
print(results_dir)

#Get mass matrix, all time integrations were being executed on the same mesh
filename = results_dir+'mesh'
f = open(filename,"rb")
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
t_lgr = np.load(f)
f.close()

ndofs_u = x_u.shape[0]
mass_1D = assemble.u_v_p1(topo_u,x_u,y_u)
mass_matrix_u = sparse.vstack([
    sparse.hstack( [mass_1D, sparse.csr_matrix((ndofs_u,ndofs_u)) ] ),
    sparse.hstack( [sparse.csr_matrix((ndofs_u,ndofs_u)), mass_1D] )
    ])

ie_s = np.arange(0,s_lgr.shape[0])
ndofs_s = max(ie_s)+1
mass_1D = assemble.u_v_p1_periodic(topo_s,s_lgr,t_lgr,ie_s)
mass_matrix_s = sparse.vstack([
    sparse.hstack( [mass_1D, sparse.csr_matrix((ndofs_s,ndofs_s)) ] ),
    sparse.hstack( [sparse.csr_matrix((ndofs_s,ndofs_s)), mass_1D] )
    ])

stiffness_1D = assemble.gradu_gradv_p1(topo_u,x_u,y_u)
stiffness_matrix_u = sparse.vstack([
    sparse.hstack( [stiffness_1D, sparse.csr_matrix((ndofs_u,ndofs_u)) ] ),
    sparse.hstack( [sparse.csr_matrix((ndofs_u,ndofs_u)), stiffness_1D] )
    ])

stiffness_1D = assemble.gradu_gradv_p1_ieq(topo_s,s_lgr,t_lgr,ie_s)
stiffness_matrix_s = sparse.vstack([
    sparse.hstack( [stiffness_1D, sparse.csr_matrix((ndofs_s,ndofs_s)) ] ),
    sparse.hstack( [sparse.csr_matrix((ndofs_s,ndofs_s)), stiffness_1D] )
    ])


input_name = results_dir+'reference'
f = open(input_name,"rb")
u_reference = np.load(f)
p_reference = np.load(f)
xs_reference = np.load(f)
ys_reference = np.load(f)
f.close()
s_reference = np.append(xs_reference, ys_reference)

N = 5

err_u_BDF1 = np.zeros(N)
err_u_BDF2 = np.zeros(N)
err_u_Theta = np.zeros(N)
err_s_BDF1 = np.zeros(N)
err_s_BDF2 = np.zeros(N)
err_s_Theta = np.zeros(N)
u_BDF1 = np.zeros((2*ndofs_u,N))
u_BDF2 = np.zeros((2*ndofs_u,N))
u_Theta = np.zeros((2*ndofs_u,N))
s_BDF1 = np.zeros((2*ndofs_s,N))
s_BDF2 = np.zeros((2*ndofs_s,N))
s_Theta = np.zeros((2*ndofs_s,N))
time_BDF1 = np.zeros(N)
time_BDF2 = np.zeros(N)
time_Theta = np.zeros(N)


for k in range(1,N+1):
    input_name = results_dir+'BDF1_'+str(k)
    f = open(input_name,"rb")
    u_BDF1[:,k-1] = np.load(f)
    p_BDF1 = np.load(f)
    xs_BDF1 = np.load(f)
    ys_BDF1 = np.load(f)
    f.close()
    s_BDF1[:,k-1] = np.append(xs_BDF1, ys_BDF1)

    # err_u_BDF1[k-1] = mth.sqrt(l2_norm(mass_matrix_u, u_BDF1 - u_reference)**2
    #                        + l2_norm(stiffness_matrix_u, u_BDF1 - u_reference)**2)
    # err_s_BDF1[k-1] = mth.sqrt(l2_norm(mass_matrix_s, s_BDF1 - s_reference)**2
    #                        + l2_norm(stiffness_matrix_s, s_BDF1 - s_reference)**2)
    err_u_BDF1[k-1] = l2_norm(mass_matrix_u, u_BDF1[:,k-1] - u_reference)
    err_s_BDF1[k-1] = l2_norm(mass_matrix_s, s_BDF1[:,k-1] - s_reference)
    # err_u_BDF1[k-1] = np.linalg.norm(u_BDF1 - u_reference, float('inf'))
    # err_s_BDF1[k-1] = np.linalg.norm(np.append(xs_BDF1, ys_BDF1) - np.append(xs_reference, ys_reference), float('inf'))

    input_name = results_dir+'BDF1_'+str(k)+'_time'
    f = open(input_name,"rb")
    step_time = np.load(f)
    sol_time = np.load(f)
    res = np.load(f)
    f.close()

    print('BDF1, max_iter, mean_iter, min_iter:')
    it = np.count_nonzero(res > 0, 1)
    print('     ', np.max(it), np.mean(it), np.min(it))
    time_BDF1[k-1] = step_time

    input_name = results_dir+'BDF2_'+str(k)
    f = open(input_name,"rb")
    u_BDF2[:,k-1] = np.load(f)
    p_BDF2 = np.load(f)
    xs_BDF2 = np.load(f)
    ys_BDF2 = np.load(f)
    f.close()
    s_BDF2[:,k-1] = np.append(xs_BDF2, ys_BDF2)

    # err_u_BDF2[k-1] = mth.sqrt(l2_norm(mass_matrix_u, u_BDF2 - u_reference)**2
    #                       + l2_norm(stiffness_matrix_u, u_BDF2 - u_reference)**2)
    # err_s_BDF2[k-1] = mth.sqrt(l2_norm(mass_matrix_s, s_BDF2 - s_reference)**2
    #                       + l2_norm(stiffness_matrix_s, s_BDF2 - s_reference)**2)
    err_u_BDF2[k-1] = l2_norm(mass_matrix_u, u_BDF2[:,k-1] - u_reference)
    err_s_BDF2[k-1] = l2_norm(mass_matrix_s, s_BDF2[:,k-1] - s_reference)
    # err_u_BDF2[k-1] = np.linalg.norm(u_BDF2 - u_reference, float('inf'))
    # err_s_BDF2[k-1] = np.linalg.norm(np.append(xs_BDF2, ys_BDF2) - np.append(xs_reference, ys_reference), float('inf'))

    input_name = results_dir+'BDF2_'+str(k)+'_time'
    f = open(input_name,"rb")
    step_time = np.load(f)
    sol_time = np.load(f)
    res = np.load(f)
    f.close()

    print('BDF2, max_iter, mean_iter, min_iter:')
    it = np.count_nonzero(res > 0, 1)
    print('     ', np.max(it), np.mean(it), np.min(it))
    time_BDF2[k-1] = step_time

    input_name = results_dir+'Theta_'+str(k)
    f = open(input_name,"rb")
    u_Theta[:,k-1] = np.load(f)
    p_Theta = np.load(f)
    xs_Theta = np.load(f)
    ys_Theta= np.load(f)
    f.close()
    s_Theta[:,k-1] = np.append(xs_Theta, ys_Theta)

    # err_u_Theta[k-1] = mth.sqrt(l2_norm(mass_matrix_u, u_Theta - u_reference)**2
    #                       + l2_norm(stiffness_matrix_u, u_Theta - u_reference)**2)
    # err_s_Theta[k-1] = mth.sqrt(l2_norm(mass_matrix_s, s_Theta - s_reference)**2
    #                       + l2_norm(stiffness_matrix_s, s_Theta - s_reference)**2)
    err_u_Theta[k-1] = l2_norm(mass_matrix_u, u_Theta[:,k-1] - u_reference)
    err_s_Theta[k-1] = l2_norm(mass_matrix_s, s_Theta[:,k-1] - s_reference)
    # err_u_Theta[k-1] = np.linalg.norm(u_Theta - u_reference, float('inf'))
    # err_s_Theta[k-1] = np.linalg.norm(np.append(xs_Theta, ys_Theta) - np.append(xs_reference, ys_reference), float('inf'))

    input_name = results_dir+'Theta_'+str(k)+'_time'
    f = open(input_name,"rb")
    step_time = np.load(f)
    sol_time = np.load(f)
    res = np.load(f)
    f.close()

    print('Theta, max_iter, mean_iter, min_iter:')
    it = np.count_nonzero(res > 0, 1)
    print('      ', np.max(it), np.mean(it), np.min(it))

    time_Theta[k-1] = step_time

print('BDF1 Error u:        '+str(err_u_BDF1))
print('Error decay BDF1 u:  '+str(np.divide(err_u_BDF1[0:N-1], err_u_BDF1[1:N])))
print('BDF2 Error u:        '+str(err_u_BDF2))
print('Error decay BDF2 u:  '+str(np.divide(err_u_BDF2[0:N-1], err_u_BDF2[1:N])))
print('Theta Error u:       '+str(err_u_Theta))
print('Error decay Theta u: '+str(np.divide(err_u_Theta[0:N-1], err_u_Theta[1:N])))


print('BDF1 Error s:        '+str(err_s_BDF1))
print('Error decay BDF1 s:  '+str(np.divide(err_s_BDF1[0:N-1], err_s_BDF1[1:N])))
print('BDF2 Error s:        '+str(err_s_BDF2))
print('Error decay BDF2 s:  '+str(np.divide(err_s_BDF2[0:N-1], err_s_BDF2[1:N])))
print('Theta Error s:       '+str(err_s_Theta))
print('Error decay Theta s: '+str(np.divide(err_s_Theta[0:N-1], err_s_Theta[1:N])))

print('average step time BDF1:  '+str(time_BDF1))
print('average step time BDF2:  '+str(time_BDF2))
print('average step time Theta: '+str(time_Theta))

rate_u_BDF1 = np.zeros(N-2)
rate_u_BDF2 = np.zeros(N-2)
rate_u_Theta = np.zeros(N-2)
rate_s_BDF1 = np.zeros(N-2)
rate_s_BDF2 = np.zeros(N-2)
rate_s_Theta = np.zeros(N-2)

# for k in range(0,N-2):
#     rate_u_BDF1[k] = np.log2(l1_norm(mass_matrix_u, u_BDF1[:,k] - u_BDF1[:,k+1]) / l1_norm(mass_matrix_u, u_BDF1[:,k+1] - u_BDF1[:,k+2]))
#     rate_u_BDF2[k] = np.log2(l1_norm(mass_matrix_u, u_BDF2[:,k] - u_BDF2[:,k+1]) / l1_norm(mass_matrix_u, u_BDF2[:,k+1] - u_BDF2[:,k+2]))
#     rate_u_Theta[k] = np.log2(l1_norm(mass_matrix_u, u_Theta[:,k] - u_Theta[:,k+1]) / l1_norm(mass_matrix_u, u_Theta[:,k+1] - u_Theta[:,k+2]))
#
#     rate_s_BDF1[k] = np.log2(l1_norm(mass_matrix_s, s_BDF1[:,k] - s_BDF1[:,k+1]) / l1_norm(mass_matrix_s, s_BDF1[:,k+1] - s_BDF1[:,k+2]))
#     rate_s_BDF2[k] = np.log2(l1_norm(mass_matrix_s, s_BDF2[:,k] - s_BDF2[:,k+1]) / l1_norm(mass_matrix_s, s_BDF2[:,k+1] - s_BDF2[:,k+2]))
#     rate_s_Theta[k] = np.log2(l1_norm(mass_matrix_s, s_Theta[:,k] - s_Theta[:,k+1]) / l1_norm(mass_matrix_s, s_Theta[:,k+1] - s_Theta[:,k+2]))
#
# print('l1 rate BDF1 u: ' + str(rate_u_BDF1))
# print('l1 rate BDF2 u: ' + str(rate_u_BDF2))
# print('l1 rate Theta u: ' + str(rate_u_Theta))
# print('l1 rate BDF1 s: ' + str(rate_s_BDF1))
# print('l1 rate BDF2 s: ' + str(rate_s_BDF2))
# print('l1 rate Theta s: ' + str(rate_s_Theta))


for k in range(0,N-2):
    rate_u_BDF1[k] = np.log2(l2_norm(mass_matrix_u, u_BDF1[:,k] - u_BDF1[:,k+1]) / l2_norm(mass_matrix_u, u_BDF1[:,k+1] - u_BDF1[:,k+2]))
    rate_u_BDF2[k] = np.log2(l2_norm(mass_matrix_u, u_BDF2[:,k] - u_BDF2[:,k+1]) / l2_norm(mass_matrix_u, u_BDF2[:,k+1] - u_BDF2[:,k+2]))
    rate_u_Theta[k] = np.log2(l2_norm(mass_matrix_u, u_Theta[:,k] - u_Theta[:,k+1]) / l2_norm(mass_matrix_u, u_Theta[:,k+1] - u_Theta[:,k+2]))

    rate_s_BDF1[k] = np.log2(l2_norm(mass_matrix_s, s_BDF1[:,k] - s_BDF1[:,k+1]) / l2_norm(mass_matrix_s, s_BDF1[:,k+1] - s_BDF1[:,k+2]))
    rate_s_BDF2[k] = np.log2(l2_norm(mass_matrix_s, s_BDF2[:,k] - s_BDF2[:,k+1]) / l2_norm(mass_matrix_s, s_BDF2[:,k+1] - s_BDF2[:,k+2]))
    rate_s_Theta[k] = np.log2(l2_norm(mass_matrix_s, s_Theta[:,k] - s_Theta[:,k+1]) / l2_norm(mass_matrix_s, s_Theta[:,k+1] - s_Theta[:,k+2]))

print('l2 rate BDF1 u: ' + str(rate_u_BDF1))
print('l2 rate BDF2 u: ' + str(rate_u_BDF2))
print('l2 rate Theta u: ' + str(rate_u_Theta))
print('l2 rate BDF1 s: ' + str(rate_s_BDF1))
print('l2 rate BDF2 s: ' + str(rate_s_BDF2))
print('l2 rate Theta s: ' + str(rate_s_Theta))

# for k in range(0,N-2):
#     rate_u_BDF1[k] = np.log2(np.linalg.norm(u_BDF1[:,k] - u_BDF1[:,k+1], float('inf')) / np.linalg.norm(u_BDF1[:,k+1] - u_BDF1[:,k+2], float('inf')))
#     rate_u_BDF2[k] = np.log2(np.linalg.norm(u_BDF2[:,k] - u_BDF2[:,k+1], float('inf')) / np.linalg.norm(u_BDF2[:,k+1] - u_BDF2[:,k+2], float('inf')))
#     rate_u_Theta[k] = np.log2(np.linalg.norm(u_Theta[:,k] - u_Theta[:,k+1], float('inf')) / np.linalg.norm(u_Theta[:,k+1] - u_Theta[:,k+2], float('inf')))
#
#     rate_s_BDF1[k] = np.log2(np.linalg.norm(s_BDF1[:,k] - s_BDF1[:,k+1], float('inf')) / np.linalg.norm(s_BDF1[:,k+1] - s_BDF1[:,k+2], float('inf')))
#     rate_s_BDF2[k] = np.log2(np.linalg.norm(s_BDF2[:,k] - s_BDF2[:,k+1], float('inf')) / np.linalg.norm(s_BDF2[:,k+1] - s_BDF2[:,k+2], float('inf')))
#     rate_s_Theta[k] = np.log2(np.linalg.norm(s_Theta[:,k] - s_Theta[:,k+1], float('inf')) / np.linalg.norm(s_Theta[:,k+1] - s_Theta[:,k+2], float('inf')))
#
# print('l inf rate BDF1 u: ' + str(rate_u_BDF1))
# print('l inf rate BDF2 u: ' + str(rate_u_BDF2))
# print('l inf rate Theta u: ' + str(rate_u_Theta))
# print('l inf rate BDF1 s: ' + str(rate_s_BDF1))
# print('l inf rate BDF2 s: ' + str(rate_s_BDF2))
# print('l inf rate Theta s: ' + str(rate_s_Theta))
