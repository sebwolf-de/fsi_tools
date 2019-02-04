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
iter_BDF1 = np.zeros(N)
iter_BDF2 = np.zeros(N)
iter_Theta = np.zeros(N)

norm_u = l2_norm(mass_matrix_u, u_reference)
norm_s = l2_norm(mass_matrix_s, s_reference)

u_x = u_reference[0:ndofs_u]
u_y = u_reference[ndofs_u:2*ndofs_u]
abs_val_u = np.sqrt(u_x**2 + u_y**2)
v = np.max(abs_val_u)
rho = 1.0
nu = 1.0
l = 1.0
re = v*rho*l/nu
print('Reynolds number estimate: ' + str(re))


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
    iter_BDF1[k-1] = np.max(it)
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
    iter_BDF2[k-1] = np.max(it)

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
    iter_Theta[k-1] = np.max(it)


print('norm of u ref: ' + str(norm_u))
print('norm of s ref: ' + str(norm_s))

err_u_BDF1_rel = np.divide(err_u_BDF1, norm_u)
err_u_BDF2_rel = np.divide(err_u_BDF2, norm_u)
err_u_Theta_rel = np.divide(err_u_Theta, norm_u)

decay_u_BDF1 = np.divide(err_u_BDF1[0:N-1], err_u_BDF1[1:N])
decay_u_BDF2 = np.divide(err_u_BDF2[0:N-1], err_u_BDF2[1:N])
decay_u_Theta = np.divide(err_u_Theta[0:N-1], err_u_Theta[1:N])

print('BDF1 Error u:        '+str(err_u_BDF1_rel))
print('Error decay BDF1 u:  '+str(decay_u_BDF1))
print('BDF2 Error u:        '+str(err_u_BDF2_rel))
print('Error decay BDF2 u:  '+str(decay_u_BDF2))
print('Theta Error u:       '+str(err_u_Theta_rel))
print('Error decay Theta u: '+str(decay_u_Theta))

err_s_BDF1_rel = np.divide(err_s_BDF1, norm_s)
err_s_BDF2_rel = np.divide(err_s_BDF2, norm_s)
err_s_Theta_rel = np.divide(err_s_Theta, norm_s)

decay_s_BDF1 = np.divide(err_s_BDF1[0:N-1], err_s_BDF1[1:N])
decay_s_BDF2 = np.divide(err_s_BDF2[0:N-1], err_s_BDF2[1:N])
decay_s_Theta = np.divide(err_s_Theta[0:N-1], err_s_Theta[1:N])

print('BDF1 Error s:        '+str(err_s_BDF1_rel))
print('Error decay BDF1 s:  '+str(decay_s_BDF1))
print('BDF2 Error s:        '+str(err_s_BDF2_rel))
print('Error decay BDF2 s:  '+str(decay_s_BDF2))
print('Theta Error s:       '+str(err_s_Theta_rel))
print('Error decay Theta s: '+str(decay_s_Theta))

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

print('')

print('$\\frac{1}{40}$ &$' + \
    '{:.2e}'.format(err_u_BDF1_rel[1]) + \
    '}$&      &$' + \
    '{:.2e}'.format(err_u_BDF2_rel[1]) + \
    '}$&      &$' + \
    '{:.2e}'.format(err_u_Theta_rel[1]) + \
    '}&$$ \\\\')

print('$\\frac{1}{80}$ &$' + \
    '{:.2e}'.format(err_u_BDF1_rel[2]) + \
    '}$&$' + '{:.2}'.format(np.log2(decay_u_BDF1[1])) +'$&$' + \
    '{:.2e}'.format(err_u_BDF2_rel[2]) + \
    '}$&$' + '{:.3}'.format(np.log2(decay_u_BDF2[1])) +'$&$' + \
    '{:.2e}'.format(err_u_Theta_rel[2]) + \
    '}$&$' + '{:.3}'.format(np.log2(decay_u_Theta[1])) +'$\\\\')

print('$\\frac{1}{160}$&$' + \
    '{:.2e}'.format(err_u_BDF1_rel[3]) + \
    '}$&$' + '{:.2}'.format(np.log2(decay_u_BDF1[2])) +'$&$' + \
    '{:.2e}'.format(err_u_BDF2_rel[3]) + \
    '}$&$' + '{:.3}'.format(np.log2(decay_u_BDF2[2])) +'$&$' + \
    '{:.2e}'.format(err_u_Theta_rel[3]) + \
    '}$&$' + '{:.3}'.format(np.log2(decay_u_Theta[2])) +'$\\\\')

print('$\\frac{1}{320}$&$' + \
    '{:.2e}'.format(err_u_BDF1_rel[4]) + \
    '}$&$' + '{:.2}'.format(np.log2(decay_u_BDF1[3])) +'$&$' + \
    '{:.2e}'.format(err_u_BDF2_rel[4]) + \
    '}$&$' + '{:.3}'.format(np.log2(decay_u_BDF2[3])) +'$&$' + \
    '{:.2e}'.format(err_u_Theta_rel[4]) + \
    '}$&$' + '{:.3}'.format(np.log2(decay_u_Theta[3])) +'$\\\\')

print('')

print('$\\frac{1}{40}$ &$' + \
    '{:.2e}'.format(err_s_BDF1_rel[1]) + \
    '}$&      &$' + \
    '{:.2e}'.format(err_s_BDF2_rel[1]) + \
    '}$&      &$' + \
    '{:.2e}'.format(err_s_Theta_rel[1]) + \
    '}$&$$ \\\\')

print('$\\frac{1}{80}$ &$' + \
    '{:.2e}'.format(err_s_BDF1_rel[2]) + \
    '}$&$' + '{:.2}'.format(np.log2(decay_s_BDF1[1])) +'$&$' + \
    '{:.2e}'.format(err_s_BDF2_rel[2]) + \
    '}$&$' + '{:.3}'.format(np.log2(decay_s_BDF2[1])) +'$&$' + \
    '{:.2e}'.format(err_s_Theta_rel[2]) + \
    '}$&$' + '{:.3}'.format(np.log2(decay_s_Theta[1])) +'$\\\\')

print('$\\frac{1}{160}$&$' + \
    '{:.2e}'.format(err_s_BDF1_rel[3]) + \
    '}$&$' + '{:.2}'.format(np.log2(decay_s_BDF1[2])) +'$&$' + \
    '{:.2e}'.format(err_s_BDF2_rel[3]) + \
    '}$&$' + '{:.3}'.format(np.log2(decay_s_BDF2[2])) +'$&$' + \
    '{:.2e}'.format(err_s_Theta_rel[3]) + \
    '}$&$' + '{:.3}'.format(np.log2(decay_s_Theta[2])) +'$\\\\')

print('$\\frac{1}{320}$&$' + \
    '{:.2e}'.format(err_s_BDF1_rel[4]) + \
    '}$&$' + '{:.2}'.format(np.log2(decay_s_BDF1[3])) +'$&$' + \
    '{:.2e}'.format(err_s_BDF2_rel[4]) + \
    '}$&$' + '{:.3}'.format(np.log2(decay_s_BDF2[3])) +'$&$' + \
    '{:.2e}'.format(err_s_Theta_rel[4]) + \
    '}$&$' + '{:.3}'.format(np.log2(decay_s_Theta[3])) +'$\\\\')

print('BE&$'+str(iter_BDF1[1])+'$&$'+str(iter_BDF1[2])+'$&$'+str(iter_BDF1[3])+'$&$'+str(iter_BDF1[4])+'$\\\\')
print('BDF2&$'+str(iter_BDF2[1])+'$&$'+str(iter_BDF2[2])+'$&$'+str(iter_BDF2[3])+'$&$'+str(iter_BDF2[4])+'$\\\\')
print('Theta&$'+str(iter_Theta[1])+'$&$'+str(iter_Theta[2])+'$&$'+str(iter_Theta[3])+'$&$'+str(iter_Theta[4])+'$\\\\')
