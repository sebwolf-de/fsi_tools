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

def fexp(f):
    return int(np.floor(np.log10(abs(f)))) if f != 0 else 0

def fman(f):
    return f/10**fexp(f)

def format_latex(f):
    return "%#.3g" % fman(f) +'\cdot 10^{'+str(fexp(f))+'}'


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
err_u_CN = np.zeros(N)
err_u_TR = np.zeros(N)
err_s_BDF1 = np.zeros(N)
err_s_BDF2 = np.zeros(N)
err_s_CN = np.zeros(N)
err_s_TR = np.zeros(N)
u_BDF1 = np.zeros((2*ndofs_u,N))
u_BDF2 = np.zeros((2*ndofs_u,N))
u_CN = np.zeros((2*ndofs_u,N))
u_TR = np.zeros((2*ndofs_u,N))
s_BDF1 = np.zeros((2*ndofs_s,N))
s_BDF2 = np.zeros((2*ndofs_s,N))
s_CN = np.zeros((2*ndofs_s,N))
s_TR = np.zeros((2*ndofs_s,N))
time_BDF1 = np.zeros(N)
time_BDF2 = np.zeros(N)
time_CN = np.zeros(N)
time_TR = np.zeros(N)
iter_BDF1 = np.zeros(N)
iter_BDF2 = np.zeros(N)
iter_CN = np.zeros(N)
iter_TR = np.zeros(N)

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
    ### BDF1/BE
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

    print('BDF1, max res: ' + str(res))
    iter_BDF1[k-1] = np.max(res)
    time_BDF1[k-1] = step_time

    ### BDF2
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

    print('BDF2, max res: ' + str(res))
    iter_BDF2[k-1] = np.max(res)
    time_BDF2[k-1] = step_time

    ### CN
    input_name = results_dir+'CN_'+str(k)
    f = open(input_name,"rb")
    u_CN[:,k-1] = np.load(f)
    p_CN = np.load(f)
    xs_CN = np.load(f)
    ys_CN= np.load(f)
    f.close()
    s_CN[:,k-1] = np.append(xs_CN, ys_CN)

    # err_u_CN[k-1] = mth.sqrt(l2_norm(mass_matrix_u, u_CN - u_reference)**2
    #                       + l2_norm(stiffness_matrix_u, u_CN - u_reference)**2)
    # err_s_CN[k-1] = mth.sqrt(l2_norm(mass_matrix_s, s_CN - s_reference)**2
    #                       + l2_norm(stiffness_matrix_s, s_CN - s_reference)**2)
    err_u_CN[k-1] = l2_norm(mass_matrix_u, u_CN[:,k-1] - u_reference)
    err_s_CN[k-1] = l2_norm(mass_matrix_s, s_CN[:,k-1] - s_reference)
    # err_u_CN[k-1] = np.linalg.norm(u_CN - u_reference, float('inf'))
    # err_s_CN[k-1] = np.linalg.norm(np.append(xs_CN, ys_CN) - np.append(xs_reference, ys_reference), float('inf'))

    input_name = results_dir+'CN_'+str(k)+'_time'
    f = open(input_name,"rb")
    step_time = np.load(f)
    sol_time = np.load(f)
    res = np.load(f)
    f.close()

    print('CN, max res: ' + str(res))
    iter_CN[k-1] = np.max(res)
    time_CN[k-1] = step_time

    ### TR
    input_name = results_dir+'TR_'+str(k)
    f = open(input_name,"rb")
    u_TR[:,k-1] = np.load(f)
    p_TR = np.load(f)
    xs_TR = np.load(f)
    ys_TR= np.load(f)
    f.close()
    s_TR[:,k-1] = np.append(xs_TR, ys_TR)

    # err_u_TR[k-1] = mth.sqrt(l2_norm(mass_matrix_u, u_TR - u_reference)**2
    #                       + l2_norm(stiffness_matrix_u, u_TR - u_reference)**2)
    # err_s_TR[k-1] = mth.sqrt(l2_norm(mass_matrix_s, s_TR - s_reference)**2
    #                       + l2_norm(stiffness_matrix_s, s_TR - s_reference)**2)
    err_u_TR[k-1] = l2_norm(mass_matrix_u, u_TR[:,k-1] - u_reference)
    err_s_TR[k-1] = l2_norm(mass_matrix_s, s_TR[:,k-1] - s_reference)
    # err_u_TR[k-1] = np.linalg.norm(u_TR - u_reference, float('inf'))
    # err_s_TR[k-1] = np.linalg.norm(np.append(xs_TR, ys_TR) - np.append(xs_reference, ys_reference), float('inf'))

    input_name = results_dir+'TR_'+str(k)+'_time'
    f = open(input_name,"rb")
    step_time = np.load(f)
    sol_time = np.load(f)
    res = np.load(f)
    f.close()

    print('TR, max res: ' + str(res))
    iter_TR[k-1] = np.max(res)
    time_TR[k-1] = step_time

print('norm of u ref: ' + str(norm_u))
print('norm of s ref: ' + str(norm_s))

err_u_BDF1_rel = np.divide(err_u_BDF1, norm_u)
err_u_BDF2_rel = np.divide(err_u_BDF2, norm_u)
err_u_CN_rel = np.divide(err_u_CN, norm_u)
err_u_TR_rel = np.divide(err_u_TR, norm_u)

decay_u_BDF1 = np.divide(err_u_BDF1[0:N-1], err_u_BDF1[1:N])
decay_u_BDF2 = np.divide(err_u_BDF2[0:N-1], err_u_BDF2[1:N])
decay_u_CN = np.divide(err_u_CN[0:N-1], err_u_CN[1:N])
decay_u_TR = np.divide(err_u_TR[0:N-1], err_u_TR[1:N])

print('BDF1 Error u:       '+str(err_u_BDF1_rel))
print('Error decay BDF1 u: '+str(decay_u_BDF1))
print('BDF2 Error u:       '+str(err_u_BDF2_rel))
print('Error decay BDF2 u: '+str(decay_u_BDF2))
print('CN Error u:         '+str(err_u_CN_rel))
print('Error decay CN u:   '+str(decay_u_CN))
print('TR Error u:         '+str(err_u_TR_rel))
print('Error decay TR u:   '+str(decay_u_TR))

err_s_BDF1_rel = np.divide(err_s_BDF1, norm_s)
err_s_BDF2_rel = np.divide(err_s_BDF2, norm_s)
err_s_CN_rel = np.divide(err_s_CN, norm_s)
err_s_TR_rel = np.divide(err_s_TR, norm_s)

decay_s_BDF1 = np.divide(err_s_BDF1[0:N-1], err_s_BDF1[1:N])
decay_s_BDF2 = np.divide(err_s_BDF2[0:N-1], err_s_BDF2[1:N])
decay_s_CN = np.divide(err_s_CN[0:N-1], err_s_CN[1:N])
decay_s_TR = np.divide(err_s_TR[0:N-1], err_s_TR[1:N])

print('BDF1 Error s:       '+str(err_s_BDF1_rel))
print('Error decay BDF1 s: '+str(decay_s_BDF1))
print('BDF2 Error s:       '+str(err_s_BDF2_rel))
print('Error decay BDF2 s: '+str(decay_s_BDF2))
print('CN Error s:         '+str(err_s_CN_rel))
print('Error decay CN s:   '+str(decay_s_CN))
print('TR Error s:         '+str(err_s_TR_rel))
print('Error decay TR s:   '+str(decay_s_TR))

print('average step time BDF1:  '+str(time_BDF1))
print('average step time BDF2:  '+str(time_BDF2))
print('average step time CN:    '+str(time_CN))
print('average step time TR:    '+str(time_TR))

rate_u_BDF1 = np.zeros(N-2)
rate_u_BDF2 = np.zeros(N-2)
rate_u_CN = np.zeros(N-2)
rate_u_TR = np.zeros(N-2)
rate_s_BDF1 = np.zeros(N-2)
rate_s_BDF2 = np.zeros(N-2)
rate_s_CN = np.zeros(N-2)
rate_s_TR = np.zeros(N-2)

for k in range(0,N-2):
    rate_u_BDF1[k] = np.log2(l2_norm(mass_matrix_u, u_BDF1[:,k] - u_BDF1[:,k+1]) / l2_norm(mass_matrix_u, u_BDF1[:,k+1] - u_BDF1[:,k+2]))
    rate_u_BDF2[k] = np.log2(l2_norm(mass_matrix_u, u_BDF2[:,k] - u_BDF2[:,k+1]) / l2_norm(mass_matrix_u, u_BDF2[:,k+1] - u_BDF2[:,k+2]))
    rate_u_CN[k] = np.log2(l2_norm(mass_matrix_u, u_CN[:,k] - u_CN[:,k+1]) / l2_norm(mass_matrix_u, u_CN[:,k+1] - u_CN[:,k+2]))
    rate_u_TR[k] = np.log2(l2_norm(mass_matrix_u, u_TR[:,k] - u_TR[:,k+1]) / l2_norm(mass_matrix_u, u_TR[:,k+1] - u_TR[:,k+2]))

    rate_s_BDF1[k] = np.log2(l2_norm(mass_matrix_s, s_BDF1[:,k] - s_BDF1[:,k+1]) / l2_norm(mass_matrix_s, s_BDF1[:,k+1] - s_BDF1[:,k+2]))
    rate_s_BDF2[k] = np.log2(l2_norm(mass_matrix_s, s_BDF2[:,k] - s_BDF2[:,k+1]) / l2_norm(mass_matrix_s, s_BDF2[:,k+1] - s_BDF2[:,k+2]))
    rate_s_CN[k] = np.log2(l2_norm(mass_matrix_s, s_CN[:,k] - s_CN[:,k+1]) / l2_norm(mass_matrix_s, s_CN[:,k+1] - s_CN[:,k+2]))
    rate_s_TR[k] = np.log2(l2_norm(mass_matrix_s, s_TR[:,k] - s_TR[:,k+1]) / l2_norm(mass_matrix_s, s_TR[:,k+1] - s_TR[:,k+2]))

print('l2 rate BDF1 u: ' + str(rate_u_BDF1))
print('l2 rate BDF2 u: ' + str(rate_u_BDF2))
print('l2 rate CN u:   ' + str(rate_u_CN))
print('l2 rate TR u:   ' + str(rate_u_TR))
print('l2 rate BDF1 s: ' + str(rate_s_BDF1))
print('l2 rate BDF2 s: ' + str(rate_s_BDF2))
print('l2 rate CN s:   ' + str(rate_s_CN))
print('l2 rate TR s:   ' + str(rate_s_TR))

spaces = ['  ', ' ', '']

print('')

print('$0.05$   &$' + \
    format_latex(err_u_BDF1_rel[1])+ \
    '$&       &$' + \
    format_latex(err_u_BDF2_rel[1]) + \
    '$&      &$' + \
    format_latex(err_u_CN_rel[1]) + \
    '$&      &$' + \
    format_latex(err_u_TR_rel[1]) + \
    '$& \\\\')

for k in [1, 2, 3]:
    print('$'+str(0.1*2**(-k-1))+'$'+spaces[k-1]+'&$' + \
        format_latex(err_u_BDF1_rel[k+1]) + \
        '$&$' + "%#.3g" % np.log2(decay_u_BDF1[k]) +'$&$' + \
        format_latex(err_u_BDF2_rel[k+1]) + \
        '$&$' + "%#.3g" % np.log2(decay_u_BDF2[k]) +'$&$' + \
        format_latex(err_u_CN_rel[k+1]) + \
        '$&$' + "%#.3g" % np.log2(decay_u_CN[k])+'$&$' + \
        format_latex(err_u_TR_rel[k+1]) + \
        '$&$' + "%#.3g" % np.log2(decay_u_TR[k]) +'$\\\\')

print('')

print('$0.05$   &$' + \
    format_latex(err_s_BDF1_rel[1])+ \
    '$&       &$' + \
    format_latex(err_s_BDF2_rel[1]) + \
    '$&      &$' + \
    format_latex(err_s_CN_rel[1]) + \
    '$&       &$' + \
    format_latex(err_s_TR_rel[1]) + \
    '$& \\\\')

for k in [1, 2, 3]:
    print('$'+str(0.1*2**(-k-1))+'$'+spaces[k-1]+'&$' + \
        format_latex(err_s_BDF1_rel[k+1]) + \
        '$&$' + "%#.3g" % np.log2(decay_s_BDF1[k]) +'$&$' + \
        format_latex(err_s_BDF2_rel[k+1]) + \
        '$&$' + "%#.3g" % np.log2(decay_s_BDF2[k]) +'$&$' + \
        format_latex(err_s_CN_rel[k+1]) + \
        '$&$' + "%#.3g" % np.log2(decay_s_CN[k])+'$&$' + \
        format_latex(err_s_TR_rel[k+1]) + \
        '$&$' + "%#.3g" % np.log2(decay_s_TR[k]) +'$\\\\')

print('')

print('BE     &$'+format_latex(iter_BDF1[1])+'$&$'+format_latex(iter_BDF1[2])+'$&$'+format_latex(iter_BDF1[3])+'$&$'+format_latex(iter_BDF1[4])+'$\\\\')
print('BDF-$2$&$'+format_latex(iter_BDF2[1])+'$&$'+format_latex(iter_BDF2[2])+'$&$'+format_latex(iter_BDF2[3])+'$&$'+format_latex(iter_BDF2[4])+'$\\\\')
print('CN     &$'+format_latex(iter_CN[1])+'$&$'+format_latex(iter_CN[2])+'$&$'+format_latex(iter_CN[3])+'$&$'+format_latex(iter_CN[4])+'$\\\\')
print('TR     &$'+format_latex(iter_TR[1])+'$&$'+format_latex(iter_TR[2])+'$&$'+format_latex(iter_TR[3])+'$&$'+format_latex(iter_TR[4])+'$\\\\')
