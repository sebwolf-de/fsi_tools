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

#Get mass matrix, all time integrations were being executed on the same mesh
filename = results_dir+'mesh'
print filename
f = file(filename,"rb")
topo_p = np.load(f)
x_p = np.load(f)
y_p = np.load(f)
topo_u = np.load(f)
x_u = np.load(f)
y_u = np.load(f)
c2f = np.load(f)
f.close()

ndofs_u = x_u.shape[0]
mass_1D = assemble.u_v_p1(topo_u,x_u,y_u)
mass_matrix = sparse.vstack([
    sparse.hstack( [mass_1D, sparse.csr_matrix((ndofs_u,ndofs_u)) ] ),
    sparse.hstack( [sparse.csr_matrix((ndofs_u,ndofs_u)), mass_1D] )
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
f.close()

N = 4
err_BDF1 = np.zeros(N)
err_BDF2 = np.zeros(N)
err_Theta = np.zeros(N)

for k in range(1,N+1):
    input_name = results_dir+'BDF1_'+str(k)
    f = file(input_name,"rb")
    u_BDF1 = np.load(f)
    p_BDF1 = np.load(f)
    f.close()

    #err_BDF1[k-1] = mth.sqrt(l2_norm(mass_matrix, u_BDF1 - u_reference)**2
    #                        + l2_norm(stiffness_matrix, u_BDF1 - u_reference)**2)
    err_BDF1[k-1] = l2_norm(mass_matrix, u_BDF1 - u_reference)

    input_name = results_dir+'BDF2_'+str(k)
    f = file(input_name,"rb")
    u_BDF2 = np.load(f)
    p_BDF2 = np.load(f)
    f.close()

    # err_BDF2[k-1] = mth.sqrt(l2_norm(mass_matrix, u_BDF2 - u_reference)**2
    #                       + l2_norm(stiffness_matrix, u_BDF2 - u_reference)**2)
    err_BDF2[k-1] = l2_norm(mass_matrix, u_BDF2 - u_reference)

    input_name = results_dir+'Theta_'+str(k)
    f = file(input_name,"rb")
    u_Theta = np.load(f)
    p_Theta = np.load(f)
    f.close()

    #err_Theta[k-1] = mth.sqrt(l2_norm(mass_matrix, u_Theta - u_reference)**2
    #                       + l2_norm(stiffness_matrix, u_Theta - u_reference)**2)
    err_Theta[k-1] = l2_norm(mass_matrix, u_Theta - u_reference)

print 'BDF1 Error:  '+str(err_BDF1)
print 'BDF2 Error:  '+str(err_BDF2)
print 'Theta Error: '+str(err_Theta)
print 'Error decay BDF1:  '+str(np.divide(err_BDF1[0:5], err_BDF1[1:6]))
print 'Error decay BDF2:  '+str(np.divide(err_BDF2[0:5], err_BDF2[1:6]))
print 'Error decay Theta: '+str(np.divide(err_Theta[0:5], err_Theta[1:6]))
