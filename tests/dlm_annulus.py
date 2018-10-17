#! /usr/bin/env python

import time
import os
import sys
import numpy as np
from scipy import sparse
import scipy.sparse.linalg as sp_la
import matplotlib.pyplot as plt
import math as mth
import json

# nicola modules
import lin_tri_mesh as lin_t3
import basis_func as shp
import assemble
import la_utils
import viewers
import geom_utils as geom
from shapely.geometry import Polygon

from preconditioner import BlockPreconditioner
from parameters_handler import ParametersHandler

def read_initial_condition(cn_time):
    #filename = './mesh/'+ph.sim_prefix
    #f = file(filename,"rb")
    #topo_p = np.load(f)
    #x_p = np.load(f)
    #y_p = np.load(f)
    #topo_u = np.load(f)
    #x_u = np.load(f)
    #y_u = np.load(f)
    #c2f = np.load(f)
    #topo_s = np.load(f)
    #sx_n = np.load(f)
    #sy_n = np.load(f)
    #s_lgr = np.load(f)
    #t_lgr = np.load(f)
    #f.close()
    filename = "./"+ph.results_directory+"/"
    filename += ph.sim_prefix + '_'+str(cn_time).zfill(4)
    f = file(filename,"rb")
    u = np.load(f)
    p = np.load(f)
    xs = np.load(f)
    ys = np.load(f)
    f.close()
    return xs,ys#,s_lgr,t_lgr,topo_s


def write_mesh():
    filename = results_dir+'mesh'#'./mesh/'+sim_prefix
    f = file(filename,"wb")
    if ph.mesh_prefix != 'thin_':
        np.save(f,topo_p)
        np.save(f,x_p)
        np.save(f,y_p)
        np.save(f,topo_u)
        np.save(f,x_u)
        np.save(f,y_u)
        np.save(f,c2f)
        np.save(f,topo_s)
        np.save(f,sx_n)
        np.save(f,sy_n)
        np.save(f,s_lgr)
        np.save(f,t_lgr)
    elif ph.mesh_prefix == 'thin_':
        np.save(f,topo_p)
        np.save(f,x_p)
        np.save(f,y_p)
        np.save(f,topo_u)
        np.save(f,x_u)
        np.save(f,y_u)
        np.save(f,c2f)
        np.save(f,topo_s)
        np.save(f,sx_n)
        np.save(f,sy_n)
        np.save(f,s_lgr)
    f.close()
    return

def assemble_blockwise_force_BDF1():#ux_n,uy_n,sx_n,sy_n):
    f_rhs_x = 1/ph.dt*Mv11.dot(ux_n)
    f_rhs_y = 1/ph.dt*Mv11.dot(uy_n)

    bc_id = np.where(y_u < delta_x/10)
    f_rhs_y[bc_id] = 0

    bc_id = np.where(y_u > 1-delta_x/10)
    f_rhs_x[bc_id] = 0
    f_rhs_y[bc_id] = 0

    bc_id = np.where(x_u > 1-delta_x/10)
    f_rhs_x[bc_id] = 0
    f_rhs_y[bc_id] = 0

    bc_id = np.where(x_u < delta_x/10)
    f_rhs_x[bc_id] = 0

    s_rhs_x = 1/ph.dt*MX11.dot(dx_n)
    s_rhs_y = 1/ph.dt*MX11.dot(dy_n)

    rhs = np.append(f_rhs_x, f_rhs_y)
    rhs = np.append(rhs, np.zeros((ndofs_p)))
    rhs = np.append(rhs, np.zeros((2*ndofs_s)))
    rhs = np.append(rhs, s_rhs_x)
    rhs = np.append(rhs, s_rhs_y)
    rhs = np.append(rhs, np.zeros(1))

    return rhs

def assemble_blockwise_matrix_BDF1():
    mat1 = sparse.hstack([A_BDF1,
                         -BT,
                         sparse.csr_matrix((ndofs_u*2,ndofs_s*2)),
                         GT,
                         sparse.csr_matrix((ndofs_u*2,1))
                         ])

    mat2 = sparse.hstack([-B,
                          sparse.csr_matrix((ndofs_p,ndofs_p)),
                          sparse.csr_matrix((ndofs_p,ndofs_s*2)),
                          sparse.csr_matrix((ndofs_p,ndofs_s*2)),
                          mean_p.transpose()
                          ])

    mat3 = sparse.hstack([sparse.csr_matrix((ndofs_s*2,ndofs_u*2)),
                          sparse.csr_matrix((ndofs_s*2,ndofs_p)),
                          FX,
                          -MXT,
                          sparse.csr_matrix((ndofs_s*2,1))
                          ])

    mat4 = sparse.hstack([-G,
                          sparse.csr_matrix((ndofs_s*2,ndofs_p)),
                          1/ph.dt*MX,
                          sparse.csr_matrix((ndofs_s*2,ndofs_s*2)),
                          sparse.csr_matrix((ndofs_s*2,1))
                          ])

    mat5 = sparse.hstack([sparse.csr_matrix((1,ndofs_u*2)),
                          mean_p,
                          sparse.csr_matrix((1,ndofs_s*2)),
                          sparse.csr_matrix((1,ndofs_s*2)),
                          sparse.csr_matrix((1,1))
                          ])

    mat = sparse.vstack([mat1,mat2,mat3,mat4,mat5])
    mat = mat.tocsr()
    return mat

def assemble_blockwise_force_BDF2():#ux_n,uy_n,ux_n_old,uy_n_old,sx_n,sy_n,sx_n_old,sy_n_old):
    f_rhs_x = (2*Mv11.dot(ux_n) - 0.5*Mv11.dot(ux_n_old))/ph.dt
    f_rhs_y = (2*Mv11.dot(uy_n) - 0.5*Mv11.dot(uy_n_old))/ph.dt

    bc_id = np.where(y_u < delta_x/10)
    f_rhs_y[bc_id] = 0

    bc_id = np.where(y_u > 1-delta_x/10)
    f_rhs_x[bc_id] = 0
    f_rhs_y[bc_id] = 0

    bc_id = np.where(x_u > 1-delta_x/10)
    f_rhs_x[bc_id] = 0
    f_rhs_y[bc_id] = 0

    bc_id = np.where(x_u < delta_x/10)
    f_rhs_x[bc_id] = 0

    #p_rhs = B.dot(u_n_old)

    l_rhs_x = (2*MX11.dot(dx_n) - 0.5*MX11.dot(dx_n_old))/ph.dt
    l_rhs_y = (2*MX11.dot(dy_n) - 0.5*MX11.dot(dy_n_old))/ph.dt

    #f_rhs_x = np.reshape(f_rhs_x,(ndofs_u))
    #f_rhs_y = np.reshape(f_rhs_y,(ndofs_u))

    rhs = np.append(f_rhs_x, f_rhs_y)
    rhs = np.append(rhs, np.zeros((ndofs_p)))
    #rhs = np.append(rhs, p_rhs)
    rhs = np.append(rhs, np.zeros((2*ndofs_s)))
    rhs = np.append(rhs, l_rhs_x)
    rhs = np.append(rhs, l_rhs_y)
    rhs = np.append(rhs, np.zeros(1))

    return rhs

def assemble_blockwise_matrix_BDF2():
    mat1 = sparse.hstack([A_BDF2,
                          -BT,
                          sparse.csr_matrix((ndofs_u*2,ndofs_s*2)),
                          GT,
                          sparse.csr_matrix((ndofs_u*2,1))
                          ])

    mat2 = sparse.hstack([-B,
                          sparse.csr_matrix((ndofs_p,ndofs_p)),
                          sparse.csr_matrix((ndofs_p,ndofs_s*2)),
                          sparse.csr_matrix((ndofs_p,ndofs_s*2)),
                          mean_p.transpose()
                          ])

    mat3 = sparse.hstack([sparse.csr_matrix((ndofs_s*2,ndofs_u*2)),
                          sparse.csr_matrix((ndofs_s*2,ndofs_p)),
                          FX,
                          -MXT,
                          sparse.csr_matrix((ndofs_s*2,1))
                          ])

    mat4 = sparse.hstack([-G,
                          sparse.csr_matrix((ndofs_s*2,ndofs_p)),
                          1.5/ph.dt*MX,
                          sparse.csr_matrix((ndofs_s*2,ndofs_s*2)),
                          sparse.csr_matrix((ndofs_s*2,1))
                          ])

    mat5 = sparse.hstack([sparse.csr_matrix((1,ndofs_u*2)),
                          mean_p,
                          sparse.csr_matrix((1,ndofs_s*2)),
                          sparse.csr_matrix((1,ndofs_s*2)),
                          sparse.csr_matrix((1,1))
                          ])

    mat = sparse.vstack([mat1,mat2,mat3,mat4,mat5])
    mat = mat.tocsr()
    return mat

def assemble_blockwise_force_Theta():#ux_n,uy_n,p_n,sx_n,sy_n,l_n):
    f_rhs_x = 1/ph.dt*Mv11.dot(ux_n) - 0.5*A11.dot(ux_n) + 0.5*BT1.dot(p_n) - 0.5*GT11.dot(l_n[0:ndofs_s])
    f_rhs_y = 1/ph.dt*Mv11.dot(uy_n) - 0.5*A11.dot(uy_n) + 0.5*BT2.dot(p_n) - 0.5*GT22.dot(l_n[ndofs_s:2*ndofs_s])

    bc_id = np.where(y_u < delta_x/10)
    f_rhs_y[bc_id] = 0

    bc_id = np.where(y_u > 1-delta_x/10)
    f_rhs_x[bc_id] = 0
    f_rhs_y[bc_id] = 0

    bc_id = np.where(x_u > 1-delta_x/10)
    f_rhs_x[bc_id] = 0
    f_rhs_y[bc_id] = 0

    bc_id = np.where(x_u < delta_x/10)
    f_rhs_x[bc_id] = 0

    p_rhs = 0.5*B.dot(u_n)

    s_rhs_x = -0.5*FX11.dot(dx_n) + 0.5*MXT11.dot(np.reshape(l_n[0:ndofs_s],(ndofs_s)))
    s_rhs_y = -0.5*FX22.dot(dy_n) + 0.5*MXT22.dot(np.reshape(l_n[ndofs_s:2*ndofs_s],(ndofs_s)))

    l_rhs_x = -1/ph.dt*MX11.dot(dx_n)
    l_rhs_y = -1/ph.dt*MX11.dot(dy_n)
    l_rhs = np.append(l_rhs_x, l_rhs_y) - np.reshape(0.5*G.dot(u_n), (2*ndofs_s))

    #f_rhs_x = np.reshape(f_rhs_x,(ndofs_u))
    #f_rhs_y = np.reshape(f_rhs_y,(ndofs_u))

    rhs = np.append(f_rhs_x, f_rhs_y)
    rhs = np.append(rhs, p_rhs)
    rhs = np.append(rhs, s_rhs_x)
    rhs = np.append(rhs, s_rhs_y)
    rhs = np.append(rhs, l_rhs)
    rhs = np.append(rhs, np.zeros(1))

    return rhs

def assemble_blockwise_matrix_Theta():
    mat1 = sparse.hstack([A_Theta,
                          -0.5*BT,
                          sparse.csr_matrix((ndofs_u*2,ndofs_s*2)),
                          0.5*GT,
                          sparse.csr_matrix((ndofs_u*2,1))
                          ])

    mat2 = sparse.hstack([0.5*B,
                          sparse.csr_matrix((ndofs_p,ndofs_p)),
                          sparse.csr_matrix((ndofs_p,ndofs_s*2)),
                          sparse.csr_matrix((ndofs_p,ndofs_s*2)),
                          mean_p.transpose()
                          ])

    mat3 = sparse.hstack([sparse.csr_matrix((ndofs_s*2,ndofs_u*2)),
                          sparse.csr_matrix((ndofs_s*2,ndofs_p)),
                          0.5*FX,
                          -0.5*MXT,
                          sparse.csr_matrix((ndofs_s*2,1))
                          ])

    mat4 = sparse.hstack([0.5*G,
                          sparse.csr_matrix((ndofs_s*2,ndofs_p)),
                          -1/ph.dt*MX,
                          sparse.csr_matrix((ndofs_s*2,ndofs_s*2)),
                          sparse.csr_matrix((ndofs_s*2,1))
                          ])

    mat5 = sparse.hstack([sparse.csr_matrix((1,ndofs_u*2)),
                          mean_p,
                          sparse.csr_matrix((1,ndofs_s*2)),
                          sparse.csr_matrix((1,ndofs_s*2)),
                          sparse.csr_matrix((1,1))
                          ])

    mat = sparse.vstack([mat1,mat2,mat3,mat4,mat5])
    mat = mat.tocsr()
    return mat

def unassemble_sol_blocks(sol):
    u_n = sol[0:2*ndofs_u]
    p_n1 = sol[2*ndofs_u:2*ndofs_u+ndofs_p]

    sx_n1 = np.zeros( sx_n.shape )
    sy_n1 = np.zeros( sy_n.shape )

    sx_n1 = sol[2*ndofs_u+ndofs_p:2*ndofs_u+ndofs_p+ndofs_s]

    sy_n1 = sol[2*ndofs_u+ndofs_p+  ndofs_s:
                           2*ndofs_u+ndofs_p+2*ndofs_s]
    return u_n,p_n1,sx_n1,sy_n1

def area_measure(xs,ys):
    area_mes = MX11 * sx_n + MX11 * sy_n
    area_mes = area_mes * np.ones(area_mes.shape)
    area_mes = np.sum(area_mes)
    return area_mes

def l2_norm(M,g):
    l2_g = M.dot(g)
    l2_g = np.dot(l2_g.transpose(),g)
    l2_g = mth.sqrt(l2_g)
    return l2_g

def write_output():
    filename = results_dir +'cn_time_'+str(cn_time).zfill(ph.time_index_digits)
    f = file(filename,"wb")
    np.save(f,u_n)
    np.save(f,p_n)
    np.save(f,sx_n)
    np.save(f,sy_n)
    f.close()
    print '--------------------------------------'
    print 'results saved to:'
    print filename
    print '--------------------------------------'
    return

def write_time():
    filename = results_dir +'time'
    f = file(filename,"wb")
    np.save(f,np.average(step_time))
    np.save(f,np.average(sol_time))
    f.close()

def eval_str_area():
    if ph.mesh_prefix == 'thin_':
        x = np.reshape(sx_n,(sx_n.shape[0],1))
        x = np.vstack([x,[[0]]])
        y = np.reshape(sy_n,(sy_n.shape[0],1))
        y = np.vstack([y,[[0]]])
        area = geom.area_evaluation(x[:,0],y[:,0])
    elif ph.mesh_prefix != 'thin_':
        area = 0
        for row in topo_s:
            x_l = sx_n[row]
            y_l = sy_n[row]
            eval_p = np.zeros((x_l.shape[0],2))
            eval_p[:,0] = x_l
            eval_p[:,1] = y_l
            poly = Polygon(tuple(eval_p.tolist()))
            area+= poly.area
    return area

def get_diffusion():
    return diffusion

def get_energy():
    return energy

def get_prefix():
    return ph.sim_prefix

np.set_printoptions(precision=4)
np.set_printoptions(suppress=True)

if len(sys.argv) > 1:
    ph = ParametersHandler(sys.argv[1])
else:
    ph = ParametersHandler('simulation_parameters_fsi.json')
ph.simulation_info()

nx_p = ph.n_delta_x
delta_x = 1./nx_p
ny_p = nx_p
delta_y = 1./ny_p
(topo_p,x_p,y_p,
 topo_u,x_u,y_u,
 c2f) = lin_t3.mesh_t3_iso_t6(nx_p,ny_p,delta_x,delta_y)

(topo_p,x_p,y_p) = lin_t3.mesh_t3_t0(nx_p,ny_p,delta_x,delta_y)

if ph.mesh_prefix != 'thin_':
    #nx = 4
    #delta_x = 1./nx
    #ny = 2
    #delta_y = 1./ny
    #(topo_s,s_lgr,t_lgr) = lin_t3.mesh_t3(nx,ny,delta_x,delta_y)
    filename = './mesh_collection/' + ph.mesh_name+'.msh'
    (topo_s,s_lgr,t_lgr) = lin_t3.load_msh(filename)

    # =============================================#
    # Square or full ellipsis initial conditions
    #
    #sx_zero = .5 * s_lgr
    #sy_zero = .5 * t_lgr
    #
    #sx_n = .25/.9 * s_lgr;
    #sy_n =.9      * t_lgr;
    # =============================================#
    R0 = .3
    R1 = .5
    ray = R0 + (s_lgr * (R1-R0))
    if ph.equilibrium_at_zero == True:
        sx_zero = np.zeros(s_lgr.shape)
        sy_zero = np.zeros(t_lgr.shape)
        sx_n = 1./1.4*(ray * np.cos(mth.pi/2 * t_lgr))
        sy_n =    1.4*(ray * np.sin(mth.pi/2 * t_lgr))
    elif ph.equilibrium_at_zero == False:
        s_lgr = ray * np.cos(mth.pi/2 * t_lgr)
        t_lgr = ray * np.sin(mth.pi/2 * t_lgr)
        sx_zero = s_lgr
        sy_zero = t_lgr
        sx_n = 1./1.4*(s_lgr)
        sy_n =    1.4*(t_lgr)
elif ph.mesh_prefix == 'thin_':
    ray = .5
    deform = 1.4
    (topo_s,sx_n,sy_n,s_lgr,ieq_s) = lin_t3.lin_str_mesh(ph.n_delta_s,ray,deform)
    sx_zero = np.zeros(s_lgr.shape)
    sy_zero = np.zeros(s_lgr.shape)

#(sx_n,sy_n) = read_initial_condition(40)
ie_s = np.arange(0,s_lgr.shape[0])

if sum(ph.stampa) !=0:
    results_dir = ph.results_directory+'/'+ph.sim_prefix+'/binary_data/'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    ph.dump_to_json(ph.results_directory+'/'+ph.sim_prefix+'/simulation_parameters.json')
    write_mesh()

ndofs_u = max(x_u.shape)
ndofs_p = max(x_p.shape) + topo_p.shape[0]
ndofs_s = max(ie_s)+1

if ph.mesh_prefix != 'thin_':
    MX11 = assemble.u_v_p1_periodic(topo_s,s_lgr,t_lgr,ie_s)
    FX11 = assemble.gradu_gradv_p1_ieq(topo_s,s_lgr,t_lgr,ie_s)
elif ph.mesh_prefix == 'thin_':
    MX11 = assemble.u_v_lin_p1(topo_s,s_lgr,ieq_s)
    FX11 = assemble.gradu_gradv_lin_p1(topo_s,s_lgr,ieq_s)

MX22 = MX11
MXT11 = MX11
MXT22 = MX11

FX11 = ph.kappa * FX11
FX22 = FX11

bc_id = np.where(sy_n < delta_x/10)
FX22 = la_utils.set_diag(FX22,bc_id)
#MX22 = la_utils.set_diag(MX22,bc_id)
MXT22 = la_utils.clear_rows(MXT22,bc_id)


bc_id = np.where(sx_n < delta_x/10)
FX11 = la_utils.set_diag(FX11,bc_id)
#MX11 = la_utils.set_diag(MX11,bc_id)
MXT11 = la_utils.clear_rows(MXT11,bc_id)

MX = sparse.vstack([
    sparse.hstack([MX11,sparse.csr_matrix((ndofs_s,ndofs_s))]),
    sparse.hstack([sparse.csr_matrix((ndofs_s,ndofs_s)),MX22])
    ])

MXT = sparse.vstack([
    sparse.hstack([MXT11,sparse.csr_matrix((ndofs_s,ndofs_s))]),
    sparse.hstack([sparse.csr_matrix((ndofs_s,ndofs_s)),MXT22])
    ])

FX = sparse.vstack([
    sparse.hstack([FX11,sparse.csr_matrix((ndofs_s,ndofs_s))]),
    sparse.hstack([sparse.csr_matrix((ndofs_s,ndofs_s)),FX22])
    ])

rows = np.arange(0,1)
vals = np.ones((1))
mp = sparse.coo_matrix((vals, (rows,rows)), shape=(1,1))

rows = np.arange(0,ndofs_p)
vals = np.ones((ndofs_p))
Mp = sparse.coo_matrix((vals, (rows,rows)), shape=(ndofs_p,ndofs_p))
Mv11 = assemble.u_v_p1(topo_u,x_u,y_u)
A11 = assemble.gradu_gradv_p1(topo_u,x_u,y_u)
A11 = A11/ph.reynolds
A11_BDF1 = A11 + Mv11/ph.dt
A11_BDF2 = A11 + Mv11*1.5/ph.dt
A11_Theta = 0.5 * A11 + Mv11/ph.dt
A22_BDF1 = A11_BDF1
A22_BDF2 = A11_BDF2
A22_Theta = A11_Theta

(BT1,BT2) = assemble.divu_p_p1_iso_p2_p1p0(topo_p,x_p,y_p,
           topo_u,x_u,y_u,c2f)

BT = sparse.vstack([BT1,BT2])
B = BT.transpose()


bc_id = np.where( y_u < delta_x/10)
A22_BDF1 = la_utils.set_diag(A22_BDF1,bc_id)
A22_BDF2 = la_utils.set_diag(A22_BDF2,bc_id)
A22_Theta = la_utils.set_diag(A22_Theta,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)

bc_id = np.where( y_u > 1-delta_x/10)
A11_BDF1 = la_utils.set_diag(A11_BDF1,bc_id)
A22_BDF1 = la_utils.set_diag(A22_BDF1,bc_id)
A11_BDF2 = la_utils.set_diag(A11_BDF2,bc_id)
A22_BDF2 = la_utils.set_diag(A22_BDF2,bc_id)
A11_Theta = la_utils.set_diag(A11_Theta,bc_id)
A22_Theta = la_utils.set_diag(A22_Theta,bc_id)
BT1 = la_utils.clear_rows(BT1,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)

bc_id = np.where( x_u > 1-delta_x/10)
A11_BDF1 = la_utils.set_diag(A11_BDF1,bc_id)
A22_BDF1 = la_utils.set_diag(A22_BDF1,bc_id)
A11_BDF2 = la_utils.set_diag(A11_BDF2,bc_id)
A22_BDF2 = la_utils.set_diag(A22_BDF2,bc_id)
A11_Theta = la_utils.set_diag(A11_Theta,bc_id)
A22_Theta = la_utils.set_diag(A22_Theta,bc_id)
BT1 = la_utils.clear_rows(BT1,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)

bc_id = np.where( x_u < delta_x/10)
A11_BDF1 = la_utils.set_diag(A11_BDF1,bc_id)
A11_BDF2 = la_utils.set_diag(A11_BDF2,bc_id)
A11_Theta = la_utils.set_diag(A11_Theta,bc_id)
BT1 = la_utils.clear_rows(BT1,bc_id)

Mv = sparse.vstack([
    sparse.hstack( [Mv11, sparse.csr_matrix((ndofs_u,ndofs_u)) ] ),
    sparse.hstack( [sparse.csr_matrix((ndofs_u,ndofs_u)), Mv11] )
    ])

A_BDF1 = sparse.vstack([
    sparse.hstack( [A11_BDF1, sparse.csr_matrix((ndofs_u,ndofs_u)) ] ),
    sparse.hstack( [sparse.csr_matrix((ndofs_u,ndofs_u)), A22_BDF1] )
    ])
A_BDF2 = sparse.vstack([
    sparse.hstack( [A11_BDF2, sparse.csr_matrix((ndofs_u,ndofs_u)) ] ),
    sparse.hstack( [sparse.csr_matrix((ndofs_u,ndofs_u)), A22_BDF2] )
    ])
A_Theta = sparse.vstack([
    sparse.hstack( [A11_Theta, sparse.csr_matrix((ndofs_u,ndofs_u)) ] ),
    sparse.hstack( [sparse.csr_matrix((ndofs_u,ndofs_u)), A22_Theta] )
    ])

BT = sparse.vstack([BT1,BT2])

mean_p = np.zeros((1,ndofs_p))
x_l = x_p[topo_p[0,0:3]]
y_l = y_p[topo_p[0,0:3]]
eval_p = np.zeros((0,2))
(phi_dx,phi_dy,phi,omega) = shp.tri_p1(x_l,y_l,eval_p)

for row in topo_p:
    mean_p[0,row] += omega * np.array([1./3.,1./3.,1./3.,1])

ux_n = np.zeros((ndofs_u))
uy_n = np.zeros((ndofs_u))
u_n = np.zeros((2*ndofs_u))
l_n = np.zeros((2*ndofs_s))
p_n = np.zeros((ndofs_p))

ux_n_old = ux_n
uy_n_old = uy_n
u_n_old = u_n
p_n_old = p_n
sx_n_old = sx_n
sy_n_old = sy_n
l_n_old = l_n

dx_n = sx_n - sx_zero
dy_n = sy_n - sy_zero
dx_n_old = dx_n
dy_n_old = dy_n

grade = np.linspace(0.0, 1.0, sum(ph.stampa))

str_area_zero = eval_str_area()

sol_time = np.array([])
step_time = np.array([])

color_id = 0
energy = []
for cn_time in range(0,len(ph.stampa)):
    step_t0 = time.time()
    if ph.mesh_prefix != 'thin_':
        (str_segments,fluid_id) = geom.fluid_intersect_mesh(topo_u,x_u,y_u,
                        topo_s,sx_n,sy_n)
        GT11 = assemble.u_s_p1_thick(x_u,y_u,topo_u,
                        s_lgr,t_lgr,
                        sx_n,sy_n,topo_s,ie_s,
                        str_segments,fluid_id)
    elif ph.mesh_prefix == 'thin_':
        t0 = time.time()
        (str_segments,fluid_id) = geom.fluid_intersect_string(topo_u,x_u,y_u,
                       topo_s,sx_n,sy_n)

        GT11 = assemble.u_s_p1(topo_u,x_u,y_u,
                        topo_s,sx_n,sy_n,s_lgr,ieq_s,
                        str_segments,fluid_id)

        t1 = time.time()
        print 'G time = ' + str((t1-t0))

    GT22 = GT11
    G11 = GT11.transpose()

    G = sparse.vstack([
            sparse.hstack([G11,sparse.csr_matrix((ndofs_s,ndofs_u))]),
            sparse.hstack([sparse.csr_matrix((ndofs_s,ndofs_u)),G11]) ])

    bc_id = np.where( y_u < delta_x/10)
    GT22 = la_utils.clear_rows(GT22,bc_id)

    bc_id = np.where( y_u > 1-delta_x/10)
    GT11 = la_utils.clear_rows(GT11,bc_id)
    GT22 = la_utils.clear_rows(GT22,bc_id)

    bc_id = np.where( x_u > 1-delta_x/10)
    GT11 = la_utils.clear_rows(GT11,bc_id)
    GT22 = la_utils.clear_rows(GT22,bc_id)

    bc_id = np.where( x_u < delta_x/10)
    GT11 = la_utils.clear_rows(GT11,bc_id)

    GT = sparse.vstack([
            sparse.hstack([GT11,sparse.csr_matrix((ndofs_u,ndofs_s))]),
            sparse.hstack([sparse.csr_matrix((ndofs_u,ndofs_s)),GT22]) ])

    if ph.time_integration == 'BDF1':
        mat = assemble_blockwise_matrix_BDF1()
        force = assemble_blockwise_force_BDF1()#ux_n,uy_n,sx_n,sy_n)
    elif ph.time_integration == 'Theta':
        mat = assemble_blockwise_matrix_Theta()
        force = assemble_blockwise_force_Theta()#,uy_n,p_n,sx_n,sy_n,l_n)
    elif ph.time_integration == 'BDF2':
        if cn_time == 0:
            mat = assemble_blockwise_matrix_Theta()
            force = assemble_blockwise_force_Theta()#ux_n,uy_n,sx_n,sy_n)
        else:
            mat = assemble_blockwise_matrix_BDF2()
            force = assemble_blockwise_force_BDF2()#ux_n,uy_n,ux_n_old,uy_n_old,sx_n,sy_n,sx_n_old,sy_n_old)

    sol_t0 = time.time()
    sol = sp_la.spsolve(mat,force)
    sol_t1 = time.time()

    ux_n_old = ux_n
    uy_n_old = uy_n
    u_n_old = u_n
    p_n_old = p_n
    sx_n_old = sx_n
    sy_n_old = sy_n
    dx_n_old = dx_n
    dy_n_old = dy_n
    l_n_old = l_n

    u_n = sol[0:2*ndofs_u]
    ux_n = sol[0:ndofs_u]
    uy_n = sol[ndofs_u:2*ndofs_u]
    p_n = sol[2*ndofs_u:2*ndofs_u+ndofs_p]
    dx_n = sol[2*ndofs_u+ndofs_p:2*ndofs_u+ndofs_p+ndofs_s]
    dy_n = sol[2*ndofs_u+ndofs_p+ndofs_s:2*ndofs_u+ndofs_p+2*ndofs_s]
    sx_n = sx_zero + dx_n
    sy_n = sy_zero + dy_n
    l_n = sol[2*ndofs_u+ndofs_p+2*ndofs_s:2*ndofs_u+ndofs_p+4*ndofs_s]

    str_area = eval_str_area()

    diffusion = str_area/str_area_zero

    p_all_zero = bool(np.all(p_n==0))
    exploded = bool(np.amax(p_n) > 1e+10)

    nrg =(l2_norm(FX,np.append(dx_n, dy_n)))**2 + (l2_norm(Mv,np.append(ux_n, uy_n)))**2
    energy.append(nrg)

    if (exploded==True or p_all_zero == True):
        diffusion = 999
    print '--------------------------------------'
    print 'cn_time   = ' + str(cn_time)
    print 't         = ' + str(cn_time*ph.dt)
    print 'diffusion = ' + str(diffusion)
    print 'energy    = ' + str(nrg)
    print 'pressure == 0? ' + str(p_all_zero)
    print 'exploded     ? ' + str(exploded)
    print '--------------------------------------'

#   if diffusion > 2:
#        break
#    elif diffusion < .8:
#        break

    if ph.stampa[cn_time] == True:
        write_output()
    step_t1 = time.time()
    print 'step time = ' + str((step_t1-step_t0))
    print 'sol  time = ' + str((sol_t1-sol_t0))

    step_time = np.append(step_time, step_t1-step_t0)
    sol_time = np.append(sol_time, sol_t1-sol_t0)

write_time()
