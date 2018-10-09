#! /usr/bin/env python

import time
import os
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
    #xs_n = np.load(f)
    #ys_n = np.load(f)
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
        np.save(f,xs_n)
        np.save(f,ys_n)
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
        np.save(f,xs_n)
        np.save(f,ys_n)
        np.save(f,s_lgr)
    f.close()
    return

def assemble_blockwise_force_BDF1(ux_n,uy_n,xs_n,ys_n):

    size = 2*ndofs_u+ndofs_p+1+4*ndofs_s
    rhs = np.zeros((size))

    f_rhs_x = 1/ph.dt*Mv11.dot(ux_n)
    f_rhs_y = 1/ph.dt*Mv11.dot(uy_n)

    bc_id = np.where( y_u < delta_x/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    bc_id = np.where( y_u > 1-delta_x/10)
    f_rhs_x[bc_id,:] = 1.
    f_rhs_y[bc_id,:] = 0

    bc_id = np.where( x_u > 1-delta_x/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    bc_id = np.where( x_u < delta_x/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    s_rhs_x = 1/ph.dt*MX11.dot(dx_n)
    s_rhs_y = 1/ph.dt*MX11.dot(dy_n)

    f_rhs_x = np.reshape(f_rhs_x,(ndofs_u))
    f_rhs_y = np.reshape(f_rhs_y,(ndofs_u))

    s = 0
    e = ndofs_u
    rhs[s:e] = f_rhs_x
    s = ndofs_u
    e = 2*ndofs_u
    rhs[s:e] = f_rhs_y

    s_rhs_x = np.reshape(s_rhs_x,(ndofs_s))
    s_rhs_y = np.reshape(s_rhs_y,(ndofs_s))

    s = 2*ndofs_u+ndofs_p+2*ndofs_s
    e = 2*ndofs_u+ndofs_p+3*ndofs_s
    rhs[s:e] = s_rhs_x

    s = 2*ndofs_u+ndofs_p+3*ndofs_s
    e = 2*ndofs_u+ndofs_p+4*ndofs_s
    rhs[s:e] = s_rhs_y

    return rhs

def assemble_blockwise_matrix_BDF1():
    mat1 = sparse.hstack([A_BDF1,-BT])
    mat1 = sparse.hstack([mat1,sparse.csr_matrix((ndofs_u*2,ndofs_s*2))])
    mat1 = sparse.hstack([mat1,GT])
    mat1 = sparse.hstack([mat1,sparse.csr_matrix((ndofs_u*2,1))])

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


def assemble_blockwise_force_BDF2(ux_n,uy_n,ux_n_old,uy_n_old,xs_n,ys_n,xs_n_old,ys_n_old):

    size = 2*ndofs_u+ndofs_p+1+4*ndofs_s
    rhs = np.zeros((size))

    f_rhs_x = (2*Mv11.dot(ux_n) - 0.5*Mv11.dot(ux_n_old))/ph.dt
    f_rhs_y = (2*Mv11.dot(uy_n) - 0.5*Mv11.dot(uy_n_old))/ph.dt

    bc_id = np.where( y_u < delta_x/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    bc_id = np.where( y_u > 1-delta_x/10)
    f_rhs_x[bc_id,:] = 1.
    f_rhs_y[bc_id,:] = 0

    bc_id = np.where( x_u > 1-delta_x/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    bc_id = np.where( x_u < delta_x/10)
    f_rhs_x[bc_id,:] = 0
    f_rhs_y[bc_id,:] = 0

    s_rhs_x = (2*MX11.dot(xs_n - xs_zero) - 0.5*MX11.dot(xs_n_old - xs_zero))/ph.dt
    s_rhs_y = (2*MX11.dot(ys_n - ys_zero) - 0.5*MX11.dot(ys_n_old - ys_zero))/ph.dt

    f_rhs_x = np.reshape(f_rhs_x,(ndofs_u))
    f_rhs_y = np.reshape(f_rhs_y,(ndofs_u))

    s = 0
    e = ndofs_u
    rhs[s:e] = f_rhs_x
    s = ndofs_u
    e = 2*ndofs_u
    rhs[s:e] = f_rhs_y

    s_rhs_x = np.reshape(s_rhs_x,(ndofs_s))
    s_rhs_y = np.reshape(s_rhs_y,(ndofs_s))

    s = 2*ndofs_u+ndofs_p+2*ndofs_s
    e = 2*ndofs_u+ndofs_p+3*ndofs_s
    rhs[s:e] = s_rhs_x

    s = 2*ndofs_u+ndofs_p+3*ndofs_s
    e = 2*ndofs_u+ndofs_p+4*ndofs_s
    rhs[s:e] = s_rhs_y

    return rhs

def assemble_blockwise_matrix_BDF2():
    mat1 = sparse.hstack([A_BDF2,-BT])
    mat1 = sparse.hstack([mat1,sparse.csr_matrix((ndofs_u*2,ndofs_s*2))])
    mat1 = sparse.hstack([mat1,GT])
    mat1 = sparse.hstack([mat1,sparse.csr_matrix((ndofs_u*2,1))])

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

def unassemble_sol_blocks(sol):
    u_n1 = sol[0:2*ndofs_u]
    p_n1 = sol[2*ndofs_u:2*ndofs_u+ndofs_p]

    xs_n1 = np.zeros( xs_n.shape )
    ys_n1 = np.zeros( ys_n.shape )

    xs_n1 = sol[2*ndofs_u+ndofs_p:2*ndofs_u+ndofs_p+ndofs_s]

    ys_n1 = sol[2*ndofs_u+ndofs_p+  ndofs_s:
                           2*ndofs_u+ndofs_p+2*ndofs_s]
    return u_n1,p_n1,xs_n1,ys_n1

def area_measure(xs,ys):
    area_mes = MX11 * xs_n + MX11 * ys_n
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
    np.save(f,u_n1)
    np.save(f,p)
    np.save(f,xs_n1)
    np.save(f,ys_n1)
    f.close()
    print '--------------------------------------'
    print 'results saved to:'
    print filename
    print '--------------------------------------'
    return

def eval_str_area():
    if ph.mesh_prefix == 'thin_':
        x = np.reshape(xs_n,(xs_n.shape[0],1))
        x = np.vstack([x,[[0]]])
        y = np.reshape(ys_n,(ys_n.shape[0],1))
        y = np.vstack([y,[[0]]])
        area = geom.area_evaluation(x[:,0],y[:,0])
    elif ph.mesh_prefix != 'thin_':
        area = 0
        for row in topo_s:
            x_l = xs_n[row]
            y_l = ys_n[row]
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

    viewers.tri_plot(s_lgr, t_lgr, topo_s)
    # =============================================#
    # Square or full ellipsis initial conditions
    #
    #xs_zero = .5 * s_lgr
    #ys_zero = .5 * t_lgr
    #
    #xs_n = .25/.9 * s_lgr;
    #ys_n =.9      * t_lgr;
    # =============================================#
    R0 = .3
    R1 = .5
    ray = R0 + (s_lgr * (R1-R0))
    if ph.equilibrium_at_zero == True:
        xs_zero = np.zeros(s_lgr.shape)
        ys_zero = np.zeros(t_lgr.shape)
        xs_n = 1./1.4*(ray * np.cos(mth.pi/2 * t_lgr))
        ys_n =    1.4*(ray * np.sin(mth.pi/2 * t_lgr))
    elif ph.equilibrium_at_zero == False:
        s_lgr = 0.25 + 0.5*s_lgr #ray * np.cos(mth.pi/2 * t_lgr)
        t_lgr = 0.25 + 0.5*t_lgr #ray * np.sin(mth.pi/2 * t_lgr)
        xs_zero = s_lgr
        ys_zero = t_lgr
        xs_n = (s_lgr)
        ys_n = (t_lgr)
elif ph.mesh_prefix == 'thin_':
    ray = .5
    deform = 1.4
    (topo_s,xs_n,ys_n,s_lgr,ieq_s) = lin_t3.lin_str_mesh(ph.n_delta_s,ray,deform)
    xs_zero = np.zeros(s_lgr.shape)
    ys_zero = np.zeros(s_lgr.shape)

#xs_zero = np.zeros(s_lgr.shape)
#ys_zero = np.zeros(t_lgr.shape)
#xs_n = np.zeros(s_lgr.shape)
#ys_n = np.zeros(s_lgr.shape)

#(xs_n,ys_n) = read_initial_condition(40)
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

bc_id = np.where( ys_n < delta_x/10)
FX22 = la_utils.set_diag(FX22,bc_id)
#MX22 = la_utils.set_diag(MX22,bc_id)
MXT22 = la_utils.clear_rows(MXT22,bc_id)


bc_id = np.where( xs_n < delta_x/10)
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
A22_BDF1 = A11_BDF1
A22_BDF2 = A11_BDF2

(BT1,BT2) = assemble.divu_p_p1_iso_p2_p1p0(topo_p,x_p,y_p,
           topo_u,x_u,y_u,c2f)

BT = sparse.vstack([BT1,BT2])
B = BT.transpose()


bc_id = np.where( y_u < delta_x/10)
A11_BDF1 = la_utils.set_diag(A11_BDF1,bc_id)
A22_BDF1 = la_utils.set_diag(A22_BDF1,bc_id)
A11_BDF2 = la_utils.set_diag(A11_BDF2,bc_id)
A22_BDF2 = la_utils.set_diag(A22_BDF2,bc_id)
BT1 = la_utils.clear_rows(BT1,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)


bc_id = np.where( y_u > 1-delta_x/10)
A11_BDF1 = la_utils.set_diag(A11_BDF1,bc_id)
A22_BDF1 = la_utils.set_diag(A22_BDF1,bc_id)
A11_BDF2 = la_utils.set_diag(A11_BDF2,bc_id)
A22_BDF2 = la_utils.set_diag(A22_BDF2,bc_id)
BT1 = la_utils.clear_rows(BT1,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)

bc_id = np.where( x_u > 1-delta_x/10)
A11_BDF1 = la_utils.set_diag(A11_BDF1,bc_id)
A22_BDF1 = la_utils.set_diag(A22_BDF1,bc_id)
A11_BDF2 = la_utils.set_diag(A11_BDF2,bc_id)
A22_BDF2 = la_utils.set_diag(A22_BDF2,bc_id)
BT1 = la_utils.clear_rows(BT1,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)

bc_id = np.where( x_u < delta_x/10)
A11_BDF1 = la_utils.set_diag(A11_BDF1,bc_id)
A22_BDF1 = la_utils.set_diag(A22_BDF1,bc_id)
A11_BDF2 = la_utils.set_diag(A11_BDF2,bc_id)
A22_BDF2 = la_utils.set_diag(A22_BDF2,bc_id)
BT1 = la_utils.clear_rows(BT1,bc_id)
BT2 = la_utils.clear_rows(BT2,bc_id)

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

BT = sparse.vstack([BT1,BT2])

mean_p = np.zeros((1,ndofs_p))
x_l = x_p[topo_p[0,0:3]]
y_l = y_p[topo_p[0,0:3]]
eval_p = np.zeros((0,2))
(phi_dx,phi_dy,phi,omega) = shp.tri_p1(x_l,y_l,eval_p)

for row in topo_p:
    mean_p[0,row] += omega * np.array([1./3.,1./3.,1./3.,1])

#filename = results_dir +'cn_time_'+str(0).zfill(ph.time_index_digits)
#f = file(filename,"rb")
#u_n_old = np.load(f)
#np.load(f)
#xs_n_old = np.load(f)
#ys_n_old = np.load(f)
#f.close()

ux_n = np.zeros((ndofs_u,1))
uy_n = np.zeros((ndofs_u,1))
#ux_n_old = u_n_old[0      :  ndofs_u]
#uy_n_old = u_n_old[ndofs_u:2*ndofs_u]
#ux_n_old = np.reshape(ux_n_old, ux_n.shape)
#uy_n_old = np.reshape(uy_n_old, uy_n.shape)

ux_n_old = ux_n
uy_n_old = uy_n
xs_n_old = xs_n
ys_n_old = ys_n

dx_n = xs_n - xs_zero
dy_n = ys_n - ys_zero

grade = np.linspace(0.0, 1.0, sum(ph.stampa))

str_area_zero = eval_str_area()

color_id = 0
energy = []
for cn_time in range(0,len(ph.stampa)):
    step_t0 = time.time()
    if ph.mesh_prefix != 'thin_':
        (str_segments,fluid_id) = geom.fluid_intersect_mesh(topo_u,x_u,y_u,
                        topo_s,xs_n,ys_n)
        GT11 = assemble.u_s_p1_thick(x_u,y_u,topo_u,
                        s_lgr,t_lgr,
                        xs_n,ys_n,topo_s,ie_s,
                        str_segments,fluid_id)
    elif ph.mesh_prefix == 'thin_':
        t0 = time.time()
        (str_segments,fluid_id) = geom.fluid_intersect_string(topo_u,x_u,y_u,
                       topo_s,xs_n,ys_n)

        GT11 = assemble.u_s_p1(topo_u,x_u,y_u,
                        topo_s,xs_n,ys_n,s_lgr,ieq_s,
                        str_segments,fluid_id)

        t1 = time.time()
        print 'G time = ' + str((t1-t0))

    GT22 = GT11
    G11 = GT11.transpose()

    G = sparse.vstack([
            sparse.hstack([G11,sparse.csr_matrix((ndofs_s,ndofs_u))]),
            sparse.hstack([sparse.csr_matrix((ndofs_s,ndofs_u)),G11]) ])

    #bc_id = np.where( y_u < delta_x/10)
    #GT11 = la_utils.clear_rows(GT11,bc_id)
    #GT22 = la_utils.clear_rows(GT22,bc_id)

    #bc_id = np.where( y_u > 1-delta_x/10)
    #GT11 = la_utils.clear_rows(GT11,bc_id)
    #GT22 = la_utils.clear_rows(GT22,bc_id)

    #bc_id = np.where( x_u > 1-delta_x/10)
    #GT11 = la_utils.clear_rows(GT11,bc_id)
    #GT22 = la_utils.clear_rows(GT22,bc_id)

    #bc_id = np.where( x_u < delta_x/10)
    #GT11 = la_utils.clear_rows(GT11,bc_id)
    #GT22 = la_utils.clear_rows(GT22,bc_id)

    GT = sparse.vstack([
            sparse.hstack([GT11,sparse.csr_matrix((ndofs_u,ndofs_s))]),
            sparse.hstack([sparse.csr_matrix((ndofs_u,ndofs_s)),GT22]) ])

    if ph.time_integration == 'BDF1':
        mat = assemble_blockwise_matrix_BDF1()
        force = assemble_blockwise_force_BDF1(ux_n,uy_n,xs_n,ys_n)
    elif ph.time_integration == 'BDF2':
        if cn_time == 0:
            mat = assemble_blockwise_matrix_BDF1()
            force = assemble_blockwise_force_BDF1(ux_n,uy_n,xs_n,ys_n)
        else:
            mat = assemble_blockwise_matrix_BDF2()
            force = assemble_blockwise_force_BDF2(
                    ux_n,uy_n,ux_n_old,uy_n_old,xs_n,ys_n,xs_n_old,ys_n_old)

    x = np.hstack([xs_n,ys_n])

    sol_t0 = time.time()
    sol = sp_la.spsolve(mat,force)
    sol_t1 = time.time()

    ux_n_old = ux_n
    uy_n_old = uy_n
    xs_n_old = xs_n
    ys_n_old = ys_n

    u = sol[0:2*ndofs_u]
    p = sol[2*ndofs_u:2*ndofs_u+ndofs_p]
    x = sol[2*ndofs_u+ndofs_p:2*ndofs_u+ndofs_p+2*ndofs_s]
    l = sol[2*ndofs_u+ndofs_p+2*ndofs_s:2*ndofs_u+ndofs_p+4*ndofs_s]

    u_n1 = u
    ux_n1 = u[0      :  ndofs_u]
    uy_n1 = u[ndofs_u:2*ndofs_u]

    dx_n1 = x[0      :  ndofs_s]
    dy_n1 = x[ndofs_s:2*ndofs_s]

    dx_n = np.reshape(dx_n1, dx_n.shape)
    dy_n = np.reshape(dy_n1, dx_n.shape)

    xs_n1 = xs_zero+dx_n
    ys_n1 = ys_zero+dy_n

    xs_n = np.reshape(xs_n1, xs_n.shape)
    ys_n = np.reshape(ys_n1, ys_n.shape)
    ux_n = np.reshape(ux_n1, ux_n.shape)
    uy_n = np.reshape(uy_n1, uy_n.shape)

    str_area = eval_str_area()

    diffusion = str_area/str_area_zero

    p_all_zero = bool(np.all(p==0))
    exploded = bool(np.amax(p) > 1e+10)

    nrg =(l2_norm(FX,x))**2 + (l2_norm(Mv,u))**2
    energy.append(nrg)

    if (exploded==True or p_all_zero == True):
        diffusion = 999
    print '--------------------------------------'
    print 'cn_time   = ' + str(cn_time)
    print 'diffusion = ' + str(diffusion)
    print 'energy    = ' + str(nrg)
    print 'pressure == 0? ' + str(p_all_zero)
    print 'exploded     ? ' + str(exploded)
    print '--------------------------------------'

    #if diffusion > 2:
    #    break
    #elif diffusion < .8:
    #    break

    if ph.stampa[cn_time] == True:
        write_output()
    step_t1 = time.time()
    print 'step time = ' + str((step_t1-step_t0))
    print 'sol  time = ' + str((sol_t1-sol_t0))
