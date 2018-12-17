#! /usr/bin/env python

import json
import matplotlib.pyplot as plt
import math as mth
import numpy as np
import os
import pyamg
from scipy import sparse
import scipy.sparse.linalg as sp_la
import sys
import time

sys.path.append('../../modules')
import assemble
import basis_func as shp
import geom_utils as geom
import la_utils
import lin_tri_mesh as lin_t3
from parameters_handler import ParametersHandler
from preconditioner import BlockPreconditioner
from shapely.geometry import Polygon
import viewers

mat_BDF1 = sparse.load_npz('matrix_BDF1.npz')
mat_BDF2 = sparse.load_npz('matrix_BDF2.npz')
mat_Theta = sparse.load_npz('matrix_Theta.npz')
f = open('rhs', 'rb')
force = np.load(f)
ndofs_u = np.load(f)
ndofs_p = np.load(f)
ndofs_s = np.load(f)
f.close

print(mat_BDF1.shape)
print(force.shape)

# mat_inv = sp_la.inv(mat_BDF1.tocsc())
# _,largest_sv,_ = sp_la.svds(mat_BDF1, k=1)
# _,inv_largest,_ = sp_la.svds(mat_inv, k=1)
# cond_BDF1 = largest_sv*inv_largest
#
# mat_inv = sp_la.inv(mat_BDF2.tocsc())
# _,largest_sv,_ = sp_la.svds(mat_BDF2, k=1)
# _,inv_largest,_ = sp_la.svds(mat_inv, k=1)
# cond_BDF2 = largest_sv*inv_largest
#
# mat_inv = sp_la.inv(mat_Theta.tocsc())
# _,largest_sv,_ = sp_la.svds(mat_Theta, k=1)
# _,inv_largest,_ = sp_la.svds(mat_inv, k=1)
# cond_Theta = largest_sv*inv_largest

# print(cond_BDF1)
print(np.linalg.cond(mat_BDF1.todense()))
# print(cond_BDF2)
print(np.linalg.cond(mat_BDF2.todense()))
# print(cond_Theta)
print(np.linalg.cond(mat_Theta.todense()))

quit()

# A_f = mat[0:2*ndofs_u, 0:2*ndofs_u]
# BT = -mat[0:2*ndofs_u, 2*ndofs_u + ndofs_p]
# LT_f = mat[0:2*ndofs_u, 2*ndofs_u + ndofs_p + 2*ndofs_s:2*ndofs_u + ndofs_p + 4*ndofs_s]
# K_s = mat[2*ndofs_u + ndofs_p:2*ndofs_u + ndofs_p + 2*ndofs_s, 2*ndofs_u + ndofs_p:2*ndofs_u + ndofs_p + 2*ndofs_s]
# LT_s = mat[2*ndofs_u + ndofs_p:2*ndofs_u + ndofs_p + 2*ndofs_s, 2*ndofs_u + ndofs_p + 2*ndofs_s:2*ndofs_u + ndofs_p + 4*ndofs_s]
# B = BT.transpose()
# L_f = LT_f.transpose()
# L_s = LT_s.transpose()
#
# f = force[0:2*ndofs_u]
# g = force[2*ndofs_u + ndofs_p:2*ndofs_u + ndofs_p + 2*ndofs_s]


mat = mat_BDF1.tocsc()
t0 = time.time()
fill_factor = 100
drop_tol = 1e-6
spilu = sp_la.spilu(mat, fill_factor=fill_factor, drop_tol=drop_tol)
M_x = lambda x: spilu.solve(x)
precond = sp_la.LinearOperator(mat.shape, M_x)
t1 = time.time()
print('preconditioner, t = ' + str(t1 - t0))
t0 = time.time()
(sol, it) = sp_la.bicgstab(mat, force, M=precond, tol=1e-8)
print(it)
t1 = time.time()
print('bicgstab, t = ' + str(t1-t0))


t0 = time.time()
x = sp_la.spsolve(mat, force)
t1 = time.time()

print('direct solver, t = ' + str(t1-t0))
