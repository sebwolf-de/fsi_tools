#! /bin/env python

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

import scipy.sparse.linalg as sp_la
import sys
import time

from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import LineString

sys.path.append('../../modules')
import assemble
import basis_func as basis
import la_utils
import lin_tri_mesh as lin_t3
import quadrature


n = 10
dx = 1./n

(topo,x,y) = lin_t3.mesh_t3(n,n,dx,dx)

K = assemble.gradu_gradv_p1(topo,x,y)
M = assemble.u_v_p1(topo,x,y)

bc_id_unten = np.where(y < dx/10)
bc_id_oben = np.where(y > 1-dx/10)
bc_id_links = np.where(x < dx/10)
bc_id_rechts = np.where(x > 1-dx/10)

K = la_utils.set_diag(K,bc_id_unten)
K = la_utils.set_diag(K,bc_id_oben)
K = la_utils.set_diag(K,bc_id_links)
K = la_utils.set_diag(K,bc_id_rechts)

c = 200
f = -c*(2-12*x+12*x**2)*(1-y)**2*y**2 - c*(2-12*y+12*y**2)*(1-x)**2*x**2

f = M.dot(f)

f[bc_id_unten] = 0
f[bc_id_oben] = 0
f[bc_id_links] = 0
f[bc_id_rechts] = 0

es = lambda x,y: np.reshape(c*(1-x)**2*x**2 * (1-y)**2*y**2, (x.shape[0]))
es_x = lambda x,y: c*(2*x-6*x**2+4*x**3)*(1-y)**2*y**2
es_y = lambda x,y: c*(1-x)**2*x**2*(2*y-6*y**2+4*y**3)
es_eval = es(x,y)


sol = sp_la.spsolve(K, f)

err_inf = np.linalg.norm(es_eval - sol, float('inf'))
err_l2 = quadrature.l2error_on_mesh(sol, es, x, y, topo)
err_h1 = np.sqrt(err_l2**2 + quadrature.h1_semi_error_on_mesh(sol, es_x, es_y, x, y, topo)**2)
print('inf',err_inf)
print('l2 ',err_l2[0])
print('h1 ',err_h1[0])

X = np.reshape(x, (n+1,n+1))
Y = np.reshape(y, (n+1,n+1))
Z = np.reshape(sol, (n+1,n+1))



fig = plt.figure()
ax = fig.gca(projection='3d')

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
m = np.max(np.abs(Z))
ax.set_zlim(-1.01*m, 1.01*m)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
