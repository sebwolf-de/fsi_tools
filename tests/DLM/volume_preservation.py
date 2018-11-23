#! /usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy import sparse
from shapely.geometry import Polygon
import sys
from vtk import *

import sys
import assemble

def read_area(k):
    reader = vtk.vtkXMLUnstructuredGridReader()
    filename = '../../../ans-ifem/ans-ifem/out/Cavity-solid-'+str(k).zfill(5)+'.vtu'
    reader.SetFileName(filename)
    reader.Update() # Needed because of GetScalarRange

    unstrGrid = reader.GetOutput()

    numcells = unstrGrid.GetNumberOfCells()
    area = 0
    for j in range(numcells):
        p_0 = (unstrGrid.GetCell(j).GetPoints().GetPoint(0)[0], unstrGrid.GetCell(j).GetPoints().GetPoint(0)[1])
        p_1 = (unstrGrid.GetCell(j).GetPoints().GetPoint(1)[0], unstrGrid.GetCell(j).GetPoints().GetPoint(1)[1])
        p_2 = (unstrGrid.GetCell(j).GetPoints().GetPoint(2)[0], unstrGrid.GetCell(j).GetPoints().GetPoint(2)[1])
        p_3 = (unstrGrid.GetCell(j).GetPoints().GetPoint(3)[0], unstrGrid.GetCell(j).GetPoints().GetPoint(3)[1])


        poly = Polygon([p_0, p_1, p_2, p_3])
        area += poly.area
    return area


def eval_str_area(k):
    input_name = results_dir+'cn_time_'+str(k).zfill(3)
    f = file(input_name,"rb")
    u = np.load(f)
    p = np.load(f)
    x_s = np.load(f)
    y_s = np.load(f)
    f.close()
    area = 0
    for row in topo_s:
        eval_p = np.zeros((3,2))
        eval_p[:,0] = x_s[row]
        eval_p[:,1] = y_s[row]
        poly = Polygon(tuple(eval_p.tolist()))
        area+= poly.area
    return area

def l2_norm(M,g):
    l2_g = M.dot(g)
    l2_g = np.dot(l2_g.transpose(),g)
    l2_g = math.sqrt(l2_g)
    return l2_g

if len(sys.argv) > 2:
    results_dir = sys.argv[1]
    n = int(sys.argv[2])
else:
    n = 800

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
x_s = np.load(f)
y_s = np.load(f)
s_lgr = np.load(f)
t_lgr = np.load(f)
f.close()


# ie_s = np.arange(0,s_lgr.shape[0])
# KS11 = assemble.gradu_gradv_p1_ieq(topo_s,s_lgr,t_lgr,ie_s)
# MF11 = assemble.u_v_p1(topo_u,x_u,y_u)
# KS = sparse.vstack([
#     sparse.hstack([KS11, sparse.csr_matrix(KS11.shape)]),
#     sparse.hstack([sparse.csr_matrix(KS11.shape), KS11])
# ])
# MF = sparse.vstack([
#     sparse.hstack([MF11, sparse.csr_matrix(MF11.shape)]),
#     sparse.hstack([sparse.csr_matrix(MF11.shape), MF11])
# ])



diffusion = np.zeros((2, n))
energy = np.zeros((n))

str_area_zero = eval_str_area(0)
deal_area_zero = read_area(0)
# dx_n = sx_n - s_lgr
# dy_n = sy_n - t_lgr
#energy[0] =(l2_norm(KS,np.append(dx_n, dy_n)))**2 + (l2_norm(MF,u))**2
for cn_time in range(1, n):
    str_area = eval_str_area(cn_time)
    diffusion[0, cn_time] = (str_area/str_area_zero - 1.)*100
    deal_area = read_area(cn_time)
    diffusion[1, cn_time] = (deal_area/deal_area_zero - 1.)*100
    # dx_n = sx_n - s_lgr
    # dy_n = sy_n - t_lgr
    #energy[cn_time] =(l2_norm(KS,np.append(dx_n, dy_n)))**2 + (l2_norm(MF,u))**2
    #print str(cn_time).zfill(3) + ' ' + str(diffusion[cn_time]) + ' ' + str(energy[cn_time])

plt.plot(np.arange(0,n), diffusion[0,:])#, diffusion[1,:])
plt.xlabel('time (s)')
plt.ylabel('volume change (%)')
plt.title('Volume preservation for the disk example (BDF1)')
plt.grid(True)
plt.show()

# plt.plot(np.arange(0,n), energy)
# plt.xlabel('time (s)')
# plt.ylabel('energy')
# plt.title('Energy preservation for the disk example (BDF1)')
# plt.grid(True)
# plt.show()
