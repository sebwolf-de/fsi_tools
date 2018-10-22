#! /usr/bin/env python

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# nicola modules
import lin_tri_mesh as lin_t3
import basis_func as shp
import assemble
import la_utils
import viewers
import geom_utils as geom

results_dir = 'results/Convergence_Analysis_Annulus/'

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

pics_dir = results_dir+'pics/'
if not os.path.exists(pics_dir):
    os.makedirs(pics_dir)

vel_dir = pics_dir+'vel/'
if not os.path.exists(vel_dir):
    os.makedirs(vel_dir)

prex_dir = pics_dir+'prex/'
if not os.path.exists(prex_dir):
    os.makedirs(prex_dir)

str_dir = pics_dir+'str/'
if not os.path.exists(str_dir):
    os.makedirs(str_dir)

N = int(np.sqrt(x_u.shape)) - 1

for k in range(1,6):
    input_name = results_dir+'BDF1_dt=1_'+str(2**k)+'_result'
    f = file(input_name,"rb")
    u = np.load(f)
    p = np.load(f)
    xs = np.load(f)
    ys = np.load(f)
    f.close()
    print '------------------------------'
    print input_name
    output_name = str_dir+'BDF1_dt=1_'+str(2**k)+'_result'
    viewers.tri_plot_tex(xs,ys,topo_s,'-b',output_name)
    print ' -> ' + output_name
    plt.close("all")
    output_name = vel_dir+'BDF1_dt=1_'+str(2**k)+'_result'
    viewers.quiver_vel(x_u,y_u,u,N,N,output_name)
    print ' -> ' + output_name
    plt.close("all")
    output_name = prex_dir+'BDF1_dt=1_'+str(2**k)+'_result'
    viewers.plot_sol_p1p0_tex(x_p,y_p,p,topo_p,output_name)
    print ' -> ' + output_name
    plt.close("all")

    input_name = results_dir+'BDF2_dt=1_'+str(2**k)+'_result'
    f = file(input_name,"rb")
    u = np.load(f)
    p = np.load(f)
    xs = np.load(f)
    ys = np.load(f)
    f.close()
    print '------------------------------'
    print input_name
    output_name = str_dir+'BDF2_dt=1_'+str(2**k)+'_result'
    viewers.tri_plot_tex(xs,ys,topo_s,'-b',output_name)
    print ' -> ' + output_name
    plt.close("all")
    output_name = vel_dir+'BDF2_dt=1_'+str(2**k)+'_result'
    viewers.quiver_vel(x_u,y_u,u,N,N,output_name)
    print ' -> ' + output_name
    plt.close("all")
    output_name = prex_dir+'BDF2_dt=1_'+str(2**k)+'_result'
    viewers.plot_sol_p1p0_tex(x_p,y_p,p,topo_p,output_name)
    print ' -> ' + output_name
    plt.close("all")

    input_name = results_dir+'Theta_dt=1_'+str(2**k)+'_result'
    f = file(input_name,"rb")
    u = np.load(f)
    p = np.load(f)
    xs = np.load(f)
    ys = np.load(f)
    f.close()
    print '------------------------------'
    print input_name
    output_name = str_dir+'Theta_dt=1_'+str(2**k)+'_result'
    viewers.tri_plot_tex(xs,ys,topo_s,'-b',output_name)
    print ' -> ' + output_name
    plt.close("all")
    output_name = vel_dir+'Theta_dt=1_'+str(2**k)+'_result'
    viewers.quiver_vel(x_u,y_u,u,N,N,output_name)
    print ' -> ' + output_name
    plt.close("all")
    output_name = prex_dir+'Theta_dt=1_'+str(2**k)+'_result'
    viewers.plot_sol_p1p0_tex(x_p,y_p,p,topo_p,output_name)
    print ' -> ' + output_name
    plt.close("all")

input_name = results_dir+'reference'
f = file(input_name,"rb")
u = np.load(f)
p = np.load(f)
xs = np.load(f)
ys = np.load(f)
f.close()
print '------------------------------'
print input_name
output_name = str_dir+'reference'
viewers.tri_plot_tex(xs,ys,topo_s,'-b',output_name)
print ' -> ' + output_name
plt.close("all")
output_name = vel_dir+'reference'
viewers.quiver_vel(x_u,y_u,u,N,N,output_name)
print ' -> ' + output_name
plt.close("all")
output_name = prex_dir+'reference'
viewers.plot_sol_p1p0_tex(x_p,y_p,p,topo_p,output_name)
print ' -> ' + output_name
plt.close("all")
