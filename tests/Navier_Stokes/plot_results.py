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

from parameters_handler import ParametersHandler

ph = ParametersHandler('simulation_parameters.json')
ph.simulation_info()

time_list = np.where(ph.stampa)[0]

results_dir = ph.results_directory+'/'+ph.sim_prefix+'/binary_data/'

filename = results_dir+'mesh'
f = file(filename,"rb")
topo_p = np.load(f)
x_p = np.load(f)
y_p = np.load(f)
topo_u = np.load(f)
x_u = np.load(f)
y_u = np.load(f)
c2f = np.load(f)
f.close()

pics_dir = ph.results_directory+'/'+ph.sim_prefix+'/pics/'
if not os.path.exists(pics_dir):
    os.makedirs(pics_dir)

vel_dir = pics_dir+'vel/'
if not os.path.exists(vel_dir):
    os.makedirs(vel_dir)

prex_dir = pics_dir+'prex/'
if not os.path.exists(prex_dir):
    os.makedirs(prex_dir)

for cn_time in time_list:
    input_name = results_dir+'cn_time_'+str(cn_time).zfill(ph.time_index_digits)
    print(input_name)
    f = file(input_name,"rb")
    u = np.load(f)
    p = np.load(f)
    f.close()
    print('------------------------------')
    print(input_name)
    #print ' -> ' + output_name
    output_name = vel_dir + 'cn_time_'+str(cn_time).zfill(ph.time_index_digits)
    viewers.quiver_vel(x_u,y_u,u,2*ph.n_delta_x,2*ph.n_delta_x,output_name)
    print(' -> ' + output_name)
    plt.close("all")
    output_name = prex_dir + 'cn_time_'+str(cn_time).zfill(ph.time_index_digits)
    viewers.plot_sol_p1p0_tex(x_p,y_p,p,topo_p,output_name)
    print(' -> ' + output_name)
    plt.close("all")
