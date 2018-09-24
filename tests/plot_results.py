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

import json

with open('simulation_parameters.json') as f:
    params = json.load(f)

base = params["delta_time_base"]
esponente = params["delta_time_negative_esponent"]
kappa = params["structure_stiffness_kappa"]
reynolds = params["reynolds_number"]

n_delta_x = params["fluid_triangulation_intervals"]
n_delta_s = params["structure_triangulation_intervals"]

dt = base*10**(-esponente)
n_times = params["time_intervals_to_be_simulated"]
no_print_intervals = params["time_intervals_in_between_printed_results"]

n_digits = 6

n_times = params["time_intervals_to_be_simulated"]
no_print_intervals = params["time_intervals_in_between_printed_results"]
stampa = []
for i in range(n_times):
    if (i%no_print_intervals==0):
        stampa.append(True)
    else:
        stampa.append(False)

time_list = np.where(stampa)[0]

#time_list = 1*np.arange(0,1)#[0]#10*np.arange(0,2)

bool_conversion = {"true_string" : True, "false_string" : False}
equilibrium_at_zero = bool_conversion[params["equilibrium_at_zero"]]

#solver_type = 'ibm_'
solver_type = 'dlm_'
mesh_prefix = 'str_'
#mesh_prefix = 'thin_'
mesh_name = mesh_prefix+str(int(n_delta_s))+'_'
#mesh_name = 'unstr_32_'
#mesh_name = 'thin_'#'unstr_'+str(int(n_delta_s))
sim_prefix = solver_type
sim_prefix += mesh_name
sim_prefix += 'dt'+str(int(base))+'em'+str(int(esponente))
sim_prefix += '_hx'+str(int(n_delta_x))+'_hs'+str(int(n_delta_s))
sim_prefix += '_k'+str(int(kappa))
sim_prefix += '_re'+str(int(reynolds))
sim_prefix += '_eq_at_zero_'+str(equilibrium_at_zero)

results_dir = 'results/'+sim_prefix+'/binary_data/'

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

pics_dir = sim_prefix+'/pics/'
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

for cn_time in time_list:
    input_name = results_dir+'cn_time_'+str(cn_time).zfill(n_digits)
    f = file(input_name,"rb")
    u = np.load(f)
    p = np.load(f)
    xs = np.load(f)
    ys = np.load(f)
    f.close()
    print '------------------------------'
    print input_name
    #print ' -> ' + output_name
    output_name = str_dir + 'cn_time_'+str(cn_time).zfill(n_digits)
    # if mesh_prefix == 'thin':
    #     viewers.plot_thin_str(xs,ys,output_name)
    # else:
    #     viewers.tri_plot_tex(xs,ys,topo_s,'-','b',output_name)
    print ' -> ' + output_name
    plt.close("all")
    output_name = vel_dir + 'cn_time_'+str(cn_time).zfill(n_digits)
    viewers.quiver_vel(x_u,y_u,u,2*n_delta_x,2*n_delta_x,output_name)
    print ' -> ' + output_name
    plt.close("all")
    output_name = prex_dir + 'cn_time_'+str(cn_time).zfill(n_digits)
    viewers.plot_sol_p1p0_tex(x_p,y_p,p,topo_p,output_name)
    print ' -> ' + output_name
    plt.close("all")
