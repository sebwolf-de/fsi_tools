# fsi_tools #

![Alt cover picture](title.png?raw=true "Title")

In this repository I share the code I produce during my master's thesis. It is based on Nicola Cavalini's work.

## Required Dependencies ##
This pieces of code requires the following libraries.* At the current stage, the standard packages coming with and Ubuntu distribution should be ok.*

1. Numpy X.X.X
2. Scipy X.X.X
3. Matplotlib 1.3.0
4. Sphinx X.X.X
5. Shapely X.X.X

Sebastian: I run Fedora and have installed the standard packages. Everything works fine 

## How do I get set up? ##

This not really "software", this more "code to be shared", so you do not really need to install it. Nevertheless it relies on several popular packages that have to installed. 

### Environment variables: ###
We need to tell python where our modules are. We need to add the ``fsi_tools/modules`` directory to the python path.  This is done adding the the following line to ``/home/username/.bashrc``:
```
export PYTHONPATH=$PYTHONPATH:/home/username/path_to_my_dir/fsi_tools/modules
```
Check that everything is all-right typing (still in the shell):
```
echo %PATH%
```
# Time to work, finally #

Once the configuration is done, it is time to get some work done. There are several working files right now:

# Navier Stokes Solver

By invoking 
```
python ns_lid_cavity.py
```
you run a simple solver for Navier-Stokes equation in the well known lid driven cavity example. The script reads `simulation_parameters_ns.json`. The `plot_results_ns.py` creates some nice outputs for the velocity and pressure solution. 

# Distributed Lagrangian Multiplier Scripts

There are two scripts to solve fluid structure interaction problems with the method of Distributed Lagrangian Multipliers.

For `dlm_annulus.py` there is an anulus in a box with homogeneous Dirichlet boundaries. At the beginning the annulus is centered and then moves according to its initial deformation, the elastic forces of the structure and the response of the fluid.

The `dlm_cavity.py` scripts simulates a square initially at rest in a lid driven cavity.

The solver reads the simulation parameters in the `simulation_parameters_fsi.json` file runs the simulation.
The same `simulation_parameters.json` is read by the `plot_results_fsi.py` script. `plot_result_fsis.py` reads the slution parameters and, collects the binary data, and writes plots of the solution

.
