# fsi_tools #

![Alt cover picture](title.png?raw=true "Title")

In this repository I share the code I produce during my master's thesis. It is based on Nicola Cavalini's work.

## Required Dependencies ##
This pieces of code requires the following libraries.

1. Numpy X.X.X
2. Scipy X.X.X
3. Matplotlib 1.3.0
4. Sphinx X.X.X
5. Shapely X.X.X
6. Rtree X.X.X

You can either install the packages via your distribution's packaging tool or via pip.

## How do I get set up? ##

This not really "software", this more "code to be shared", so you do not really need to install it. Nevertheless it relies on several popular packages that have to installed.

### Environment variables: ###
We need to tell python where our modules are. We need to add the ``fsi_tools/modules`` directory to the python path.  This is done adding the the following line to ``/home/username/.bashrc``:
```
export PYTHONPATH=$PYTHONPATH:/home/username/path_to_my_dir/fsi_tools/modules
```

# Time to work, finally #

Once the configuration is done, it is time to get some work done. There are several working files right now:

# Basics

In the folder

```
test/Basics
```

there are several scripts, that are based on basic numerical algorithms.

# Navier Stokes Solver

This should be a solver for Navier-Stokes equation only. It is still experimental.

By invoking
```
python ns_master.py
```
you run a solver for Navier-Stokes equation. The script reads `simulation_parameters_ns.json`. The `plot_results_ns.py` creates some nice outputs for the velocity and pressure solution. With `export_to_vtk.py` you can export the results to view them in an external program e.g. paraview.

# Distributed Lagrangian Multiplier Scripts

Right now, there are two examples implemented: An initially deformed annulus, that moves back to its undeformed position, a small disk floating in a lid driven cavity.

There are two scripts to solve fluid structure interaction problems with the method of Distributed Lagrangian Multipliers.

`dlm_semiimplicit.py` uses a semiimplicit time stepping method, `dlm_fixpoint.py` is a fully implicit method, which uses a fixpoint method to solve the nonlinear equations.

The solver reads the simulation parameters in the `simulation_parameters_fsi.json` file, where you can define material parameters and solver parameters.

`plot_results.py` creates plots, `export_to_vtk.py` exports the results for an external viewer.

There are also some scripts for evaluations and convergence analysis included.
