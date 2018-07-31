# fsi_tools #

This code that I am sharing isn't meant to be, or ether become, a library. On the other hand I needed a modern way to share my code, this is why I created this page. I said this code can't become a library. "Unfortunately" this code started as my Python playground, and I would have never bet a penny that it would have taken me this far. This code isn't enough "thought in advance" to become a library. It is good enough to be shared with students and to prototype. This was my first python code, but #nonsonopoicosipirla.

# Anaconda

[Anaconda](https://www.anaconda.com/download/#linux) solves all dependencies, if you successfully install Anaconda you shouldn't bother
next paragraph.


## Required Dependencies ##
This pieces of code requires the following libraries.* At the current stage, the standard packages coming with and Ubuntu distribution should be ok.*

1. Numpy X.X.X
2. Scipy X.X.X
3. Matplotlib 1.3.0
4. Sphinx X.X.X

Nicola: I have compiled the source for most of them. You are invited to contribute with the version numbers that work for you.

### How to check installation and version: ###

Invoke the python interpreter, type ``python`` on the terminal. You should get something like:

```
Python 2.7.6 (default, Mar 22 2014, 22:59:56)
[GCC 4.8.2] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>>

```
The ``>>>`` is now your command line. If you want to have fun type ``import antigravity``. Now let's get back to work:

```
>>> import numpy as np
>>> np.__version__
```
```

>>> import scipy as sp
>>> sp.__version__
```
```
>>> import matplotlib as mp
>>> mp.__version__
```
```
>>> import sphinx as sx
>>> sx.__version__
```
## How do I get set up? ##

This not really "software", this more "code to be shared", so you do not really need to install it. Nevertheless it relies on several popular packages that have to installed. Let's be more specific.

### How to get the code###
Assuming you have **git** installed (if not see the first point of dependences), and you have a bitbucket account, you can get the code  with:
```
git clone https://your_bitbucket_account_name@bitbucket.org/nicolaalessandro/fsi_tools.git
```

### How to install and upgrade the dependencies: ###
  * git: the code is in a "git repository" so it can be easily shared and contributed. You should have it installed by default. If not:
```
sudo apt-get install git
```
   * numpy: numerical package for python.
```
sudo apt-get install python-numpy
```
   * SciPy: sparse numerical package.
```
sudo apt-get install python-scipy
```
   * matplotlib > 1.3x unfortunately i use triangles, and triplot is available from this version on.
```
sudo apt-get install python-matplotlib
```
   * Sphinx is required to build documentation.
```
sudo apt-get install python-sphinx
```

### Environment variables: ###
We need to tell python where our modules are. We need to add the ``fsi_tools/modules`` directory to the python path.  This is done adding the the following line to ``/home/username/.bascrc``:
```
export PYTHONPATH=$PYTHONPATH:/home/username/path_to_my_dir/fsi_tools/modules
```
For those brave who still run windows in 2015, you can set python path (I copy pasted from [here](https://docs.python.org/2/using/windows.html)) from the shell (that I donno even how to open):
```
set PYTHONPATH=%PYTHONPATH%;C:\path_to_my_dir\fsi_tools\modules
```
Check that everything is all-right typing (still in the shell):
```
echo %PATH%
```
# Time to work, finally #

Once the configuration is done, it is time to get some work done. Fist build the documentation.

### Building the Documentation ###

* Assuming you are in ``fsi_tools`` directory, move to doc, ``cd doc``.
* Type ``make html``.
* Go to ``build/html/`` directory and open ``index.html`` with your browser.

### Running a test ###
The tests are documented in the documentation :)

* Go to ``fsi_tools/tests`` directory.
* Type ``python ex_laplace.py``.

### Everyday workflow with git ###
The user should get some basic git knowledge [him self](http://git-scm.com/book/en/v2). On top of my head I am writing the basic daily git workflow. First you want to update your repository, this is done pulling the latest modifications from the origin and branch master.
```
#!shell
git pull origin master
```
In case you add a file and you want to have it in the git repo, you should add it:
```
#!shell
git add filename
```
Once you have done some changes it is good idea to inspect them:
```
#!shell
git status
```
If you are satisfied with what you did you should commit your work:
```
#!shell
git commit -a -m "briefly describe your work"
```
And if you want to share your commits you can push them on the repo:
```
#!shell
git push origin master
```
