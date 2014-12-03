# fsi_tools #

This code that I am sharing isn't meant to be or ether become a library. On the other hand I needed a modern way to share my code this is why I created this page. I said this code can't become a library. "Unfortunately" this code started as my Python playground and I would have never bet a penny that it would have taken me this far. This code isn't enough "thought in advance" to become a library. It is good enough to be shared with students and to prototype. This was my first python code, but #nonsonopoicosipirla.  

## How do I get set up? ##

This not really "software", this more "code to be shared", so you do not really need to install it. Nevertheless it relies on several popular packages that have to installed. Let's be more specific. 

### How to get the code###
Assuming you have **git** installed (if not see the first point of dependences) you can get the code  with: 
```
#!shell
git clone https://nicolaalessandro@bitbucket.org/nicolaalessandro/fsi_tools.git
```

### Dependencies: ###
  * git: the code is in a "git repository" so it can be easily shared and contributed. You should have it installed by default. If not: 
```
#!shell
sudo apt-get install git
```
   * numpy: numerical package for python. 
```
#!shell
sudo apt-get install python-numpy
```
   * SciPy: sparse numerical package. 
```
#!shell
sudo apt-get install python-scipy
```
   * matplotlib. plotting. 
```
#!shell
sudo apt-get install python-matplotlib
```

### Building the Documentation ###

* Assuming you are in ``fsi_tools`` directory, move to doc, ``cd doc``.
* Type ``make html``.
* Go to ``build/html/`` directory and open ``index.html`` with your browser. 

### Running a test ###

* Go to ``fsi_tools/tests`` directory.
* Type ``python ex_laplace.py``.