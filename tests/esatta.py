import numpy as np

def sol_esatta(x,y):
	esa=x*(1-x)*y*(1-y)
	return esa

def dx_sol_esatta(x,y):
	dx_esa=(1-2*x)*y*(1-y)
	return dx_esa

def dy_sol_esatta(x,y):
	dy_esa=(1-2*y)*x*(1-x)
	return dy_esa
