"""Generating and Reading a Mesh. 

.. moduleauthor:: Nicola Cavallini

"""
import lin_tri_mesh as lin_t3
import viewers

if __name__== '__main__':

    nx = 4
    delta_x = 1./nx
    ny = nx
    delta_y = 1./ny
    (topo,x,y) = lin_t3.grid_t3(nx,ny,delta_x,delta_y)

    viewers.tri_plot(x,y,topo)
    
    (topo,x,y) = lin_t3.load_msh('../gmsh_apps/structerd.msh')
    viewers.tri_plot(x,y,topo)
    
    (topo,x,y) = lin_t3.load_msh('../gmsh_apps/structerd_sym.msh')
    viewers.tri_plot(x,y,topo)
    
    (topo,x,y) = lin_t3.load_msh('../gmsh_apps/unstr_square.msh')
    viewers.tri_plot(x,y,topo)