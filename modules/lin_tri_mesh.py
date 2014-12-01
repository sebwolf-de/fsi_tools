import numpy as np
import math as mth

# ================================================
# define a linear triangular uniform triangulation
# ================================================

def mesh_t3(nx,ny,delta_x,delta_y):
    topo = np.zeros((2*nx*ny,3),dtype=int)
    x = np.linspace(0,nx*delta_x,nx+1)
    y = np.linspace(0,ny*delta_y,ny+1)
    (x,y) = np.meshgrid(x,y)
    x = np.concatenate(x)
    y = np.concatenate(y)
    i_el = 0
    for i in range(ny):
        for j in range(nx):
            nodes = np.array([ i    * (nx + 1 ) + j,
                              (i+1) * (nx + 1 ) + j + 1,
                              (i+1) * (nx + 1 ) + j])
            topo[i_el][:] = nodes
            i_el = i_el+1
            nodes = np.array([ i * (nx + 1) + j,
                               i * (nx + 1) + j + 1,
                            (i+1)* (nx + 1) + j + 1])
            topo[i_el][:] = nodes
            i_el = i_el+1
    
    return topo, x, y

def mesh_t3_t0(nx,ny,delta_x,delta_y):
    (topo, x, y) = mesh_t3(nx,ny,delta_x,delta_y)
    tmp = np.arange((nx+1)*(ny+1),(nx+1)*(ny+1)+2*nx*ny)
    tmp = np.reshape(tmp,(2*nx*ny,1))
    topo = np.hstack([topo,tmp])
    return topo, x, y


def grid_t3(nx,ny,delta_x,delta_y):
    topo = np.zeros((2*nx*ny,3),dtype=int)
    x = np.linspace(0,nx*delta_x,nx+1)
    y = np.linspace(0,ny*delta_y,ny+1)
    (x,y) = np.meshgrid(x,y)
    x = np.concatenate(x)
    y = np.concatenate(y)
    i_el = 0
    for i in range(ny):
        for j in range(nx):
            nodes = np.array([ i * (nx + 1 ) + j,
                          i * (nx + 1 ) + j + 1,
                        (i+1)*(nx + 1 ) + j])
            
            topo[i_el][:] = nodes
            i_el = i_el+1
            nodes = np.array([ i * (nx + 1 ) + j + 1,
                        (i+1)*(nx + 1 ) + j + 1,
                        (i+1)*(nx + 1 ) + j])
            topo[i_el][:] = nodes
            i_el = i_el+1
    
    return topo, x, y

def grid_t3_t0(nx,ny,delta_x,delta_y):
    topo = np.zeros((2*nx*ny,4),dtype=int)
    x = np.linspace(0,nx*delta_x,nx+1)
    y = np.linspace(0,ny*delta_y,ny+1)
    (x,y) = np.meshgrid(x,y)
    x = np.concatenate(x)
    y = np.concatenate(y)
    i_el = 0
    for i in range(ny):
        for j in range(nx):
            nodes = np.array([ i * (nx + 1 ) + j,
                          i * (nx + 1 ) + j + 1,
                        (i+1)*(nx + 1 ) + j,
                        i_el+(nx+1)*(ny+1)])
            
            topo[i_el][:] = nodes
            i_el = i_el+1
            nodes = np.array([ i * (nx + 1 ) + j + 1,
                        (i+1)*(nx + 1 ) + j + 1,
                        (i+1)*(nx + 1 ) + j,
                        i_el+(nx+1)*(ny+1)])
            topo[i_el][:] = nodes
            i_el = i_el+1
    
    return topo, x, y


# ================================================

def mesh_t3_iso_t6(nx,ny,delta_x,delta_y):

    (tt3,xt3,yt3) = mesh_t3(nx,ny,delta_x,delta_y)

    (tt6,xt6,yt6) = mesh_t3(2*nx,2*ny,delta_x/2,delta_y/2)

    corse_to_fine = np.zeros((2*nx*ny,4),dtype=int)

    for i in range(ny):
        for j in range(nx):
            lr_coarse = (2*nx)*i+2*j
            fine = np.zeros((1,4))
            fine[0][0] = (j*4)+8*nx*i+1
            fine[0][1] = (j*4)+8*nx*i+2
            fine[0][2] = (j*4)+8*nx*i+3
            fine[0][3] = (j*4)+8*nx*i+4*nx+3
            corse_to_fine[lr_coarse] = fine
            ul_coarse = (2*nx)*i+2*j+1
            fine[0][0] = (j*4)+8*nx*i
            fine[0][1] = (j*4)+8*nx*i+4*nx
            fine[0][2] = (j*4)+8*nx*i+4*nx+1
            fine[0][3] = (j*4)+8*nx*i+4*nx+2
            corse_to_fine[ul_coarse] = fine
    return tt3,xt3,yt3,tt6,xt6,yt6,corse_to_fine

def grid_t3_iso_t6(nx,ny,delta_x,delta_y):

    (tt3,xt3,yt3) = grid_t3(nx,ny,delta_x,delta_y)

    (tt6,xt6,yt6) = grid_t3(2*nx,2*ny,delta_x/2,delta_y/2)

    corse_to_fine = np.zeros((2*nx*ny,4),dtype=int)

    for i in range(ny):
        for j in range(nx):
            ll_coarse = (2*nx)*i+2*j
            fine = np.zeros((1,4))
            fine[0][0] = (j*4)+8*nx*i
            fine[0][1] = (j*4)+8*nx*i+1
            fine[0][2] = (j*4)+8*nx*i+2
            fine[0][3] = (j*4)+8*nx*i+4*nx
            corse_to_fine[ll_coarse] = fine
            ur_coarse = (2*nx)*i+2*j+1
            fine[0][0] = (j*4)+8*nx*i+3
            fine[0][1] = (j*4)+8*nx*i+4*nx+1
            fine[0][2] = (j*4)+8*nx*i+4*nx+2
            fine[0][3] = (j*4)+8*nx*i+4*nx+3
            corse_to_fine[ur_coarse] = fine
    
    return tt3,xt3,yt3,tt6,xt6,yt6,corse_to_fine

# ================================================

def ibm_mesh(node_s,ray,deform):

    s_lgr = np.linspace(0,1.-1./node_s,node_s)
    s_lgr = 2*mth.pi * ray * s_lgr

    x_str = .5 + ray * np.cos(s_lgr/ray+mth.pi/4)
    y_str = .5 + (ray-deform) * np.sin(s_lgr/ray+mth.pi/4)

    XY_str = np.zeros((s_lgr.shape[0],2))

    XY_str[:,0] = x_str
    XY_str[:,1] = y_str
    return s_lgr, XY_str

def lin_str_mesh(node_s,ray,deform):

    s_lgr = np.linspace(0,1,node_s+1)
    
    s_lgr = mth.pi/2 * ray * s_lgr
    x_str = ray * np.cos(s_lgr/ray)
    y_str = ray * np.sin(s_lgr/ray)

    x_str = 1./deform*(x_str)
    y_str =    deform*(y_str)
    

    topo = np.zeros((node_s,2), dtype='int')
    ieq = np.zeros((node_s+1,1), dtype='int')

    for cnt in range(0,node_s):
        topo[cnt,0] = cnt
        topo[cnt,1] = cnt+1
        #if cnt == node_s-1:
        #    ieq[cnt] = cnt
        #    ieq[cnt+1] = 0
        #else:
        ieq[cnt] = cnt
        ieq[cnt+1] = cnt+1

    return topo, x_str, y_str, s_lgr, ieq

def load_msh(filename):
    f = open ( filename , 'r')

    x = np.array([])
    y = np.array([])
    topo = np.array([], dtype=int)

    for line in f:
        #print line[0]
        if line[0]=='$':
            print 'non fare un cippa'
        else:
            l = map(float,line.split())
            if len(l) == 4:
                x = np.append(x,l[1])
                y = np.append(y,l[2])
                #print 'ciao'
            elif len(l) == 8:
                row = l[5:8]
                row = np.array(row, dtype=int)
                topo = np.append(topo,row)
                #print row
            #    
            #print len(l)
    topo = np.reshape(topo,(len(topo)/3,3))
    topo = topo-1
    r_id = 0 
    for row in topo:
        ck =      (x[row[1]]-x[row[0]])*(y[row[2]]-y[row[0]])
        ck = ck - (x[row[2]]-x[row[0]])*(y[row[1]]-y[row[0]])
        if ck < 0:
            topo[r_id,:] = np.array([[row[0],row[2],row[1]]])
        r_id+=1
    return topo,x,y
