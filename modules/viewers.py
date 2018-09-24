from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np



def tri_plot(x,y,topo):
    """
    This is a super simple function to plot a triangular mesh.
    It is also an interesting exercize to see how to loop over the
    mesh elementzs. ``x`` and ``y`` are one dimensional arrays containg the nodes
    coordinates. On each row of the ``topo`` matrix is stored the
    dof (degree of freedom) with support on the row-th element. Notice that
    we are looping over the first dimension of ``topo``, so we only need to
    write a standard ``python`` for loop:

    .. code:: python

        for row in topo:

    The ``numpy`` command ``hstack`` is used to attach the last node to the
    first one. This is done to "close" our triangles.

    .. code:: python

       row = np.hstack([row,row[0]])

    We use the ``[`` ``]`` brackets to acces the nodes and define the local coords.

    .. code:: python

       x_l = x[row]
       y_l = y[row]

    Now we only need to plot the local coords for every triangle:

    .. code:: python

       plt.plot(x_l, y_l,'-b',linewidth=2)

    """
    for row in topo:
        row = np.hstack([row,row[0]])
        x_l = x[row]
        y_l = y[row]
        plt.plot(x_l, y_l,'-b',linewidth=2)


    plt.show()

    return


def plot_sol_p1(x,y,z,topo):

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.plot_trisurf(x, y, z, cmap=cm.jet,vmin=min(z), vmax=max(z), linewidth=0.2)

    #for row in topo:
    #    x_l = x[row]
    #    y_l = y[row]
    #    z_l = z[row]
    #    ax.plot_trisurf(x_l, y_l, z_l, cmap=cm.jet,
    #                    vmin=min(z), vmax=max(z), linewidth=0.2)
    #ax.view_init(90, 0)
    #ax.set_zlim([-20,20])
    plt.show()
    return

def plot_sol_p1p0(x,y,z,topo):

    fig = plt.figure()
    ax = fig.gca(projection='3d', aspect='equal')

    #ax.plot_trisurf(x, y, z, cmap=cm.jet,vmin=min(z), vmax=max(z), linewidth=0.2)

    for row in topo:
        x_l = x[row[0:3]]
        y_l = y[row[0:3]]
        z_l = z[row[0:3]]+z[row[3]]
        max_z = max(z_l)
        ax.plot_trisurf(x_l, y_l, z_l, cmap=cm.jet,
                        vmin=-110, vmax=110, linewidth=0.)
    ax.view_init(70,40)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    #ax.axis('equal')
    plt.show()
    return
def plot_sol_p1p0_tex(x,y,z,topo,filename):

    fig = plt.figure()
    ax = fig.gca(projection='3d', aspect='equal')

    #ax.plot_trisurf(x, y, z, cmap=cm.jet,vmin=min(z), vmax=max(z), linewidth=0.2)

    for row in topo:
        x_l = x[row[0:3]]
        y_l = y[row[0:3]]
        z_l = z[row[0:3]]+z[row[3]]
        max_z = max(z_l)
        ax.plot_trisurf(x_l, y_l, z_l, cmap=cm.jet,
                        vmin=-20, vmax=20., linewidth=0.)
    font_size = 18
    #ax.view_init(50,30)
    ax.set_xlabel(r'$x$',fontsize=font_size)
    ax.set_ylabel(r'$y$',fontsize=font_size)
    ax.set_zlabel(r'$p$',fontsize=font_size)
    ax.set_xticks([0,0.357,.5,1])
    ax.set_xticklabels([r'$0$', r'$0.357$', r'$0.5$',r'$1$'],
                       fontsize=font_size)
    ax.set_yticks([0,.5,.7,1])
    ax.set_yticklabels([r'$0$', r'$0.5$', r'$0.7$',r'$1$'],
                       fontsize=font_size)
    ax.set_zticks([-20,0,20])
    ax.set_zticklabels(
        [r'$-20$',r'$0$',r'$20$'],
        fontsize=font_size)
    ax.set_ylim([0,1])
    ax.set_xlim([0,1])
    ax.set_zlim([-20,20])
    ax.view_init(80,60)
    ax.axis('equal')

    fig.savefig(filename+'_prex.png')
    return

def tri_plot_nodes_num(x,y,topo,line):
    i_el = 0
    fig = plt.figure()
    ax = fig.gca(aspect='equal')
    for row in topo:
        nds = np.hstack([row,row[0]])
        x_l = x[nds]
        y_l = y[nds]
        x_g = sum(x[row])/3
        y_g = sum(y[row])/3
        #plt.text(x_g,y_g,str(i_el))
        plt.plot(x_l, y_l,line,linewidth=2)
        i_el +=1

    node_id = 0
    for xn,yn in zip(x,y):
        ax.text(xn,yn,str(node_id))
        node_id +=1


    #plt.axis([min(x)-.1*max(x), 1.1*max(x),
    #          min(y)-.1*max(y), 1.1*max(y)])

    #plt.savefig(filename)

    return

def tri_plot_tex(x,y,topo,line,filename):
    fig = plt.figure()

    ax = fig.gca(aspect='equal')
    for row in topo:
        nds = np.hstack([row,row[0]])
        x_l = x[nds]
        y_l = y[nds]
        plt.plot(x_l, y_l,line,linewidth=1)

    ax.set_ylim([0,1])
    ax.set_xlim([0,1])
    font_size = 18
    ax.set_xlabel(r'$x$',fontsize=font_size)
    ax.set_ylabel(r'$y$',fontsize=font_size)
    ax.set_xticks([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])
    ax.set_xticklabels([r'$0$',r'$0.1$',r'$0.2$',r'$0.3$',r'$0.4$' ,r'$0.5$',r'$0.6$',r'$0.7$',r'$0.8$',r'$0.9$',r'$1$'],
                       fontsize=font_size)
    ax.set_yticks([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])
    ax.set_yticklabels([r'$0$',r'$0.1$',r'$0.2$',r'$0.3$',r'$0.4$' ,r'$0.5$',r'$0.6$',r'$0.7$',r'$0.8$',r'$0.9$',r'$1$'],
                       fontsize=font_size)

    ax.grid(True)
    fig.savefig(filename+'.png')

    #plt.savefig(filename)

    return

def plot_thin_str(x,y,filename):
    fig = plt.figure()
    ax = fig.gca(aspect='equal')
    plt.plot(x,y)

    ax.set_ylim([0,1])
    ax.set_xlim([0,1])


    fig.savefig(filename+'.png')
    return

def plot_ibm(XY_str):
    x = XY_str[:,0]
    y = XY_str[:,1]
    plt.plot(x,y)
    plt.axis('equal')
    return

def quiver_vel(x,y,u,nx,ny,filename):
    ndofs = u.shape[0]

    ux = u[0:ndofs/2]
    uy = u[ndofs/2:ndofs]

    #print u.shape
    #print ux.shape
    #print uy.shape

    x = np.reshape(x , (nx+1,ny+1))
    y = np.reshape(y , (nx+1,ny+1))
    ux = np.reshape(ux , (nx+1,ny+1))
    uy = np.reshape(uy , (nx+1,ny+1))

    fig = plt.figure()

    plt.quiver(x,y,ux,uy)
    plt.axis([-.1, 1.1, -.1, 1.1])
    fig.savefig(filename+'.png')

    return

def streamlines_str_plot(x,y,u,nx,ny,xs,ys,topo_s,filename):

    ndofs = u.shape[0]

    ux = u[0:ndofs/2]
    uy = u[ndofs/2:]

    x = np.reshape(x , (nx+1,ny+1))
    y = np.reshape(y , (nx+1,ny+1))
    ux = np.reshape(ux , (nx+1,ny+1))
    uy = np.reshape(uy , (nx+1,ny+1))

    fig = plt.figure()
    ax = fig.gca(aspect='equal')

    speed = np.sqrt(ux*ux+uy*uy)

    lw = 5*speed/speed.max()
    plt.streamplot(x,y,ux,uy,color=speed,linewidth=lw)

    for row in topo_s:
        nds = np.hstack([row,row[0]])
        x_l = xs[nds]
        y_l = ys[nds]
        plt.plot(x_l, y_l,color='0.',linewidth=1)

    ax.set_ylim([0,1])
    ax.set_xlim([0,1])
    font_size = 18
    ax.set_xlabel(r'$x$',fontsize=font_size)
    ax.set_ylabel(r'$y$',fontsize=font_size)
    ax.set_xticks([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])
    ax.set_xticklabels([r'$0$',r'$0.1$',r'$0.2$',r'$0.3$',r'$0.4$' ,r'$0.5$',r'$0.6$',r'$0.7$',r'$0.8$',r'$0.9$',r'$1$'],
                       fontsize=font_size)
    ax.set_yticks([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])
    ax.set_yticklabels([r'$0$',r'$0.1$',r'$0.2$',r'$0.3$',r'$0.4$' ,r'$0.5$',r'$0.6$',r'$0.7$',r'$0.8$',r'$0.9$',r'$1$'],
                       fontsize=font_size)

    ax.grid(True)
    plt.show()
    fig.savefig(filename+'.png')

    return

def vel_str_plot(x,y,u,nx,ny,xs,ys,topo_s,filename):
    ndofs = u.shape[0]

    ux = u[0:ndofs/2]
    uy = u[ndofs/2:ndofs]

    #print u.shape
    #print ux.shape
    #print uy.shape

    x = np.reshape(x , (nx+1,ny+1))
    y = np.reshape(y , (nx+1,ny+1))
    ux = np.reshape(ux , (nx+1,ny+1))
    uy = np.reshape(uy , (nx+1,ny+1))

    fig = plt.figure()
    ax = fig.gca(aspect='equal')

    plt.quiver(x,y,ux,uy,color='0.5')
    #plt.axis([-.1, 1.1, -.1, 1.1])
    #fig.savefig(filename+'.png')

    #    fig = plt.figure()

    for row in topo_s:
        nds = np.hstack([row,row[0]])
        x_l = xs[nds]
        y_l = ys[nds]
        plt.plot(x_l, y_l,'-b',linewidth=1)

    ax.set_ylim([-.1,1.1])
    ax.set_xlim([-.1,1.1])
    font_size = 18
    ax.set_xlabel(r'$x$',fontsize=font_size)
    ax.set_ylabel(r'$y$',fontsize=font_size)
    ax.set_xticks([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])
    ax.set_xticklabels([r'$0$',r'$0.1$',r'$0.2$',r'$0.3$',r'$0.4$' ,r'$0.5$',r'$0.6$',r'$0.7$',r'$0.8$',r'$0.9$',r'$1$'],
                       fontsize=font_size)
    ax.set_yticks([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])
    ax.set_yticklabels([r'$0$',r'$0.1$',r'$0.2$',r'$0.3$',r'$0.4$' ,r'$0.5$',r'$0.6$',r'$0.7$',r'$0.8$',r'$0.9$',r'$1$'],
                       fontsize=font_size)

    ax.grid(True)
    plt.show()
    fig.savefig(filename+'.png')

    return
