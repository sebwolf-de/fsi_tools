"""Basis functions evaluation module.

.. moduleauthor:: Nicola Cavallini <andrew@invalid.com>

"""

import numpy as np

def tri_p1(x,y,eval_p):
    """
    Linear shape function on triangles, namely p1.

    Input:
    
    x : one dimensional array of triangle vertices x coords.\n
    y : one dimensional array of triangle vertices y coords.\n
    eval_p: (n,2) array of the n evaluation points. first 
            column indicates x-coord, second y-coord.\n
    
    Output:
    
    dx_phi : the three x-derivatives.\n
    dy_phi : the three y-derivatives.\n
    phi    : (n,3) array of the three shape funtions ath the n eval points.\n
    surf_e : the triangle area.\n
    
    Notice: all the quantities are computed on the current element
    
    """
    # =================================================|
    # given the triangle nodes                         |
    # i write this subroutine to evaluate              |
    # the basis (test) function drivatives and surface |
    # implied in the variational formulation           |
    # =================================================|
    # picewise linear basis function:                  |
    # -------------------------------------------------|
    # psi_j(X) = a_j * x + b_j * y +  c_j              |
    # -------------------------------------------------|
    # psi_1(X_1) = 1                                   |
    # psi_1(X_2) = 0                                   |
    # psi_1(X_3) = 0                                   |
    #                                                  |
    # similar for nodes 2,3                            |
    #                                                  |
    # a_j = der_psi_x_j                                |
    # b_j = der_psi_y_j                                |
    # j = 1...3                                        |
    # -------------------------------------------------|
    # follows the variables definition:                |
    # -------------------------------------------------|
    # ----------------------------------|
    # these variables mean i'm a coward |
    # ----------------------------------|
    #
    a_j = np.zeros((1,3), dtype = 'd')
    b_j = np.zeros((1,3), dtype = 'd')
    c_j = np.zeros((1,3), dtype = 'd')
    #
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    #
    y1 = y[0]
    y2 = y[1]
    y3 = y[2]
    #
    #
    surf_e = 1./2. * abs( x1*y3 - x1*y2 + x2*y1 - x2 * y3 + x3*y2 - x3*y1 )
    # ----------------------------------------------------------------------------|
    # pis_1 => psi_1(x1,y1) = 1                                                   |
    # ----------------------------------------------------------------------------|
    a_j[0,0] = -((-y2 + y3)/(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3))
    #
    b_j[0,0] = -((-x2 + x3)/(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3))
    #
    c_j[0,0] = -((-(x3*y2) + x2*y3)/(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3))
    #
    # ----------------------------------------------------------------------------|
    # pis_2 => psi_2(x2,y2) = 1                                                   |
    # ----------------------------------------------------------------------------|
    #
    x1 = x[1]
    x2 = x[2]
    x3 = x[0]
    #
    y1 = y[1]
    y2 = y[2]
    y3 = y[0]
    #
    a_j[0,1] = -((-y2 + y3)/(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3))
    #
    b_j[0,1] = -((-x2 + x3)/(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3))
    #
    c_j[0,1] = -((-(x3*y2) + x2*y3)/(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3))
    #
    # ----------------------------------------------------------------------------|
    # pis_3 => psi_3(x3,y3) = 1                                                   |
    # ----------------------------------------------------------------------------|
    #
    x1 = x[2]
    x2 = x[0]
    x3 = x[1]
    #
    y1 = y[2]
    y2 = y[0]
    y3 = y[1]
    #
    a_j[0,2] = -((-y2 + y3)/(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3))
    #
    b_j[0,2] = -((-x2 + x3)/(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3))
    #
    c_j[0,2] = -((-(x3*y2) + x2*y3)/(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3))
    #
    #
    dx_phi = a_j
    dy_phi = b_j
    #
    phi = np.zeros((0,3))
    #
    #print 'shape eval'
    #print a_j.shape
    #print b_j.shape
    #print c_j.shape
    #for point in eval_p:
    #    print point
    #    #f = a_j * point[0] + b_j*point[1] + c_j
    #    f =  point * np.array([a_j,b_j]) + c_j
    #    phi = np.vstack((phi,f))
    #
    #
    #print '---- start test---'
    #print phi
    coeffs = np.vstack([a_j,b_j,c_j])
    eval_p = np.array(eval_p)
    #print eval_p.shape
    #xq = eval_p[:,0]
    #yq = eval_p[:,1]
    nqp = eval_p.shape[0]
    coords = np.hstack([eval_p,
                        np.ones((nqp,1))])
    phi = np.dot(coords,coeffs)
    #print prova - phi
    #coords = np.hstack([eval_p,np.ones(nqp,1)])#np.hstack(reshape(eval_p[:,0],(,)))
    #print '---- end   test---'
    #
    #
    return dx_phi,dy_phi,phi,surf_e


def lin_p1(s,eval_p):
    """
    Linear shape function on lines, namely p1.

    Input:
    
    s : one dimensional array of line vertices coords.\n
    exal_p: (n,1) array of the n evaluation points.\n
    
    Output:
    
    a_j    : the two derivatives.\n
    phi    : (n,2) array of the two shape funtions ath the n eval points.\n
    omega  : the line length.\n
    
    """
    s0 = s[0]
    s1 = s[1]
    a_j = np.zeros((1,2))
    b_j = np.zeros((1,2))

    a_j[0,0] = 1/(s0-s1)
    a_j[0,1] = 1/(s1-s0)

    b_j[0,0] =  s1/(s1-s0)
    b_j[0,1] = -s0/(s1-s0)

    omega = s1-s0

    phi = np.zeros((0,2))

    for point in eval_p:
        f = a_j * point[0] + b_j
        phi = np.vstack((phi,f))
        
    return a_j, phi, omega
