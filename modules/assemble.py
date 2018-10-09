"""Assemble Utilities.

.. moduleauthor:: Nicola Cavallini

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from shapely.geometry import Point

# nicola modules
import la_utils
import geom_utils as geom
import basis_func as basis

def local_p1_p1_mass_tri(x_l,y_l):
    eval_p = np.vstack([
        np.reshape(x_l,(1,3)),
        np.reshape(y_l,(1,3)),
        np.ones((1,3))])
    eval_p = np.dot(eval_p,bari_coords)
    eval_p = eval_p[0:2,:].transpose()
    (dx,dy,phi,omega) = basis.tri_p1(x_l,y_l,eval_p)
    lm = np.zeros((3,3))
    for i in np.arange(0,3):
        col = np.reshape(phi[i,:],(3,1))
        row = np.reshape(phi[i,:],(1,3))
        lm += omega/3*np.dot(col,row)
    return lm


def local_fluid_str_coupling(tri,xu_l,yu_l,s_l,t_l,tri_map):
    # tri = numpy array shape (3,2)
    # giving the triangle coords
    tri = np.hstack([tri,np.ones((3,1))]).transpose()
    #print '***'
    tri_ref = np.dot(tri_map,tri)
    area = np.linalg.det(tri_ref)/2
    #print tri_ref
    tri = np.dot(tri,bari_coords)
    eval_p = np.dot(np.eye(2,3),tri).transpose()
    (dx,dy,phi_u,omega) = basis.tri_p1(xu_l,yu_l,eval_p)
    eval_p = np.hstack([eval_p,np.ones((3,1))]).transpose()
    eval_p = np.dot(tri_map,eval_p)
    eval_p = eval_p[0:2,:].transpose()
    #print eval_p.shape
    (dx,dy,psi_s,omega) = basis.tri_p1(s_l,t_l,eval_p)
    lm = np.zeros((3,3))
    for i in np.arange(0,3):
        col = np.reshape(phi_u[i,:],(3,1))
        row = np.reshape(psi_s[i,:],(1,3))
        lm += area/3*np.dot(col,row)
    return lm

def local_mass_matrix_tri(tri,xu_l,yu_l,xs_l,ys_l):
    # tri = numpy array shape (3,2)
    # giving the triangle coords
    tri = np.hstack([tri,np.ones((3,1))]).transpose()
    #area = np.linalg.det(tri)/2
    tri = np.dot(tri,bari_coords)
    eval_p = np.dot(np.matlib.eye(2,3),tri).transpose()
    (dx,dy,phi_u,omega) = basis.tri_p1(xu_l,yu_l,eval_p)
    (dx,dy,psi_s,omega) = basis.tri_p1(xs_l,ys_l,eval_p)
    lm = np.zeros((3,3))
    for i in np.arange(0,3):
        col = np.reshape(phi_u[i,:],(3,1))
        row = np.reshape(psi_s[i,:],(1,3))
        lm += omega/3*np.dot(col,row)
    return lm

def u_v_p1_periodic(topo,x,y,ie):

    ndofs = max(ie)+1

    A = sparse.csr_matrix((ndofs,ndofs))

    for row in topo:
        x_l = x[row]
        y_l = y[row]
        entries = ie[row]
        #eval_p = np.zeros((0,2))
        eval_p = np.hstack([x_l.reshape(3,1),y_l.reshape(3,1)])
        local_matrix = local_p1_p1_mass_tri(x_l,y_l)
        #local_matrix = local_mass_matrix_tri(eval_p,x_l,y_l,x_l,y_l)
        #(phi_dx,phi_dy,phi,omega) = basis.tri_p1(x_l,y_l,eval_p)
        #local_matrix = omega/3. * np.eye(3)
        [r,c] = np.meshgrid(entries,entries)
        r = np.concatenate(r)
        c = np.concatenate(c)
        vals = np.concatenate(local_matrix)
        tmp = sparse.coo_matrix((vals, (r,c)), shape=(ndofs,ndofs))
        A = A+tmp
    return A

def u_v_p1(topo,x,y):

    ndofs = max(x.shape)

    A = sparse.csr_matrix((ndofs,ndofs))

    for row in topo:
        x_l = x[row]
        y_l = y[row]
        local_matrix = local_p1_p1_mass_tri(x_l,y_l)
        #eval_points = np.zeros((0,2))
        #(phi_dx,phi_dy,phi,omega) = basis.tri_p1(x_l,y_l,eval_points)
        #local_matrix = omega/3. * np.eye(3)
        [r,c] = np.meshgrid(row,row)
        r = np.concatenate(r)
        c = np.concatenate(c)
        vals = np.concatenate(local_matrix)
        tmp = sparse.coo_matrix((vals, (r,c)), shape=(ndofs,ndofs))
        A = A+tmp
    return A

def diagonal_mass_matrix(topo_s,s_lgr,t_lgr):
    ndofs = s_lgr.shape[0]
    entries = np.zeros((ndofs))

    for row in topo_s:
        x_l = s_lgr[row]
        y_l = t_lgr[row]
        eval_points = np.zeros((0,2))
        (phi_dx,phi_dy,phi,omega) = basis.tri_p1(x_l,y_l,eval_points)
        entries[row] += omega/3.
        #print omega

    rows = np.arange(0,s_lgr.shape[0])

    M = sparse.coo_matrix((entries, (rows,rows)), shape=(ndofs,ndofs))
    return M

def u_v_p1_inv_diag(topo,x,y):

    ndofs = max(x.shape)

    A = sparse.csr_matrix((ndofs,ndofs))

    for row in topo:
        x_l = x[row]
        y_l = y[row]
        #local_matrix = local_p1_p1_mass_tri(x_l,y_l)
        eval_points = np.zeros((0,2))
        (phi_dx,phi_dy,phi,omega) = basis.tri_p1(x_l,y_l,eval_points)
        local_matrix = omega/3. * np.eye(3)
        [r,c] = np.meshgrid(row,row)
        r = np.concatenate(r)
        c = np.concatenate(c)
        vals = np.concatenate(local_matrix)
        tmp = sparse.coo_matrix((vals, (r,c)), shape=(ndofs,ndofs))
        A = A+tmp

    vals = A.data
    vals = np.power(vals,-1)


    (rows,cols) =  A.nonzero()
    A = sparse.coo_matrix((vals, (rows,cols)), shape=(ndofs,ndofs))
    return A

def u_v_p1_1d_inv_diag(topo,x):

    ndofs = max(x.shape)

    A = sparse.csr_matrix((ndofs,ndofs))

    for row in topo:
        x_l = x[row]
        omega = np.fabs(x_l[1]-x_l[0])
        #local_matrix = local_p1_p1_mass_tri(x_l,y_l)
        local_matrix = omega/2. * np.eye(2)
        [r,c] = np.meshgrid(row,row)
        r = np.concatenate(r)
        c = np.concatenate(c)
        vals = np.concatenate(local_matrix)
        tmp = sparse.coo_matrix((vals, (r,c)), shape=(ndofs,ndofs))
        A = A+tmp

    vals = A.data
    vals = np.power(vals,-1)

    (rows,cols) =  A.nonzero()
    A = sparse.coo_matrix((vals, (rows,cols)), shape=(ndofs,ndofs))
    return A


##def gradu_gradv_p1(topo,x,y):
##
##    ndofs = max(x.shape)
##
##    A = sparse.csr_matrix((ndofs,ndofs))
##
##    for row in topo:
##        x_l = x[row]
##        y_l = y[row]
##        eval_points = np.zeros((0,2))
##        (phi_dx,phi_dy,phi,omega) = basis.tri_p1(x_l,y_l,eval_points)
##        dx_j = phi_dx
##        dx_i = phi_dx.transpose()
##        dy_j = phi_dy
##        dy_i = phi_dy.transpose()
##        local_matrix = omega*(np.dot(dx_i,dx_j)+np.dot(dy_i,dy_j))
##        [r,c] = np.meshgrid(row,row)
##        r = np.concatenate(r)
##        c = np.concatenate(c)
##        vals = np.concatenate(local_matrix)
##        tmp = sparse.coo_matrix((vals, (r,c)), shape=(ndofs,ndofs))
##        A = A+tmp
##    return A

def gradu_gradv_p1_ieq(topo,x,y,ieq):

    ndofs = max(ieq)+1

    (rows,cols)= la_utils.get_sparsity_pattern_ieq(topo,ieq)

    values = np.zeros(rows.shape)

    for row in topo:
        x_l = x[row]
        y_l = y[row]
        eval_points = np.zeros((0,2))
        (phi_dx,phi_dy,phi,omega) = basis.tri_p1(x_l,y_l,eval_points)
        dx_j = phi_dx
        dx_i = phi_dx.transpose()
        dy_j = phi_dy
        dy_i = phi_dy.transpose()
        entries = ieq[row]
        #print x_l
        #print y_l
        #print row
        #print entries
        local_matrix = omega*(np.dot(dx_i,dx_j)+np.dot(dy_i,dy_j))
        values = la_utils.add_local_to_global_coo(rows,cols,values,
                            entries,entries,local_matrix)

    A = sparse.coo_matrix((values,(rows,cols)),shape=(ndofs,ndofs))
    A.tocsr()

    return A

def gradu_gradv_p1(topo,x,y):
    """
    Assembling the Laplace operator. The function name resambles the
    operator gradtient of the trial functionctions, multiplied the gradient of
    the test functions. Assuming :math:`P_1` elements on the :math:`K` trinagle we have:

    .. math::
       \int_K \mathrm{grad}(u)_j \cdot \mathrm{grad}(v)_i = \mathrm{Area}(K)\ \mathrm{grad}(u_j) \cdot \mathrm{grad}(v_i)

    In the code snippet, we can see that, by default, derivatives are represented a 1x3 row.
    The ``dx_i`` and ``dy_i`` components are transpose. So we only need the matrix product ``np.dot``
    to have the local stffness matrix. The values for thederivatives are obtained
    form the ``tri_p1`` function in the `basis_func` module, `check its documentation`_.

    .. _check its documentation: ./basis_func.html

    .. code:: python

        dx_j = phi_dx
        dx_i = phi_dx.transpose()
        dy_j = phi_dy
        dy_i = phi_dy.transpose()
        local_matrix = omega*(np.dot(dx_i,dx_j)+np.dot(dy_i,dy_j))


    Input:

    ``x``, ``y``, ``topo`` : the nodes coordinates and the connectivity.

    Output:

    ``A`` : The spasre matrix A representing the discretized operator.\n

    """
    ndofs = max(x.shape)

    (rows,cols)= la_utils.get_sparsity_pattern(topo)

    values = np.zeros(rows.shape)

    for row in topo:
        x_l = x[row]
        y_l = y[row]
        eval_points = np.zeros((0,2))

        (phi_dx,phi_dy,phi,omega) = basis.tri_p1(x_l,y_l,eval_points)
        dx_j = phi_dx
        dx_i = phi_dx.transpose()
        dy_j = phi_dy
        dy_i = phi_dy.transpose()
        local_matrix = omega*(np.dot(dx_i,dx_j)+np.dot(dy_i,dy_j))
        values = la_utils.add_local_to_global_coo(rows,cols,values,
                            row,row,local_matrix)

    A = sparse.coo_matrix((values,(rows,cols)),shape=(ndofs,ndofs))
    #plt.spy(A)
    #plt.show()
    A.tocsr()

    return A

def u_gradv_w_p1(topo, x, y, u_x, u_y):

        ndofs = max(x.shape)

        A11 = sparse.csr_matrix((ndofs,ndofs))
        A12 = sparse.csr_matrix((ndofs,ndofs))
        A21 = sparse.csr_matrix((ndofs,ndofs))
        A22 = sparse.csr_matrix((ndofs,ndofs))

        for row in topo:
            x_l = x[row]
            y_l = y[row]
            local_matrix = np.zeros((3,3))
            local_mass_matrix = local_p1_p1_mass_tri(x_l,y_l)

            eval_points = np.zeros( (3,2) )
            eval_points[:,0] = x_l.transpose()
            eval_points[:,1] = y_l.transpose()

            (w_dx,w_dy,w_l,omega_w) = basis.tri_p1(x_l,y_l,eval_points)
            (v_dx,v_dy,v_l,omega_v) = basis.tri_p1(x_l,y_l,eval_points)


            #int_w_omega = np.zeros((1,3))
            #for k in range(0,3):
            #    int_w_omega[0,k] += omega_w/3 * sum(w_l[:,k])

            #print "-----------"
            #print local_mass_matrix
            #print u_x[row]
            #print v_dx
            #print w_dx

            local_matrix = np.dot(u_x[row].transpose(), local_mass_matrix)
            local_matrix = np.dot(v_dx.transpose(), local_matrix)
            A11 = la_utils.add_local_to_global(A11,local_matrix,row,row)

            local_matrix = np.dot(u_y[row].transpose(), local_mass_matrix)
            local_matrix = np.dot(v_dy.transpose(), local_matrix)
            A12 = la_utils.add_local_to_global(A12,local_matrix,row,row)

            local_matrix = np.dot(u_x[row].transpose(), local_mass_matrix)
            local_matrix = np.dot(v_dx.transpose(), local_matrix)
            A21 = la_utils.add_local_to_global(A21,local_matrix,row,row)

            local_matrix = np.dot(u_y[row].transpose(), local_mass_matrix)
            local_matrix = np.dot(v_dy.transpose(), local_matrix)
            A22 = la_utils.add_local_to_global(A22,local_matrix,row,row)

        #plt.spy(A11)
        #plt.show()
        A11 = 0.5*(A11-A11.transpose())
        A12 = 0.5*(A12-A21.transpose())
        A21 = 0.5*(A12-A21.transpose())
        A22 = 0.5*(A22-A22.transpose())

        return A11, A12, A21, A22

def divu_p_p1_iso_p2_p1(topo_p,x_p,y_p,
           topo_u,x_u,y_u,c2f):

    ndofs_u = max(x_u.shape)
    ndofs_p = max(x_p.shape)

    B1 = sparse.csr_matrix((ndofs_u,ndofs_p))
    B2 = sparse.csr_matrix((ndofs_u,ndofs_p))

    el_id = 0
    for nd_p in topo_p:
        #print("====================")
        #print nd_p
        fine_els = c2f[el_id,:]
        x_lp = x_p[nd_p]
        y_lp = y_p[nd_p]
        local_matrix = np.zeros((3,3))
        for el in fine_els:
            nd_u = topo_u[el,:]
            x_lu = x_u[nd_u]
            y_lu = y_u[nd_u]
            eval_p = np.zeros( (3,2) )
            eval_p[:,0] = x_lu.transpose()
            eval_p[:,1] = y_lu.transpose()
            (u_dx,u_dy,u,omega_u) = basis.tri_p1(x_lu,y_lu,eval_p)
            (p_dx,p_dy,p,omega_p) = basis.tri_p1(x_lp,y_lp,eval_p)
            #print("--------------------")
            #print nd_u
            #print u_dx
            int_q_omega = np.zeros((1,3))
            #print sum(p[:,0])
            for k in range(0,3):
                int_q_omega[0,k] += omega_u/3 * sum(p[:,k])
            #print 12*int_q_omega
            local_matrix = np.dot(u_dx.transpose(),
                                  int_q_omega)
            #print local_matrix
            B1 = la_utils.add_local_to_global(B1,local_matrix,nd_u,nd_p)
            local_matrix = np.dot(u_dy.transpose(),
                                  int_q_omega)
            B2 = la_utils.add_local_to_global(B2,local_matrix,nd_u,nd_p)
            #print B1.todense()
        el_id += 1

    return B1, B2

def divu_p_p1_iso_p2_p1p0(topo_p,x_p,y_p,
           topo_u,x_u,y_u,c2f):

    ndofs_u = max(x_u.shape)
    ndofs_p = max(x_p.shape) + topo_p.shape[0]

    B1 = sparse.csr_matrix((ndofs_u,ndofs_p))
    B2 = sparse.csr_matrix((ndofs_u,ndofs_p))

    el_id = 0
    for nd_p in topo_p:
        #print("====================")
        #print nd_p
        fine_els = c2f[el_id,:]
        x_lp = x_p[nd_p[0:3]]
        y_lp = y_p[nd_p[0:3]]
        local_matrix = np.zeros((3,4))
        for el in fine_els:
            nd_u = topo_u[el,:]
            x_lu = x_u[nd_u]
            y_lu = y_u[nd_u]
            eval_p = np.zeros( (3,2) )
            eval_p[:,0] = x_lu.transpose()
            eval_p[:,1] = y_lu.transpose()
            (u_dx,u_dy,u,omega_u) = basis.tri_p1(x_lu,y_lu,eval_p)
            (p_dx,p_dy,p,omega_p) = basis.tri_p1(x_lp,y_lp,eval_p)
            #print("--------------------")
            #print nd_u
            #print u_dx
            int_q_omega = np.zeros((1,4))
            #print sum(p[:,0])
            for k in range(0,3):
                int_q_omega[0,k] += omega_u/3 * sum(p[:,k])
            int_q_omega[0,3] = omega_u
            #print 12*int_q_omega
            local_matrix = np.dot(u_dx.transpose(),
                                  int_q_omega)
            #print local_matrix
            B1 = la_utils.add_local_to_global(B1,local_matrix,nd_u,nd_p)
            local_matrix = np.dot(u_dy.transpose(),
                                  int_q_omega)
            B2 = la_utils.add_local_to_global(B2,local_matrix,nd_u,nd_p)
            #print B1.todense()
        el_id += 1

    return B1, B2

def ibm_force(XY_str,s_lgr,topo_u,x_u,y_u,point_in_tri):
    node_s = XY_str.shape[0]
    force_x = np.zeros((x_u.shape[0],1))
    force_y = np.zeros((x_u.shape[0],1))

    Xs = XY_str[:,0]
    Ys = XY_str[:,1]
    for cnt in range(0,node_s):
        if cnt == 0:
            backward = node_s-1
            middle = cnt
            forward = cnt+1
            ds = s_lgr[forward]-s_lgr[middle]
        if cnt == node_s-1:
            backward = cnt-1
            middle = cnt
            forward = 0
            ds = s_lgr[middle]-s_lgr[backward]
        else:
            backward = cnt-1
            middle = cnt
            forward = cnt+1
            ds = s_lgr[forward]-s_lgr[middle]
        Ds_X1 = (Xs[forward] - Xs[middle])/ds
        Ds_X0 = (Xs[middle] - Xs[backward])/ds
        Ds_Y1 = (Ys[forward] - Ys[middle])/ds
        Ds_Y0 = (Ys[middle] - Ys[backward])/ds
        stiff_x = Ds_X1-Ds_X0
        stiff_y = Ds_Y1-Ds_Y0
        nds = topo_u[point_in_tri[cnt]]
        x_l = x_u[nds]
        y_l = y_u[nds]
        eval_p = np.zeros((1,2))
        eval_p[0,:] = XY_str[cnt,:]
        (phi_dx,phi_dy,phi,omega) = basis.tri_p1(x_l,y_l,eval_p)
        force_x[nds] += stiff_x * phi.transpose()
        force_y[nds] += stiff_y * phi.transpose()
    return force_x, force_y

def u_s_p1(topo_u,x_u,y_u,
          topo_s,x_s,y_s,s_lgr,ieq_s,
          str_segments,fluid_id):

   r = max(ieq_s)+1
   #print r

   GT = sparse.csr_matrix((x_u.shape[0],r[0]))

   str_iel = 0
   for str_el in str_segments:
       s_dofs = ieq_s[topo_s[str_iel,:]]
       s_l = s_lgr[topo_s[str_iel,:]]
       el_list = fluid_id[str_iel]
       #print '================='
       #print s_l
       #print el_list
       iel = 0
       for segment in str_el:
           f_id = el_list[iel]
           u_dofs = topo_u[f_id,:]
           l = list(segment.coords)
           sp = geom.get_reference_coords(topo_s,x_s,y_s,s_lgr,str_iel,l)
           (ds_psi,psi,omega) = basis.lin_p1(s_l,sp)
           x_ul = x_u[topo_u[f_id,:]]
           y_ul = y_u[topo_u[f_id,:]]
           p0 = Point(l[0])
           p1 = Point(l[1])
           eval_p = np.zeros((2,2))
           eval_p[0,0] = p0.x
           eval_p[0,1] = p0.y
           eval_p[1,0] = p1.x
           eval_p[1,1] = p1.y
           (phi_dx,phi_dy,phi,omega) = basis.tri_p1(x_ul,y_ul,eval_p)
           #print '-----------------'
           #print f_id
           #print s_dofs
           #print u_dofs
           #print sp
           #print psi
           for i in range(0,2):#loop over quadrature points
               cln = np.zeros((3,1))
               row = np.zeros((1,2))
               cln[:,0] = phi[i,:]
               row[0,:] = psi[i,:]
               local_matrix = .5 * segment.length * np.dot(cln,row)
               GT = la_utils.add_local_to_global(GT,local_matrix,u_dofs,s_dofs)
               #print '***'
               #print row
               #print cln
               #print local_matrix*8
               #print '***'
           #
           #print segment.length
           #print x_ul
           #print f_id
           #print psi
           #print phi
           iel += 1
           #
       str_iel+=1
       #break

   return GT

def u_s_p1_thick(x_u,y_u,topo_u,
                s_lgr,t_lgr,
                x_str,y_str,topo_s,ie_s,
                str_segments,fluid_id):

   #(rows,cols) = la_utils.fluid_str_sparsity_pattern(
   #topo_u,topo_s,ie_s,fluid_id)

   #print rows
   #print cols

   righe = x_u.shape[0]
   colonne = max(ie_s)+1

   GT = sparse.csr_matrix((righe,colonne))

   #values = np.zeros(rows.shape)

   str_id = 0
   for chunks in str_segments:
       #print '======================'
       #print 'els = ' + str(str_id)
       nds = topo_s[str_id,:]
       xs_l = x_str[nds]
       ys_l = y_str[nds]
       s_l = s_lgr[nds]
       t_l = t_lgr[nds]
       tri_map = geom.tri_lin_map(xs_l,ys_l,s_l,t_l)
       #print s_l
       #print t_l
       ies_l = ie_s[nds]
       chunk_id = 0
       for poly in chunks:
           if poly.area>1e-10:
               elf = fluid_id[str_id][chunk_id]
               #print 'elf = ' + str(elf)
               ndf = topo_u[elf,:]
               xu_l = x_u[ndf]
               yu_l = y_u[ndf]
               #print xu_l
               #print yu_l
               triangles = geom.triangulate(poly)
               local_matrix = np.zeros((3,3))
               for tri in triangles:
                   #print '----------------------'
                   tmp = np.array(list(tri.exterior.coords)[0:3])
                   #print tmp
                   lm = local_fluid_str_coupling(tmp,xu_l,yu_l,s_l,t_l,tri_map)
                   #lm = local_mass_matrix_tri(tmp,xu_l,yu_l,xs_l,ys_l)
                   local_matrix += lm
               #print local_matrix*24*6*6
               GT = la_utils.add_local_to_global(GT,local_matrix,ndf,ies_l)
               #values = la_utils.add_local_to_global_coo(rows,cols,values,
               #                ndf,ies_l,local_matrix)
                   #break
                   #print eval_p
                   #values = la_utils.add_local_to_global_coo(rows,cols,values,
                   #            ndf,ies_l,lm)
                   #print 'triangle'
               #for tri in triangles:
               #    print len(tri.exterior.coords)
               #print chunk_id
           chunk_id+=1
           #break
           #break
       str_id += 1
       #break
   #MT = sparse.coo_matrix((values,(rows,cols)),shape=(righe,colonne))
   GT.tocsr()
   #print GT
   #print MT
   return GT

def u_v_lin_p1(topo_s,s_lgr,ieq_s):
    r = max(ieq_s)+1
    M = sparse.csr_matrix((r[0],r[0]))

    for s_nds in topo_s:
        s_dofs = ieq_s[s_nds]
        s_l = s_lgr[s_nds]
        ds = s_l[1]-s_l[0]
        eval_p = np.zeros((2,1))
        #print eval_p
        #print s_l
        eval_p[:,0] = s_l
        (ds_psi,psi,omega) = basis.lin_p1(s_l,eval_p)
        for i in range(0,2):#loop over quadrature points
            cln = np.zeros((2,1))
            row = np.zeros((1,2))
            cln[:,0] = psi[i,:]
            row[0,:] = psi[i,:]
            local_matrix = .5 * ds * np.dot(cln,row)
            M = la_utils.add_local_to_global(M,local_matrix,s_dofs,s_dofs)
            #print '-------------'
            #print local_matrix
            #
        #print local
        #print psi
        #break
    return M

def gradu_gradv_lin_p1(topo_s,s_lgr,ieq_s):
    r = max(ieq_s)+1
    FX = sparse.csr_matrix((r[0],r[0]))

    #ds = s_lgr[1]-s_lgr[0]
    for s_nds in topo_s:
        s_dofs = ieq_s[s_nds]
        s_l = s_lgr[s_nds]
        ds = s_l[1]-s_l[0]
        eval_p = np.zeros((2,1))
        eval_p[:,0] = s_l
        (ds_psi,psi,omega) = basis.lin_p1(s_l,eval_p)
        local_matrix = ds*np.dot(
            ds_psi.transpose(),
            ds_psi)
        FX = la_utils.add_local_to_global(FX,local_matrix,s_dofs,s_dofs)
    return FX

bari_coords = np.array([[2./3.,1./6.,1./6.],
                        [1./6.,2./3.,1./6.],
                        [1./6.,1./6.,2./3.]])
