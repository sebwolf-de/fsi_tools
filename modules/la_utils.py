import numpy as np
from scipy import sparse
from scipy import *
import scipy.sparse.linalg as sp_la
from numpy import linalg as la
import math as mth
import matplotlib.pyplot as plt

from preconditioner import NullPreconditioner

def linspace_wide(real,npts,width):
    # given the vector "real"
    # define equalliy spaced
    # spaced vector "r" going form
    # min(real)-width*(range) to max(real) + width*range
    # range = max - min
    r = np.linspace( np.amin(real) - width*(np.amax(real)-np.amin(real)), np.amax(real) + width*(np.amax(real)-np.amin(real)),npts)
    return r

def char_poly(real,imag,npts,width):
    # given real and imaginary parts of the
    # eigenvalues, it returns the characterstic polynomial
    # evaluated along pts points equally spaced going form
    # min(real)-width*(range) to max(real) + width*range
    # range = max - min, on the real axis. same thing on the
    # imaginary one.

    r = linspace_wide(real,npts,width)
    i = linspace_wide(imag,npts,width)
    r,i = np.meshgrid(r,i)

    eigs = np.vectorize(np.complex)(real,imag)
    z = np.vectorize(np.complex)(r,i)

    poly = np.ones(z.shape)
    poly = np.vectorize(np.complex)(poly)

    for ritz in eigs:
        poly *= z - ritz

    return r,i,poly


def real_img_eigs_parts(Nu):
    eigenvals, eigenvecs = la.eig(Nu.todense())
    imag_val = eigenvals.imag
    real_val = eigenvals.real
    return real_val, imag_val

def rand_sparse_matrix(mu,sigma,shape):
    # dense matrix in sparse format with random entries
    # used for exercies

    m = shape[0]
    n = shape[1]

    c_id = np.arange(0,n,dtype=np.int32)
    c_id = np.reshape(c_id,(1,n))

    c = np.dot(np.ones((m,1),dtype=np.int32),c_id)
    c = c.flatten()

    r_id = np.arange(0,m,dtype=np.int32)
    r_id = np.reshape(r_id,(m,1))

    r = np.dot(r_id,np.ones((1,n),dtype=np.int32))
    r = r.flatten()

    v = np.random.normal(mu,sigma,c.shape)

    A = sparse.coo_matrix((v, (r, c)), shape)

    return A

def arnoldi_iteration(A,n_iter,P=None):

    m = A.shape[0]

    if P is None:
        P = NullPreconditioner(A.shape)

    offsets = np.arange(-1,n_iter)

    n_el = n_iter
    for i in np.arange(1,n_iter):
        n_el += n_iter-i

    data = array([np.arange(1,n_iter+1)]).repeat(n_iter+1,axis=0)
    H = sparse.dia_matrix((data,offsets),(n_iter+1,n_iter),dtype=np.float64)

    c_id = np.arange(0,n_iter+1,dtype=np.int32)
    c_id = np.reshape(c_id,(1,n_iter+1))

    c = np.dot(np.ones((m,1),dtype=np.int32),c_id)
    c = c.flatten()

    r_id = np.arange(0,m,dtype=np.int32)
    r_id = np.reshape(r_id,(m,1))

    r = np.dot(r_id,np.ones((1,n_iter+1),dtype=np.int32))
    r = r.flatten()

    Q = sparse.coo_matrix((np.zeros(r.shape), (r,c)), (m,n_iter+1),dtype=np.float64)

    b = np.zeros((m,1),dtype=np.float64)
    b[0,0] = 1

    Q = Q.tolil()
    H = H.tolil()

    Q[:,0] = b/la.norm(b)

    for n in np.arange(0,n_iter):
        qn = Q[:,n].todense()
        qn = P.matvec(qn)
        v = A.dot(qn)
        for j in np.arange(0,n+1):
            qj = Q[:,j].todense()
            hjn = np.dot(qj.transpose(),v)[0,0]
            H[j,n] = hjn
            v = v - hjn*qj
        H[n+1,n] = la.norm(v,2)
        Q[:,n+1] = v/la.norm(v,2)
        ##print('arnoldi iter: ' + str(n) + ', of: ' + str(n_iter))

    ##print '----------------------'
    return H[:-1,:],Q[:,:-1]

def multiply_inverse(A,Bt,log=False):

    v = np.zeros((A.shape[0],0))
    Bt = Bt.tolil()

    # if log == True:
    #     line = 'multiply inverse, 0 % done.'
    #     ##print((line), end=' ')


    for i in range(0,Bt.shape[1]):
        f = Bt[:,i]
        s = sp_la.spsolve(A.tocsr(),f)
        s = np.reshape(s,(Bt.shape[0],1))
        v = np.hstack([v,s])
    #     if log == True:
    #         r = float(i)/float(Bt.shape[1])
    #         r = int(100*r)
    #         #print(('\r' * len(line)), end=' ')
    #         line = 'multiply inverse, ' + str(r)
    #         ##print output
    #         line += ' % done.'
    #         #print((line), end=' ')
    #
    # #if log == True:
    #     #print()

    return v

def delete_row_csr(mat, i):
    #http://stackoverflow.com/questions/13077527/is-there-a-numpy-delete-equivalent-for-sparse-matrices
    mat = mat.tocsr()
    n = mat.indptr[i+1] - mat.indptr[i]
    if n > 0:
        mat.data[mat.indptr[i]:-n] = mat.data[mat.indptr[i+1]:]
        mat.data = mat.data[:-n]
        mat.indices[mat.indptr[i]:-n] = mat.indices[mat.indptr[i+1]:]
        mat.indices = mat.indices[:-n]
    mat.indptr[i:-1] = mat.indptr[i+1:]
    mat.indptr[i:] -= n
    mat.indptr = mat.indptr[:-1]
    mat._shape = (mat._shape[0]-1, mat._shape[1])
    return mat

def delete_cln_csr(mat, i):
    mat = mat.transpose()
    mat = delete_row_csr(mat.tocsr(), i)
    mat = mat.transpose()
    return mat

def eye(n_dofs):
    rows = np.arange(0,n_dofs)
    uno = np.ones(n_dofs)
    A = sparse.coo_matrix((uno, (rows,rows)), shape=(n_dofs,n_dofs))
    return A

def delete_cln_list(X,lista):
    for i in lista:
        X = delete_cln_csr(X,i)
    return X

def delete_row_list(X,lista):
    for i in lista:
        X = delete_row_csr(X,i)
    return X

def l2_norm(M,g):
    l2_g = M.dot(g)
    l2_g = np.dot(l2_g.transpose(),g)
    l2_g = mth.sqrt(l2_g)
    return l2_g

def ismember(a, b):
    #http://stackoverflow.com/questions/7448554/replicating-the-indices-result-of-matlabs-ismember-function-in-numpy
    # tf = np.in1d(a,b) # for newer versions of numpy
    tf = np.array([i in b for i in a])
    u = np.unique(a[tf])
    index = np.array([(np.where(b == i))[0][-1] if t else 0 for i,t in zip(a,tf)])
    return tf, index


def get_connectivity_matrix(topo):
    dof_per_el = topo.shape[1]
    nels = topo.shape[0]
    ##print dof_per_el
    ##print nels
    rows = np.reshape(topo,(dof_per_el*nels))
    n_dofs = max(rows)+1
    cols = np.arange(0,nels)
    cols = np.reshape(cols,(cols.shape[0],1))
    cols = np.dot(cols,np.ones((1,3)))
    cols = np.reshape(cols,(dof_per_el*nels))
    cols = cols.astype(int)
    vals = np.ones(cols.shape)

    E = sparse.coo_matrix((vals,(rows,cols)),shape=(n_dofs,nels))
    return E

def get_connectivity_matrix_ieq(topo,ieq):
    dof_per_el = topo.shape[1]
    nels = topo.shape[0]
    ##print dof_per_el
    ##print nels
    rows = np.reshape(topo,(dof_per_el*nels))
    ##print 'rows = '
    ##print rows
    rows = ieq[rows]
    ##print rows
    n_dofs = max(rows)+1
    cols = np.arange(0,nels)
    cols = np.reshape(cols,(cols.shape[0],1))
    cols = np.dot(cols,np.ones((1,3)))
    cols = np.reshape(cols,(dof_per_el*nels))
    cols = cols.astype(int)
    vals = np.ones(cols.shape)

    E = sparse.coo_matrix((vals,(rows,cols)),shape=(n_dofs,nels))
    return E

def fluid_str_sparsity_pattern(topo_u,topo_s,ie_s,fluid_id):
    str_id = 0
    rows = np.zeros((1,0),dtype=int)
    cols = np.zeros((1,0),dtype=int)
    nel_f = topo_u.shape[0]
    ndofs_s = max(ie_s) +1
    for chunks in fluid_id:
        for elf in chunks:
            cln = ie_s[topo_s[str_id,:]]
            cln = np.reshape(cln,(1,3))
            cols = np.hstack([cols,cln])
            row  = elf*np.ones(cln.shape)
            rows = np.hstack([rows,row])
        str_id += 1
    values = np.ones(rows.shape)
    ##print nel_f
    ##print ndofs_s
    ##print max(cols[0])
    ##print max(rows[0])
    ETS = sparse.coo_matrix((values[0],(rows[0],cols[0])),shape=(nel_f,ndofs_s))
    EU = get_connectivity_matrix(topo_u)
    A = EU.dot(ETS)
    (rows,cols) = A.nonzero()
##    str_id = 0
##    rows = np.zeros((1,0),dtype=int)
##    cols = np.zeros((1,0),dtype=int)
##    for chunks in fluid_id:
##        ies_l = ie_s[topo_s[str_id,:]]
##        c = np.reshape(ies_l,(1,3))
##        ##print 'str eq ='
##        ##print ies_l
##        for fel in chunks:
##            r = np.reshape(topo_u[fel,:],(1,3))
##            rows = np.hstack([rows,r])
##            cols = np.hstack([cols,c])
##            ##print 'fluid eq ='
##            ##print topo_u[fel,:]
##        str_id += 1
##
##    entries = np.zeros(rows.shape,dtype = [('row',int),('cln',int)])
##    entries['row'] = rows
##    entries['cln'] = cols
##
##    entries = np.unique(entries)
##
##    rows = entries['row']
##    cols = entries['cln']
    return rows,cols

def get_sparsity_pattern_ieq(topo,ieq):

    E = get_connectivity_matrix_ieq(topo,ieq)

    A = E.transpose()
    A = E.dot(A)

    (rows,cols) = A.nonzero()
    return rows,cols

def get_sparsity_pattern(topo):

    E = get_connectivity_matrix(topo)

    A = E.transpose()
    A = E.dot(A)

    (rows,cols) = A.nonzero()
    return rows,cols

def get_mixed_sparsity_pattern(topo_u,topo_p):

    Eu = get_connectivity_matrix(topo_u)
    Ep = get_connectivity_matrix(topo_p)

    A = Ep.transpose()
    A = Eu.dot(A)

    (rows,cols) = A.nonzero()
    return rows,cols

def add_local_to_global_coo(rows,cols,values,
                            row,col,local_matrix):
    # first I concatenate local data
    [c,r] = np.meshgrid(row,col)
    r = np.concatenate(r)
    c = np.concatenate(c)
    local_matrix = np.concatenate(local_matrix)
    # I find the global entries corresponding to
    # the local data.
    tr = np.in1d(rows,row)
    tc = np.in1d(cols,col)
    vid = np.where((tr == True) & (tc == True))[0]
    # vid is the index of the entries in the values vector
    # that correspond to the local matrix entries.
    # Unfortunately the order in wich this values are
    # stored is not the samas in wich we build the local matrix.
    gr = rows[vid]
    gc = cols[vid]
    # we need to find a permutation index that map
    # the order of the local entries in the same order as
    # they are stored in the global data.
    entry = np.zeros(r.shape[0],dtype = [('row',int),('cln',int)])
    entry['row'] = r
    entry['cln'] = c
    # users have to control that entry and gerntry refer
    # to the same global matrix entries, but the order
    # in wich they are stored is not the same.
    gentry = np.zeros(r.shape[0],dtype = [('row',int),('cln',int)])
    gentry['row'] = gr
    gentry['cln'] = gc
    # entry and gentry are numpy stuctures. I use the analog
    # of matlab "ismember" to get the mapping form the
    # contruction order to the allocation order.
    (tf,index) = ismember(entry,gentry)
    ##print '----------'
    ##print vid
    ##print index
    ##print entry
    ##print gentry
    ##print '----------'
    # I finally permute the local matrix
    permuted_local_matrix = np.zeros(local_matrix.shape)
    permuted_local_matrix[index] = local_matrix
    # and add the local values to the global ones
    values[vid] = values[vid] + permuted_local_matrix
    return values

def add_local_to_global(B1,local_matrix,nd_u,nd_p):
    [r,c] = np.meshgrid(nd_u,nd_p)
    r = np.concatenate(r)
    c = np.concatenate(c)
    vals = np.concatenate(local_matrix.transpose())
    tmp = sparse.coo_matrix((vals, (r,c)), shape=B1.shape)
    B1 = B1 + tmp
    return B1

def get_diag_block(A,id_start,id_end):
    rows = A.shape[0]
    cols = A.shape[1]
    blk_size = id_end - id_start
    uno = np.ones((1,rows))
    uno = sparse.dia_matrix((uno,id_start), shape = (blk_size,rows))
    A = uno.dot(A)
    A = A.dot(uno.transpose())
    #ndofs = A.shape[0]
    #blk_size = id_end - id_start
    #uno = np.ones((1,blk_size))
    #uno = sparse.dia_matrix((uno,0), shape = (blk_size,ndofs))
    #A = uno.dot(A)
    #A = A.dot(uno.transpose())
    return A

def get_rows(A,id_start,id_end):
    rows = A.shape[0]
    cols = A.shape[1]
    blk_size = id_end - id_start
    uno = np.ones((1,rows))
    uno = sparse.dia_matrix((uno,id_start), shape = (blk_size,rows))
    A = uno.dot(A)
    return A

def get_cols(A,id_start,id_end):
    A = A.transpose()
    A = get_rows(A,id_start,id_end)
    A = A.transpose()
    return A

def set_diag(A,bc_id):
    ndofs = A.shape[0]
    #diago = A.diagonal()
    #Remove rows from A
    uno = np.ones((1,ndofs))
    uno[:,bc_id] = 0
    uno = sparse.dia_matrix((uno,0), shape = (ndofs,ndofs))
    A = uno.dot(A)
    #Set diagonals to 1
    uno = np.zeros((1,ndofs))
    uno[:,bc_id] = 1#diago[bc_id]
    uno = sparse.dia_matrix((uno,0), shape = (ndofs,ndofs))
    A = A+uno

    return A

def clear_rows(A,bc_id):
    ndofs = A.shape[0]
    uno = np.ones((1,ndofs))
    uno[:,bc_id] = 0
    uno = sparse.dia_matrix((uno,0), shape = (ndofs,ndofs))
    A = uno.dot(A)
    return A

def clear_cols(A,bc_id):
    ndofs = A.shape[1]
    uno = np.ones((1,ndofs))
    uno[:,bc_id] = 0
    uno = sparse.dia_matrix((uno,0), shape = (ndofs,ndofs))
    A = A.dot(uno)
    return A

def eliminate_col(BT1,i):
    BT1 = BT1.tocsr()
    BT1 = sparse.hstack([BT1[:,:i], BT1[:,i+1:]])
    BT1 = BT1.tocsr()
    return BT1

def eliminate_row(BT1,i):
    BT1 = BT1.tocsr()
    BT1 = sparse.vstack([BT1[:i, :], BT1[i+1:, :]])
    BT1 = BT1.tocsr()
    return BT1
