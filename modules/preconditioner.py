from scipy.sparse.linalg import LinearOperator
import scipy.sparse.linalg as sp_la
import numpy as np

class IdealPreconditioner(LinearOperator):
    def mult(self,v):
        v = sp_la.spsolve(self.A,v)
        return v
    def __init__(self,A):
        LinearOperator.__init__(self,A.shape,self.mult)
        self.A = A
        return


class NullPreconditioner(LinearOperator):
    def mult(self,v):
        return v
    def __init__(self,shape):
        LinearOperator.__init__(self,shape,self.mult)
        return


class BlockPreconditioner(LinearOperator):
    def mult(self,v):

        size = 0
        for ndof in self.size_blk:
            size+=ndof
        
        assert size == self.shape[0], \
               "The size of the global matrix %r,\
do not match the size of the Preconditioner %r, check you\
set the correct blocks." % (self.shape[0],size)
        
        #print v.shape
        #print self.shape
        v = np.reshape(v,self.shape[0])
        
        rhs_list = []
        sol_list = []
        offset = 0
        for i in range(0,self.n_blk):
            rhs_list.append(v[offset:offset+self.size_blk[i]])
            sol_list.append([])
            offset+=self.size_blk[i]
        
        r_id = 0
        for row in self.pc:
            rhs_id = self.n_blk-1-r_id
            #print 'riga'
            #print rhs_id
            f = rhs_list[rhs_id]
            #print rhs_id
            c_id = 0
            for cln in row[:-1]:
                #print 'colonna'
                sol_id = self.n_blk-1-c_id
                if type(cln) != bool:
                    #print 'sol_id = '+ str(sol_id)
                    f -= cln.dot(sol_list[sol_id])
                    #f -= cln.dot(sol_list[self.n_blk-1-r_id])
                c_id +=1
            mat = self.pc[r_id][r_id].tocsr()
            #print mat
            #print f.shape[0]
            sol = sp_la.spsolve(mat,f)
            sol_list[rhs_id] = sol
            r_id +=1

        v = np.zeros((self.shape[0]))
        offset = 0
        for i in range(0,self.n_blk):
            v[offset:offset+self.size_blk[i]] = sol_list[i]
            offset+=self.size_blk[i]

        self.n_iter +=1
        res = self.evaluate_residual(v)
        print 'iter = '+str(self.n_iter)+', residual = '+str(res)
        return v
    def __init__(self,shape,n_blk):
        self.n_iter = 0
        self.shape = shape
        LinearOperator.__init__(self,self.shape,self.mult)
        self.pc = []
        self.size_blk = []
        self.n_blk = n_blk

        assert shape[0] == shape[1], \
               "Prec Shape Must be squared, %r rows  not equal to %r cols" % (shape[0],shape[1])
        
        for i in range(0,n_blk):
            row = (i+1)*[False]
            #for j in range(0,i+1):
            #    row.append([])
            self.pc.append(row)
        
        for i in range(0,n_blk):
            self.size_blk.append([])
        
        return
    def set_block_ij(self,i,j,matrix):
        if i == j :
            self.size_blk[i] = matrix.shape[0]
        i = self.n_blk - i - 1
        j = self.n_blk - j - 1
        self.pc[i][j] = matrix
        return
    def print_info(self):
        print self.pc
        print self.size_blk
        return
    def set_operators(self,mat,f):
        self.mat_g = mat
        self.f_g = f
        return
    def evaluate_residual(self,v):
        res = self.mat_g.dot(v)-self.f_g
        res = np.sqrt(np.sum(np.power(res,2)))
        return res

            
    

class Preconditioner(LinearOperator):
    def print_report(self,flag):
        print '----------------------------------------'
        print 'solution flag  = '+str(flag)
        print 'iterations number = '+str(self.n_iter)
        print '----------------------------------------'
        return
    def mv(self,v):
        # a do nothing preconditioner
        self.n_iter +=1
        print 'iter = ' +str(self.n_iter)
        return v
    def set_block_00(self,A):
        self.A = A.tocsr()
        self.size_00 = A.shape[0]
        return
    def set_block_11(self,Mp):
        self.Mp = Mp.tocsr()
        self.size_11 = Mp.shape[0]
        return
    def get_niter(self):
        return self.n_iter

class Ideal(Preconditioner):
    def full(self,v):
        v = sp_la.spsolve(self.A,v)
        self.n_iter +=1
        print 'iter #' +str(self.n_iter)
        return v
    def __init__(self,shape):
        LinearOperator.__init__(self,shape,self.full)
        self.n_iter = 0

class StkPreconditioner(Preconditioner):
    def block_diagonal(self,v):
        n_prex = self.size_11
        n_vel = self.size_00
        f_vel = v[0:n_vel]
        f_prex = v[n_vel:]
        s_vel = sp_la.spsolve(self.A,f_vel)
        s_prex = sp_la.spsolve(self.Mp,f_prex)
        v[0:n_vel] = s_vel
        v[n_vel:] = s_prex
        self.n_iter +=1
        res = self.evaluate_residual(v)
        print 'iter = '+str(self.n_iter)+', residual = '+str(res)
        return v
    def __init__(self,shape):
        LinearOperator.__init__(self,shape,self.block_diagonal)
        self.n_iter = 0
    def set_operators(self,mat,f):
        self.mat_g = mat
        self.f_g = f
        return
    def evaluate_residual(self,v):
        res = self.mat_g.dot(v)-self.f_g
        res = np.sqrt(np.sum(np.power(res,2)))
        return res


class NavierPreconditioner(Preconditioner):
    def upper_tri(self,v):
        #print "ndo cazzo sto"
        n_prex = self.size_11
        n_vel = self.size_00
        f_vel = v[0:n_vel]
        f_prex = v[n_vel:]
        f_prex = np.reshape(f_prex,(self.Mp.shape[0],1))
        s_prex = sp_la.spsolve(self.Mp,f_prex)
        print 'pressione risolta'
        print s_prex.shape
        s_prex = np.reshape(s_prex,(self.Mp.shape[0],1))
        f_vel = f_vel - self.Bt.dot(s_prex)
        print f_vel.shape
        s_vel = sp_la.spsolve(self.A,f_vel)
        s_vel = np.reshape(s_vel,(self.A.shape[0],1))
        #
        v[0:n_vel] = s_vel
        v[n_vel:] = s_prex
        print 'iter = ' +str(self.n_iter)
        self.n_iter +=1
        return v
    def __init__(self,shape):
        LinearOperator.__init__(self,shape,self.upper_tri)
        self.n_iter = 0
    def set_block_01(self,Bt):
        self.Bt = Bt.tocsr()
        return
