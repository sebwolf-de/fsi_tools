#import PetscBinaryIO
from scipy.sparse import *
from scipy import *
import numpy as np

def write_accuracy_table(filename,mesh_ref,p_error_l2,prex_l2_order,\
    u_error_l2,vel_l2_order):
    # write the accuracy table 
    # now it is prepared for L2 prex and L2 velocity
    # prex_l2_order should be dimenasion-1 with respect 
    # p_error_l2, same thing for velocity

    midrule = '\\cmidrule{1-1}'
    midrule += '\\cmidrule{3-5}'
    midrule += '\\cmidrule{7-9}'
    
    out_file = open(filename, "w")
    riga = '\\begin{tabular}{cc cc cc cc cc}\\toprule'
    out_file.write(riga)
    out_file.write("\n")
    riga = '$h_x$ & & $||p-p_h||_{L^2}$ & & $L^2$-rate & & $||\mathbf{u}-\mathbf{u}_h||_{L^2}$ & & $L^2$-rate & \\\\'
    riga += midrule
    out_file.write(riga)
    out_file.write("\n")
    riga = '$1/'+str(int(mesh_ref[0]))+'$& & '
    riga += "%7.4f" % (p_error_l2[0])
    riga += '&  &'
    riga += ' - '
    riga += '&  &'
    riga += "%7.4f" % (u_error_l2[0])
    riga += '&  &'
    riga += ' - '
    riga += '&  \\\\'
    riga += midrule
    out_file.write(riga)
    out_file.write("\n")
    for line in np.arange(1,mesh_ref.shape[0]):
        riga = '$1/'+str(int(mesh_ref[line]))+'$& & '
        riga += "%7.4f" % (p_error_l2[line])
        riga += '&  &'
        riga += "%7.2f" % (prex_l2_order[line-1])
        riga += '&  &'
        riga += "%7.4f" % (u_error_l2[line])
        riga += '&  &'
        riga += "%7.2f" % (vel_l2_order[line-1])
        riga += '&  \\\\'
        if line != mesh_ref.shape[0]-1:
            riga += midrule
        else:
            riga += '\\bottomrule'
        out_file.write(riga)
        out_file.write("\n")
    out_file.write('\\end{tabular}\n')
    out_file.close()
    return

def read_dealii_vector(filename):
    header = ''
    with open(filename,'r') as f:
        byte = f.read(1)
        #print byte
        #print f.read(1)
        #print f.read(1)
        while (byte!='['):        
            header += byte
            byte = f.read(1)
        #header = header[:-1]
        print(header)
        max_len = int(header)
        #brakets = f.read(1)
        val = np.fromfile(f, dtype=np.float64, count=max_len)
    return val

def read_dealii_sparsity_pattern(filename):
    header = ''
    with open(filename,'r') as f:
        byte = f.read(1)
        while (byte!=']'):        
            header += byte
            byte = f.read(1)
        header = header[1:]
        data = [int(n) for n in header.split()]
        [max_dim,\
        rows,\
        cols,\
        max_vec_len,\
        max_row_length,\
        compressed,\
        store_diagonal_first_in_row] = data    
        brakets = f.read(1)
        rowstart = np.fromfile(f, dtype=np.uint64, count=rows+1)
        brakets = f.read(2)
        columns =  np.fromfile(f, dtype=np.uint64, count=max_vec_len)
    return rows,cols,rowstart, columns
    
def read_dealii_matrix_values(filename):
    header = ''
    with open(filename,'r') as f:
        byte = f.read(1)
        while (byte!=']'):        
            header += byte
            byte = f.read(1)
        header = header[1:]
        max_len = int(header)
        brakets = f.read(1)
        val = np.fromfile(f, dtype=np.float64, count=max_len)
    return val
    
def read_dealii_matrix(filename):
    
    (rows,cols,rowstart, columns) = read_dealii_sparsity_pattern(filename+'_sp')
    val = read_dealii_matrix_values(filename+'_vl')
    
    matrix = csr_matrix((val, columns, rowstart), (rows, cols))
    
    return matrix


def petsc_to_scipy_mat(filename):
    io = PetscBinaryIO.PetscBinaryIO()
    objects = io.readBinaryFile(filename)
    shape = objects[0][0]
    data = objects[0][1][2]
    indices = objects[0][1][1]
    indptr = objects[0][1][0]
    A = csr_matrix((data, indices, indptr), shape)
    return A

def read_bc_map(filename):
    dof_name = filename + '_dof'
    f = open(dof_name,"r")
    bc_dof = np.fromfile(f,dtype=np.int32)
    f.close()

    val_name = filename + '_val'
    f = open(val_name,"r")
    bc_val = np.fromfile(f,dtype=np.float)
    f.close()
    return bc_dof, bc_val

def petsc_to_numpy_vec(filename):
    io = PetscBinaryIO.PetscBinaryIO()
    objects = io.readBinaryFile(filename)
    v = np.array(objects[0])
    return v
