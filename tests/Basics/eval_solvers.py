#! /bin/env/python

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sp_la
import time

for k in range(0,7):
    n = 10**k

    A = sp.spdiags(np.random.rand(3, n), [-1, 0, 1], n, n).tocsc()
    B = sp.spdiags(np.random.rand(1, n), [0], n,n).tocsc()
    C = sp.hstack([sp.vstack([A, B]),
                   sp.vstack([B.transpose(), sp.csc_matrix((n,n))])],
                   "csc")
    b = np.random.rand(2*n)

    t0 = time.time()
    perm = sp.csgraph.reverse_cuthill_mckee(C, symmetric_mode = False)
    C_p = C[perm,:][:,perm]
    b_p = b[perm]

    M2 = sp_la.spilu(C_p)
    M_x = lambda x: M2.solve(x)
    M = sp_la.LinearOperator((2*n, 2*n), M_x)

    (x, it) = sp_la.bicgstab(C_p, b_p, M=M, tol=1e-8)
    t1 = time.time()
    x = x[np.argsort(perm)]

    print('-------')
    print('n = ' + str(n))
    print('res = ' + str(np.linalg.norm(C.dot(x) - b)))
    print('it = ' + str(it))
    print('t = ' + str(t1 - t0))
    print()
