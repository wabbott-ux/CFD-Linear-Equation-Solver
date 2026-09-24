import numpy as np
from scipy.sparse import csc_array
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_array
from numpy.random import rand


A = lil_array((1000000, 1000000))
A[0, :100] = rand(100)
A.setdiag(rand(100))

A = A.tocsr()
b = rand(1000000)
x = spsolve(A, b)

        
