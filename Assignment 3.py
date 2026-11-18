### Finite Volume Solver
### Will Abbott, 1013944706, Assignment 3

# Importing necessary functions from modules
import numpy as np
from scipy.sparse import coo_array
from scipy.sparse import csr_array
from scipy.sparse.linalg import bicg
from scipy.sparse.linalg import norm
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_array
from numpy.random import rand
from random import randint
import time
import math as math
import matplotlib.pyplot as plt

print("Program Start")
t_start = time.time()

L_x = 1
L_y = 1

N_x = 3
N_y = 3

delta_x = 1/N_x
delta_y = 1/N_y

gamma = 20
phi_B = 100
phi_b = 0

u_x = 2
u_y = 2

D = gamma/delta_x



