
import sys
from os.path import dirname
sys.path.append(dirname('../src/'))

import numpy as np

from cablenets import solve, plot

L = 2
k = 20

nelems = 100

q=1

nodes  = np.zeros((nelems+1,3))
for i in range(nelems+1):
    nodes[i,:] = [i*L/nelems, 0.0, 0.0]

connec = np.zeros((nelems,2),dtype=int)
for i in range(nelems):
    connec[i,:] = [i, i+1 ]

ks_vec = np.ones((nelems))* k

disp_mat = np.array([ [          0, 0   , 0, 0   ],
                      [   nelems/2, L*.5, 0, L*.05],
                      [     nelems, L   , 0, 0   ]])

fext_mat  = np.zeros((nelems-1,4))
for i in range(nelems-1):
    fext_mat[i,:] = [ i+1, 0, 0.0, -0.01] # node fx fy fz

nodes_def, normal_forces = solve( nodes, connec, ks_vec, disp_mat, fext_mat )

plot( nodes, connec, nodes_def, normal_forces )

