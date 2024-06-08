
import sys
from os.path import dirname
sys.path.append(dirname('../src/'))

import numpy as np

from cablenets import solve, plot

L = 2
k = 20

nelems = 10

nodes  = np.zeros((nelems+1,3))
for i in range(nelems+1):
    print(i)
    nodes[i,:] = [i*L/nelems, 0.0, 0.0]

connec = np.zeros((nelems,2),dtype=int)
for i in range(nelems):
    print(i)
    connec[i,:] = [i, i+1 ]

print(nodes)

print(connec)

ks_vec = np.ones((nelems))* k

disp_mat = np.array([ [0, 0,0, 0], [nelems, L,0,0]])

fext_mat  = np.zeros((nelems-1,4))
for i in range(nelems-1):
    fext_mat[i,:] = [ i+1, 0, 0.0, -1.0] # node fx fy fz
print("f", fext_mat)

nodes_def, normal_forces = solve( nodes, connec, ks_vec, disp_mat, fext_mat )

print("fdfadas", np.size(nodes_def))

plot( nodes, connec, nodes_def, normal_forces )

