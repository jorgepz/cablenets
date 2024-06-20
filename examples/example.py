
# hanging cable

import sys
from os.path import dirname
sys.path.append(dirname('../src/'))

import numpy as np
from math import pi

from cablenets import solve, plot

# scalar parameters
L      = 2 
youngs = np.array([ 20 ])
areas  = np.array([ pi*.01**2 ])
nelems = 20

nodes  = np.zeros((nelems+1,3))
for i in range(nelems+1):
    nodes[i,:] = [i*L/nelems, 0.0, 0.0]

connec = np.zeros((nelems,4),dtype=int)
for i in range(nelems):
    connec[i,:] = [0, 0, i, i+1 ]

disp_mat = np.array([ [          0, 0   , 0, 0   ],
                      [   nelems/2, L*.5, L*.2, L*.1],
                      [     nelems, L   , 0, 0   ]])

fext_mat  = np.zeros((nelems-1,4))
for i in range(nelems-1):
    fext_mat[i,:] = [ i+1, 0, 0.0, -1.0e-3] # node fx fy fz

nodes_def, normal_forces = solve( nodes, connec, youngs, areas, disp_mat, fext_mat )

plot( nodes, connec, nodes_def, normal_forces )

