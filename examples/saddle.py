# ------------------------
# saddle cable net example
# ------------------------
# square-shaped cable net with two high and two low points

# add path of code
import sys
from os.path import dirname
sys.path.append(dirname('../src/'))

# import modules
from cablenets import solve, plot
import numpy as np

# scalar parameters
L   = 2 # length of the edge
nn  = 3 # number nodes per edge
E   = 1 # young modulus
A   = 1 # cross-section area

nc = nnodes-1 # number cells per edge

nelems = nc*nc*4 + nc*2
nnodes = nn*nn

# node coordinates matrix
nodes  = np.zeros((nnodes,3))
for i in range(nn):
    for j in range(nn):
        ind_node = i*nn
        nodes[ ind_node, :] = [i*L/nelems, 0.0, 0.0]

# young and area vectors
youngs = np.array([ E ])
areas  = np.array([ A ])

# connectivity matrix
connec = np.zeros((nelems,4),dtype=int)
for i in range(nelems):
    connec[i,:] = [0, 0, i, i+1 ]

disp_mat = np.array([ [          0, 0   , 0, 0   ],
                      [   nelems/2, L*.5, 0, L*.1],
                      [     nelems, L   , 0, 0   ]])

fext_mat  = np.zeros((nelems-1,4))
for i in range(nelems-1):
    fext_mat[i,:] = [ i+1, 0, 0.0, -1.0e-3] # node fx fy fz

nodes_def, normal_forces = solve( nodes, connec, youngs, areas, disp_mat, fext_mat )

plot( nodes, connec, nodes_def, normal_forces )

