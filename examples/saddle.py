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
from math import pi, sqrt

# scalar parameters
L       = 2        # length of the edge
nn      = 10       # number nodes per edge
E       = 200      # young modulus (GPa)
d       = .01      # cross-section diameter
gamma   = 785000*9.81/1e9 # density (kg/m3*m/s2/1e9: GN/m3)

nc = nn-1 # number cells per edge

nelems = nc*nc*2 + nc*4
nnodes = nn*nn

# node coordinates matrix
nodes  = np.zeros((nnodes,3))
for j in range(nn):
    for i in range(nn):
        ind_node = i + j*nn
        print("node: ", ind_node)
        nodes[ ind_node, :] = [i*L/nc, j*L/nc, 0.0]

# young and area vectors
youngs = np.array([ E ])
A      = pi*d**2/4
areas  = np.array([ A ])

# connectivity matrix
connec = np.zeros((nelems,4),dtype=int)
# cross cables
for j in range(nc):
    for i in range(nc):
        ind_cell = i + j*nc
        print("ind cell", ind_cell, " i ", i, " j ", j)
        connec[ind_cell*2  , : ] = [0, 0, j*nn+i  , (j+1)*nn+i+1 ]
        connec[ind_cell*2+1, : ] = [0, 0, j*nn+i+1, (j+1)*nn+i   ]
ind_cell = 2*nc*nc-1
for i in range(nc):
    print("i", i)
    ind_cell = ind_cell+1
    print("ind_cell", ind_cell)
    connec[ind_cell  , : ] = [0, 0, i  , i+1 ]
    ind_cell = ind_cell+1
    connec[ind_cell  , : ] = [0, 0, i+ nn*(nn-1)  , i+1+nn*(nn-1)]
    ind_cell = ind_cell+1
    connec[ind_cell  , : ] = [0, 0, i*nn  , (i+1)*nn ]
    ind_cell = ind_cell+1
    connec[ind_cell  , : ] = [0, 0, (i+1)*nn-1  , (i+2)*nn-1 ]

print(nodes)
print(connec)

pos_z = L*.5

disp_mat = np.array([ [          0,    0, 0, pos_z   ],
                      [       nn-1,    L, 0,     0   ],
                      [  nn*(nn-1),    0, L,     0   ],
                      [    nn*nn-1,    L, L, pos_z   ]])

print("disp_mat",  disp_mat)
lref = L * sqrt(2)/(2*nc)
fext_mat  = np.zeros((nnodes,4))
for i in range(nnodes):
    fext_mat[i,:] = [ i, 0, 0.0, -gamma*areas[0]*lref*2 ] # node fx fy fz
# fext_mat  = np.array([[1, 0.0, 0 , 0]])

nodes_def, normal_forces = solve( nodes, connec, youngs, areas, disp_mat, fext_mat )

plot( nodes, connec, nodes_def, normal_forces )

