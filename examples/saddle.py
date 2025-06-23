# ------------------------
# saddle cable net example
# ------------------------
# square-shaped cable net with two high and two low points

# add path of code
import sys
from os.path import dirname
sys.path.append(dirname('../src/'))

# import modules
import cablenets as cn
import numpy as np
from math import pi, sqrt

# scalar parameters
L       = 2        # length of the edge
nn      = 12       # number nodes per edge
E       = 200      # young modulus (GPa)
d       = .01      # cross-section diameter
gamma   = 785000*9.81/1e9 # density (kg/m3*m/s2/1e9: GN/m3)

# create materials and areas vectors
linear_material = cn.Material(0, 'linear', E )
areas  = np.array([ pi*d**2/4 ])

nc = nn-1 # number cells per edge

nelems = nc*nc*2 + nc*4
nnodes = nn*nn

# node coordinates matrix
nodes  = np.zeros((nnodes,3))
for j in range(nn):
    for i in range(nn):
        ind_node = i + j*nn
        nodes[ ind_node, :] = [i*L/nc, j*L/nc, 0.0]

# connectivity matrix
connec = np.zeros((nelems,4),dtype=int)
# cross cables
for j in range(nc):
    for i in range(nc):
        ind_cell = i + j*nc
        connec[ind_cell*2  , : ] = [0, 0, j*nn+i  , (j+1)*nn+i+1 ]
        connec[ind_cell*2+1, : ] = [0, 0, j*nn+i+1, (j+1)*nn+i   ]

ind_cell = 2*nc*nc-1
for i in range(nc):
    ind_cell = ind_cell+1
    connec[ind_cell  , : ] = [0, 0, i  , i+1 ]
    ind_cell = ind_cell+1
    connec[ind_cell  , : ] = [0, 0, i+ nn*(nn-1)  , i+1+nn*(nn-1)]
    ind_cell = ind_cell+1
    connec[ind_cell  , : ] = [0, 0, i*nn  , (i+1)*nn ]
    ind_cell = ind_cell+1
    connec[ind_cell  , : ] = [0, 0, (i+1)*nn-1  , (i+2)*nn-1 ]

pos_z = L*.5

disp_mat = np.array([ [          0,    0, 0, pos_z   ],
                      [       nn-1,    L, 0,     0   ],
                      [  nn*(nn-1),    0, L,     0   ],
                      [    nn*nn-1,    L, L, pos_z   ]])

lref = L * sqrt(2)/(2*nc)
fext_mat  = np.zeros((nnodes,4))
for i in range(nnodes):
    fext_mat[i,:] = [ i, 0, 0.0, -gamma*areas[0]*lref*0.5 ] # node fx fy fz

model = cn.Model(nodes, connec, [linear_material], areas, disp_mat, fext_mat )
analy_sett = cn.AnalySettings()

# solve
nodes_defA, normal_forces, reactions, solu = cn.solve( model, analy_sett )

cn.plot( nodes, connec, nodes_defA, normal_forces, False, "A.png" )

# print(nodes_defA[:,2])

pos_z = L*.25

disp_mat = np.array([ [          0,    0, 0, pos_z   ],
                      [       nn-1,    L, 0,     0   ],
                      [  nn*(nn-1),    0, L,     0   ],
                      [    nn*nn-1,    L, L, pos_z   ]])

model = cn.Model(nodes, connec, [linear_material], areas, disp_mat, fext_mat )

# solve
nodes_defB1, normal_forces, reactions, solu = cn.solve( model, analy_sett )

cn.plot( nodes, connec, nodes_defB1, normal_forces, False, "B1.png" )

pos_z = L*.5

disp_mat = np.array([ [          0,    0, 0, pos_z   ],
                      [       nn-1,    L, 0,     0   ],
                      [  nn*(nn-1),    0, L,     0   ],
                      [    nn*nn-1,    L, L, pos_z   ]])
model = cn.Model(nodes, connec, [linear_material], areas, disp_mat, fext_mat )

# solve second part
nodes_defB2, normal_forces, reactions, solu = cn.solve( model, analy_sett, initsol=solu )

cn.plot( nodes, connec, nodes_defB2, normal_forces, False, "B2.png" )

dif_def = nodes_defA - nodes_defB2

# print("diference ", dif_def)
# print("total ", nodes_defA)

print(" norm diference ", np.linalg.norm( dif_def, ord=np.inf) )

print(" norm total ", np.linalg.norm( nodes_defA, ord=np.inf) )
