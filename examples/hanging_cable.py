# hanging cable
import sys
from os.path import dirname
sys.path.append(dirname('../src/'))
import numpy as np
from math import pi
import cablenets as cn

# input scalar parameters
L       = 2     # length
E       = 200e9 # young modulus (Pa)
d       = .01   # cross-section diameter (m)
nelems  = 20    # number of mesh elements
density = 7850  # density (kg/m3)

# create materials and areas vectors
linear_material = cn.Material(0, 'linear', E )
areas  = np.array([ pi*d**2/4 ])

# create nodes matrix
nodes  = np.zeros((nelems+1,3))
for i in range(nelems+1):
    nodes[i,:] = [i*L/nelems, 0.0, 0.0]

# define connectivity
connec = np.zeros((nelems,4),dtype=int)
for i in range(nelems):
    connec[i,:] = [0, 0, i, i+1 ]

# deformed node positions
disp_mat = np.array([ [          0, 0   , 0, 0   ],
                      [     nelems, 0.8*L   , 0, 0   ]])

# external forces
fext_mat  = np.zeros((nelems-1,4))
for i in range(nelems-1):
    fext_mat[i,:] = [ i+1, 0, 0.0, -density*areas[0]*L*9.81/nelems ] # node fx fy fz

model = cn.Model(nodes, connec, [linear_material], areas, disp_mat, fext_mat )
analy_sett = cn.AnalySettings()

# solve
nodes_def, normal_forces, reactions = cn.solve( model, analy_sett )

# plot
cn.plot( nodes, connec, nodes_def, normal_forces )

