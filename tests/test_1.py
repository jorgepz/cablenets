# test problem
# 2 element cable net
import sys
sys.path.append('./src')
import numpy as np
from cablenets import solve, plot

# scalar parameters
L = 2
youngs = np.array([2])
areas = np.array([1])
nelems = 2

nodes = np.zeros((nelems+1,3))
for i in range(nelems+1):
    nodes[i,:] = [i*L/nelems, 0.0, 0.0]

connec = np.zeros((nelems,4),dtype=int)
for i in range(nelems):
    connec[i,:] = [0, 0, i, i+1 ]

disp_mat = np.array([ [          0, 0   , 0, 0   ],
                      [     nelems, L   , 0, 0   ]])

fext_mat      = np.zeros((1,4))
fext_mat[0,:] = [ 1, 1.0, 0.0, 0.0 ] # node fx fy fz

nodes_def, normal_forces = solve( nodes, connec, youngs, areas, disp_mat, fext_mat )

plot( nodes, connec, nodes_def, normal_forces, False )

def test_normal_force():
    assert ( abs(normal_forces[1]) < 1e-6) and ( abs(normal_forces[0]-1)<1e-6 )
