# test problem
# 2 element cable net
import sys
sys.path.append('./src')
import numpy as np
from cablenets import solve, plot

# np.set_printoptions(threshold=sys.maxsize)

# scalar parameters
L = 2.0
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

nodes_def_pri, normal_forces_pri = solve( nodes, connec, youngs, areas, disp_mat, fext_mat, "primal" )

print(" nodes def", nodes_def_pri)
print(" normal forc ", normal_forces_pri)

nodes_def_dua, normal_forces_dua = solve( nodes, connec, youngs, areas, disp_mat, fext_mat, "dual" )

print(" nodes def", nodes_def_dua)
print(" normal forc ", normal_forces_dua)
plot( nodes, connec, nodes_def_dua, normal_forces_dua, False )

def test_normal_force_primal():
    assert ( abs(normal_forces_pri[1]) < 1e-6) and ( abs(normal_forces_pri[0]-1)<1e-6 )

def test_normal_force_dual():
    assert ( abs(normal_forces_dua[1]) < 1e-6) and ( abs(normal_forces_dua[0]-1)<1e-6 )
