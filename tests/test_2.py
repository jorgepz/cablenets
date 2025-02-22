# test problem
# 2 element cable net
import sys
sys.path.append('./src')
import numpy as np
from cablenets import solve, plot

# scalar parameters
L = 1
youngs = np.array([2])
areas = np.array([1])

nodes = np.zeros((4,3))
nodes[0,:] = [0.0, 0.0, 0.0]
nodes[1,:] = [ -L, L, 0.0]
nodes[2,:] = [0.0, L, 0.0]
nodes[3,:] = [  L, L, 0.0]

connec = np.zeros((3,4),dtype=int)
connec[0,:] = [0, 0, 0, 2 ]
connec[1,:] = [0, 0, 1, 2 ]
connec[2,:] = [0, 0, 3, 2 ]

disp_mat = np.array([ [     0,     0   , L*.5, 0   ],
                      [     1, -L*.5   , L*.75, 0   ],
                      [     3, L*1.25   , L*.75, 0   ]])

fext_mat      = np.zeros((1,4))
fext_mat[0,:] = [ 2, 0.0, 0.0, 0.0 ] # node fx fy fz

nodes_def, normal_forces = solve( nodes, connec, youngs, areas, disp_mat, fext_mat )

print(nodes_def)
print(normal_forces)

plot( nodes, connec, nodes_def, normal_forces, False )

def test_normal_force():
    assert ( abs(normal_forces[1]) < 1e-6) and ( abs(normal_forces[0])<1e-6 )
