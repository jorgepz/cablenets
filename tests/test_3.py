# test problem
# 2 element cable net
import sys
sys.path.append('./src')
import numpy as np
from cablenets import solve, plot
import math

# scalar parameters
L = 1.0
youngs = np.array([2])
areas = np.array([1])

nodes = np.zeros((3,3))
nodes[0,:] = [0.0, 0.0, 0.0]
nodes[1,:] = [ -L,   L, 0.0]
nodes[2,:] = [ -L,  -L, 0.0]

connec = np.zeros((2,4),dtype=int)
connec[0,:] = [0, 0, 0, 1 ]
connec[1,:] = [0, 0, 0, 2 ]

disp_mat = np.array([ [     1, -L   , L, 0   ],
                      [     2, 0.0   , -L, 0.0   ]])

fext_mat      = np.zeros((1,4))
fext_mat[0,:] = [ 0, 0.0, -1.0, 0.0 ] # node fx fy fz

nodes_def, normal_forces = solve( nodes, connec, youngs, areas, disp_mat, fext_mat )

Aadd = np.array( [ [0,0,1,0,0,0,0,0,0,0,0], [0,0,0,0,1,0,0,0,0,0,0] ] )
badd = np.array( [0,0])

nodes_def, normal_forces = solve( nodes, connec, youngs, areas, disp_mat, fext_mat, "primal", A=Aadd, b=badd )

print(nodes_def)
print(normal_forces)

plot( nodes, connec, nodes_def, normal_forces, False )

k_verif = youngs[0]*areas[0]/(math.sqrt(2)*L)
numeric_cos = ( -nodes_def[0,1]+L) / ( normal_forces[0]/k_verif + (math.sqrt(2)*L) )

numeric_vertical_ext_force = normal_forces[0]*numeric_cos

def test_normal_force():
    assert ( abs(normal_forces[1]) < 1e-6) and ( abs(numeric_vertical_ext_force - fext_mat[0,2] )<1e-6 )
