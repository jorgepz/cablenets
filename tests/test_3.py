# test problem
# 2 element cable net
import sys
sys.path.append('./src')
import numpy as np
import cablenets as cn
import math

# scalar parameters
L = 1.0
E = 2.0 # young modulus
A = 1.0 # cross-section area

linear_material = cn.Material(0, 'linear', E )
areas = np.array([A])

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

model = cn.Model(nodes, connec, [linear_material], areas, disp_mat, fext_mat )
analy_sett = cn.AnalySettings()

nodes_def, normal_forces, reactions, solu = cn.solve( model, analy_sett )

Aadd = np.array( [ [0,0,1,0,0,0,0,0,0,0,0], [0,0,0,0,1,0,0,0,0,0,0] ] )
badd = np.array( [0,0])

nodes_def, normal_forces, reactions, solu = cn.solve( model, analy_sett, A=Aadd, b=badd )

print(nodes_def)
print(normal_forces)

cn.plot( nodes, connec, nodes_def, normal_forces, False )

k_verif = E*A/(math.sqrt(2)*L)
numeric_cos = ( -nodes_def[0,1]+L) / ( normal_forces[0]/k_verif + (math.sqrt(2)*L) )

numeric_vertical_ext_force = normal_forces[0]*numeric_cos

def test_normal_force():
    assert ( abs(normal_forces[1]) < 1e-6) and ( abs(numeric_vertical_ext_force - abs(fext_mat[0,2]) )<1e-4 )
