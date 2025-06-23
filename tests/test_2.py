# test problem
# 2 element cable net
import sys
sys.path.append('./src')
import numpy as np
import cablenets as cn

# scalar parameters
L = 1
E = 2.0 # young modulus
A = 1.0 # cross-section area

linear_material = cn.Material(0, 'linear', E )
areas = np.array([A])

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

model = cn.Model(nodes, connec, [linear_material], areas, disp_mat, fext_mat )
analy_sett = cn.AnalySettings()

nodes_def, normal_forces, reactions, solu = cn.solve( model, analy_sett )

print(nodes_def)
print(normal_forces)

cn.plot( nodes, connec, nodes_def, normal_forces, False )

def test_normal_force():
    assert ( abs(normal_forces[1]) < 1e-6) and ( abs(normal_forces[0])<1e-6 )
