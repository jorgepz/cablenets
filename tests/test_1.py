# test problem
import sys
sys.path.append('./src')

import numpy as np
import cablenets as cn

nodes  = np.array([  [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0] ])

connec = np.array([[0, 1 ], [1,2]])

ks_vec = np.ones((2))

disp_mat = np.array([ [ 0, 0   , 0, 0   ],
            [ 2, 2, 0, 0   ]])

fext_mat = np.array([[ 1, 1, 0.0, 0]]) # node fx fy fz

nodes_def, normal_forces = cn.solve( nodes, connec, ks_vec, disp_mat, fext_mat )

def test_normal_force():
    assert ( abs(normal_forces[1]) < 1e-7) and ( abs(normal_forces[0]-1)<1e-7 )
