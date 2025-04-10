# test problem
# 2 element cable net
import sys
sys.path.append('./src')
import numpy as np
import cablenets as cn

# scalar parameters
L = 2.0 
E = 2.0 # young modulus
A = 1.0 # cross-section area

linear_material = cn.Material(0, 'linear', E )
areas = np.array([A])

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

# create model
model = cn.Model(nodes, connec, [linear_material], areas, disp_mat, fext_mat )

# create analysett
analy_sett = cn.AnalySettings()

# primal case
nodes_def_pri, normal_forces_pri, reactions = cn.solve( model, analy_sett )

# dual case
analy_sett.set_pd_flag("dual")
nodes_def_dua, normal_forces_dua, reactions = cn.solve( model, analy_sett )

print(" nodes def primal", nodes_def_pri)
print(" normal forc primal", normal_forces_pri)

print(" nodes def dual", nodes_def_dua)
print(" normal forc dual", normal_forces_dua)
cn.plot( nodes, connec, nodes_def_dua, normal_forces_dua, False )

def test_normal_force_primal():
    assert ( abs(normal_forces_pri[1]) < 1e-6) and ( abs(normal_forces_pri[0]-1)<1e-6 )

def test_normal_force_dual():
    assert ( abs(normal_forces_dua[1]) < 1e-6) and ( abs(normal_forces_dua[0]-1)<1e-6 )
