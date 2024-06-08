
import numpy as np
import socp_cable_solver as scs

Lx = 2
Ly = 2
k = 2
nelems = 6
nnodes = 4

nodes  = np.array([ [  0,  0, 0 ],
                    [ Lx,  0, 0 ],
                    [  0, Ly, 0 ] ])

connec = np.array([ [0,1],
                    [1,2],
                    [2,0] ] ,dtype=int)

print("nodes", nodes)
print("conec", connec)

ks_vec = np.ones((nelems))* k

disp_mat = np.array([ [0, 0,      0 , 0],
                      [1, Lx*1.1, 0  , 0] ])

fext_mat  = np.array([ [2, 1,1,0 ]]) # node fx fy fz

print("f", fext_mat)

nodes_def, normal_forces = scs.solve_socp( nodes, connec, ks_vec, disp_mat, fext_mat )
print("normal forces", normal_forces)
scs.plot( nodes, connec, nodes_def, normal_forces )
