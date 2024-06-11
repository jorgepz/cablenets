# example 
import sys
sys.path.append('../src')

import numpy as np
import cablenets as cn

Lx = 2
Ly = 2
k = 2

nCellsX = 2
nCellsY = 2
nNodesX = nCellsX
nNodesY = nCellsY

nElems = nCellsX

nodes  = np.zeros( (nNodesX * nNodesY, 3) )
for i in range( nNodesX ):
  for j in range( nNodesY ):
    nodes[i+j*nNodesX,:] = [i*Lx/nCellsX, i*Ly/nCellsY, 0.0]

connec = np.zeros([ [0,1],
                    [1,2],
                    [2,0] ])

connec = np.array([ [0,1],
                    [1,2],
                    [2,0] ])

print("nodes", nodes)
print("conec", connec)

ks_vec = np.ones(np.shape(connec)[0])* k

disp_mat = np.array([ [0, 0,      0 , 0],
                      [1, Lx*1.1, 0  , 0] ])

fext_mat  = np.array([ [2, 1,1,0 ]]) # node fx fy fz

print("f", fext_mat)

nodes_def, normal_forces = cn.solve( nodes, connec, ks_vec, disp_mat, fext_mat )
print("normal forces", normal_forces)

cn.plot( nodes, connec, nodes_def, normal_forces )
