'''
MIT License

Copyright (c) 2024 Jorge Pérez Zerpa

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import os

# load external dependencies
import numpy as np
from cvxopt import matrix, spmatrix, spdiag, solvers
import matplotlib.colors as colo
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# load cablenets functions
from cablenets.cablenets import _assemble_B, _assemble_d_or_p_vec, _assemble_P_and_q, _assemble_G, _remove_loads_in_fixed_nodes

#
# variables are x: [q,v,r] and s
#
def solve( nodes, connec, youngs, areas, disp_mat, fext_mat ):

    print( "\n=== Welcome to cablenets ===\n" )
    nnodes = np.size( nodes, 0 )
    nelem  = np.size( connec, 0 )
    print( "nodes", nnodes, "nelem: ", nelem,  )

    I = spdiag( matrix(1,(1,3)))
    B = _assemble_B(nodes, connec, I)    

    # assemble d and p
    d_vec, dofs_d = _assemble_d_or_p_vec(disp_mat)    
    p_vec, dofs_p = _assemble_d_or_p_vec(fext_mat)    

    p_vec, dofs_p = _remove_loads_in_fixed_nodes(p_vec, dofs_p, d_vec, dofs_d)

    dofs_d = dofs_d.tolist()
    dofs_p = dofs_p.tolist()
    
    BTp = B[:,dofs_p].trans()
    BTd = B[:,dofs_d].trans()

    n_dofs_d = len( dofs_d )
    n_dofs_p = len( dofs_p )

    cvxP, cvxq = _assemble_P_and_q(nodes, connec, youngs, areas, n_dofs_d, nelem, d_vec)    

    # primal-dual equality constraints
    cvxG = _assemble_G( n_dofs_d, nnodes, nelem )
    cvxh = matrix(0.0,(1,(1+3)*nelem)).trans()

    # primal equality constraints
    cvxA = matrix([
    [ spmatrix([],[],[], (n_dofs_p, nelem )), spmatrix([],[],[], (n_dofs_d, nelem )) ],
    [ BTp, BTd ],
    [ spmatrix( [],[],[], (n_dofs_p, n_dofs_d )), -spdiag( matrix(1.0, (1,n_dofs_d)) ) ]
    ]) 
    cvxb = matrix( [ matrix( np.array(p_vec) ) , matrix(0.0, (n_dofs_d, 1)) ] )

    # cone set
    cvxdims = {'l': 0, 'q': [4]*nelem , 's': []}

    solu = solvers.coneqp( cvxP, cvxq, cvxG, cvxh, cvxdims, cvxA, cvxb )   

    y = solu['y']
    x = solu['x']
    qs = x[0:nelem]
    vs = x[(nelem):(nelem+3*nelem)]
    nodes_def = np.zeros( np.size( y ))
    nodes_def[dofs_p] = -y[0:len(dofs_p)]
    nodes_def[dofs_d] = -y[len(dofs_p):(len(y))]
    nodes_def = np.reshape(nodes_def, (nnodes,3))
    normal_forces = (np.array(qs))
    normal_forces = np.reshape(normal_forces, (nelem,1))

    return nodes_def, normal_forces

# 
# 
#
def plot(nodes, connec, nodes_def, normal_forces, bool_show = True ):
    nnodes = np.size( nodes, 0 )
    nelem  = np.size( connec, 0 )
    connec = connec[:,2:4]

    print("nodes def", np.shape(nodes_def))
    print("type ", type(normal_forces))
    print("shape ", np.shape(normal_forces))

    max_normal_force = normal_forces.max()
    min_normal_force = 0# normal_forces.min()

    normali = colo.Normalize(vmin=min_normal_force, vmax=max_normal_force)

    m = cm.ScalarMappable(norm=normali, cmap='rainbow') # color

    fig = plt.figure()
    ax = fig.add_subplot( projection='3d')

    for ele in range( nelem ):
        ini_node, end_node = connec[ele, :]
        # if ele==0:
        #     legR='reference'
        #     legD='deformed'
        # else:
        #     legR = ''
        #     legD = ''

        ax.plot(    [nodes[ini_node,0], nodes[end_node,0]],
                    [nodes[ini_node,1], nodes[end_node,1]],
                 zs=[nodes[ini_node,2], nodes[end_node,2]], c='lightgray',linestyle='--')

        ax.plot(    [nodes_def[ini_node,0], nodes_def[end_node,0]],
                    [nodes_def[ini_node,1], nodes_def[end_node,1]],
                 zs=[nodes_def[ini_node,2], nodes_def[end_node,2]], c=m.to_rgba(normal_forces[ele]))
    # ax.legend()
    ax.axis('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    fig.colorbar(m, ax=ax)
    if bool_show:
        plt.show()
