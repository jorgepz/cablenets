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

# import cvx
from cvxopt import matrix, spmatrix, spdiag, solvers

# import dependencies
import numpy             as np
import matplotlib.colors as colo
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# input nparray matrices
def _assemble_B(nodes, connec, I ):
    connec = connec[:,2:4]
    nnodes = np.size( nodes, 0 )
    nelem  = np.size( connec, 0 )
    vals = []
    Is   = []
    Js   = []
    print("assembling matrix B ...")
    for ele in range( nelem ):
        print(ele)
        ini_node, end_node = connec[ele, :]
        Is.extend( range( ele*3     +0, (ele+1)*3      ) )
        Js.extend( range( ini_node*3+0, (ini_node+1)*3 ) )
        vals.extend( (-1,-1,-1) )
        
        Is.extend( range( ele*3     +0, (ele+1)*3      ) )
        Js.extend( range( end_node*3+0, (end_node+1)*3 ) )
        vals.extend( (1,1,1) )
    print("done.\n")
    return spmatrix( vals, Is, Js, (3*nelem, 3*nnodes))

#
# input nparray matrix
#
def _assemble_d_or_p_vec( mat ):
    n_nodes = np.size( mat, 0 )
    vec  = []
    dofs = []
    # assemble
    for i in range( n_nodes ):
        this_node_dofs = list( range( int(mat[i,0]*3), int((mat[i,0]+1)*3) ) )
        dofs.extend( this_node_dofs )
        vec.extend( mat[i,1:] )
    vec  = np.array( vec  )
    dofs = np.array( dofs, dtype=int )
    return vec, dofs

def _assemble_G(n_dofs_d, nnodes, nelem ):
    vals = []
    Is   = []
    Js   = []
    print("assembling matrix G ...")
    for ele in range( nelem ):
        # q
        Is.extend( [ (ele)*4 ] )
        Js.extend( [  ele    ] )
        vals.extend( [-1] )
        
        # v
        Is.extend( range(       ele*4+1,        ele*4+4 ) )
        Js.extend( range( nelem+ele*3  , nelem+(ele+1)*3      ) )
        vals.extend( (-1.0,-1.0,-1.0) )
    print(" done.\n")
    return spmatrix( vals, Is, Js, ((1+3)*nelem, (1+3)*nelem + n_dofs_d))

def _assemble_P_and_q(nodes, connec, youngs, areas, n_dofs_d, nelem, d_vec ):
    youngs = youngs[ connec[:,0]]
    areas  = areas [ connec[:,1]]
    connec = connec[:,2:4]

    vals = []
    Is   = []
    Js   = []
    lengths = np.empty((nelem))

    print("assembling matrix P and q ...")
    for ele in range( nelem ):
        # q
        Is.extend(   [ ele ] )
        Js.extend(   [ ele ] )
        lengths[ele] = np.linalg.norm( nodes[ connec[ele,1],:] - nodes[ connec[ele,0],:])
        k = youngs[ele]*areas[ele]/lengths[ele]
        vals.extend( [1/k] )
    print(" done.\n")
    P = spmatrix( vals, Is, Js, ((1+3)*nelem+n_dofs_d, (1+3)*nelem+n_dofs_d) )
    q = matrix( [ [matrix(lengths).trans()], [matrix(0.0,(1,3*nelem))], [ matrix( -np.array(d_vec) ).trans() ] ] ).trans()
    print("P", P)
    return P, q


def _remove_loads_in_fixed_nodes(p_vec, dofs_p, d_vec, dofs_d):
    filtered_dofs_p = dofs_p
    filtered_p_vec  = p_vec
    for indd in dofs_d:
        inde = np.where( filtered_dofs_p == indd )
        if len(inde[0])>0:
            filtered_p_vec  = np.delete( filtered_p_vec,  inde[0])
            filtered_dofs_p = np.delete( filtered_dofs_p, inde[0])

    return filtered_p_vec, filtered_dofs_p 


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
    p_vec_aux, dofs_p_aux = _assemble_d_or_p_vec(fext_mat)    

    p_vec  = np.zeros( 3*nnodes )
    dofs_p = np.arange( 3*nnodes )
    print("dofs d", dofs_d)
    print("type dofs p aux", type(dofs_p_aux))
    print("dofs p aux", dofs_p_aux.shape )
    print("vec p aux", p_vec_aux.shape )
    print("vec p", p_vec.shape )
    p_vec[dofs_p_aux] = p_vec[dofs_p_aux] + p_vec_aux

    print("dofs p", dofs_p)

    p_vec  = np.delete(p_vec, dofs_d)
    dofs_p = np.delete(dofs_p,dofs_d) 
    
    print("dofs p post remove", dofs_p)
    print("pvec", p_vec)

    p_vec, dofs_p = _remove_loads_in_fixed_nodes(p_vec, dofs_p, d_vec, dofs_d)

    dofs_d = dofs_d.tolist()
    dofs_p = dofs_p.tolist()
    
    print("dofs d", dofs_d)
    print("dofs p", dofs_p)

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
    nodes_def = np.zeros( (nnodes*3,1 ))
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

    # colormap = 'YlGnBu'
    colormap = 'rainbow'
    # colormap = 'Set1'
    
    nnodes = np.size( nodes, 0 )
    nelem  = np.size( connec, 0 )
    connec = connec[:,2:4]

    print("nodes def", np.shape(nodes_def))
    print("type ", type(normal_forces))
    print("shape ", np.shape(normal_forces))

    max_normal_force = normal_forces.max()
    min_normal_force = 0# normal_forces.min()

    normali = colo.Normalize(vmin=min_normal_force, vmax=max_normal_force)

    m = cm.ScalarMappable(norm=normali, cmap=colormap) # color

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
    ax.set_title('normal forces')
    fig.colorbar(m, ax=ax)
    if bool_show:
        plt.show()
