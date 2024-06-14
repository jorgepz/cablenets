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
