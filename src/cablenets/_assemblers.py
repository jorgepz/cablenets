
from cvxopt import matrix, spmatrix, spdiag, solvers, sparse
import numpy             as np


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

# ==================================================================================

# ==================================================================================

def _assemble_G_primal(nnodes, nelem, B ):
    vals = []
    Is   = []
    Js   = []
    print("assembling matrix G ...")
    for ele in range( nelem ):
        # -1*c
        Is.extend( [ (ele)*4 ] )
        Js.extend( [  ele    ] )
        vals.extend( [-1] )
    G_c = spmatrix( vals, Is, Js, (4*nelem, nelem))

    vals = []
    Is   = []
    Js   = []
    
    for ele in range( nelem ):
        Is.extend( range( ele*4+1,        ele*4+4 ) )
        Js.extend( range( ele*3+0  , ele*3+3      ) )
        vals.extend( (-1.0,-1.0,-1.0) )

    print(" ij", Is, Js)
    G_aux = spmatrix( vals, Is, Js, (4*nelem, 3*nelem))
    G_phi = G_aux * B
    print("gc size", G_c.size)
    print("gaux size", G_aux.size)
    print("gphi size", G_phi.size)
    print(" done.\n")
    return sparse( [[G_c], [G_phi] ] )

def _assemble_G_dual(n_dofs_d, nnodes, nelem ):
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

# ==================================================================================

# ==================================================================================

def _assemble_P_and_q_primal(nodes, connec, youngs, areas, nnodes, n_dofs_d, nelem, p_vec_zero_reactions ):
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
        vals.extend( [k] )
    print(" done.\n")
    P = spmatrix( vals, Is, Js, (nelem+3*nnodes, nelem+ 3*nnodes) )
    q = matrix( [ [matrix(0.0,(1,nelem))], [ matrix( -p_vec_zero_reactions ).trans() ] ] ).trans()

    return P, q

def _assemble_P_and_q_dual(nodes, connec, youngs, areas, n_dofs_d, nelem, d_vec ):
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
