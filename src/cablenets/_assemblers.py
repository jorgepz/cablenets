
from cvxopt import matrix, spmatrix, spdiag, solvers, sparse
import numpy             as np


def _assemble_B(nodes, connec, I ):
    connec = connec[:,2:4]
    nnodes = np.size( nodes, 0 )
    nelem  = np.size( connec, 0 )
    vals = []
    Is   = []
    Js   = []
    print("assembling matrix B ...",end='')
    for ele in range( nelem ):
        ini_node, end_node = connec[ele, :]
        Is.extend( range( ele*3     +0, (ele+1)*3      ) )
        Js.extend( range( ini_node*3+0, (ini_node+1)*3 ) )
        vals.extend( (-1,-1,-1) )
        
        Is.extend( range( ele*3     +0, (ele+1)*3      ) )
        Js.extend( range( end_node*3+0, (end_node+1)*3 ) )
        vals.extend( (1,1,1) )
    print("done.")
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

def _assemble_G_and_h_primal(nnodes, nelem, B, nodes, connec, isbil=False ):
    vals = []
    Is   = []
    Js   = []
    print("assembling matrix G and vector h ...",end='')

    if isbil:

        # matrix G
        for ele in range( nelem ):

            Is.extend( [ ele ] )
            Js.extend( [ ele ] )
            vals.extend( [-1] )

            Is.extend( [ ele+  nelem ] )
            Js.extend( [ ele+2*nelem    ] )
            vals.extend( [1] )

            Is.extend( [ (ele)*4+2*nelem ] )
            Js.extend( [  ele    ] )
            vals.extend( [-1] )

            Is.extend( [ (ele)*4+2*nelem ] )
            Js.extend( [  ele   +2*nelem ] )
            vals.extend( [-1] )

        G_c = spmatrix( vals, Is, Js, ((2+4)*nelem, 3*nelem))

        vals = []
        Is   = []
        Js   = []

        for ele in range( nelem ):
            Is.extend( range( 2*nelem+ele*4+1,  2*nelem + ele*4+4 ) )
            Js.extend( range( ele*3+0  , ele*3+3      ) )
            vals.extend( (-1.0,-1.0,-1.0) )

        G_aux = spmatrix( vals, Is, Js, ((2+4)*nelem, 3*nelem))
        G_phi = G_aux * B

        cvxG = sparse( [[G_c], [G_phi] ] )

        # vector h
        cvxh = matrix(0.0,(1,(2+1+3)*nelem)).trans()
        for ele in range( nelem ):
            ele_length = np.linalg.norm( nodes[ connec[ele,3],:] - nodes[ connec[ele,2],:])
            cvxh[nelem+ele]= .1 
            cvxh[2*nelem+ ele*4+0]= ele_length
    

    else:

        # matrix G
        for ele in range( nelem ):
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

        G_aux = spmatrix( vals, Is, Js, (4*nelem, 3*nelem))
        G_phi = G_aux * B

        # vector h
        cvxh = matrix(0.0,(1,(1+3)*nelem)).trans()
        for ele in range( nelem ):
            ele_length = np.linalg.norm( nodes[ connec[ele,3],:] - nodes[ connec[ele,2],:])
            cvxh[ele*4+0]= ele_length
    
        cvxG = sparse( [[G_c], [G_phi] ] )
    print(" done.")


    return cvxG, cvxh




def _assemble_G_dual(n_dofs_d, nnodes, nelem ):
    vals = []
    Is   = []
    Js   = []
    print("assembling matrix G ...",end='')
    for ele in range( nelem ):
        # q
        Is.extend( [ (ele)*4 ] )
        Js.extend( [  ele    ] )
        vals.extend( [-1] )
        
        # v
        Is.extend( range(       ele*4+1,        ele*4+4 ) )
        Js.extend( range( nelem+ele*3  , nelem+(ele+1)*3      ) )
        vals.extend( (-1.0,-1.0,-1.0) )
    print(" done.")
    return spmatrix( vals, Is, Js, ((1+3)*nelem, (1+3)*nelem + n_dofs_d))

# ==================================================================================

# ==================================================================================

def _assemble_P_and_q_primal(nodes, connec, materials, areas, nnodes, n_dofs_d, nelem, p_vec_zero_reactions, isbil ):
    
    youngs = np.zeros( (len(materials)))
    youngs_p = np.zeros( (len(materials)))
    sy = np.zeros( (len(materials)))

    if isbil:
        for ind in range(len(materials)):
            print(materials[ind])
            youngs[ind] = materials[ind].E
            youngs_p[ind] = materials[ind].Ep
            sy[ind] = materials[ind].sy
    else:
        for ind in range(len(materials)):
            print(materials[ind])
            youngs[ind] = materials[ind].E

    youngs = youngs[ connec[:,0]]
    if isbil:
        youngs_p = youngs_p[ connec[:,0]]
    sy = sy[ connec[:,0]]
    areas  = areas [ connec[:,1]]
    connec = connec[:,2:4]

    vals = []
    Is   = []
    Js   = []
    lengths = np.empty((nelem))

    print("assembling matrix P and q ...",end='')
    for ele in range( nelem ):

        if isbil:
            # q
            Is.extend(   [ 2*nelem+ele ] )
            Js.extend(   [ 2*nelem+ele ] )
            lengths[ele] = np.linalg.norm( nodes[ connec[ele,1],:] - nodes[ connec[ele,0],:])
            k = youngs[ele]*areas[ele]/lengths[ele]
            vals.extend( [k] )

            # q
            Is.extend(   [ ele ] )
            Js.extend(   [ ele ] )
            kp = youngs_p[ele]*areas[ele]/lengths[ele]
            vals.extend( [kp] )


        else:
            # q
            Is.extend(   [ ele ] )
            Js.extend(   [ ele ] )
            lengths[ele] = np.linalg.norm( nodes[ connec[ele,1],:] - nodes[ connec[ele,0],:])
            k = youngs[ele]*areas[ele]/lengths[ele]
            vals.extend( [k] )


    if isbil:
        print(" done.")
        P = spmatrix( vals, Is, Js, (3*nelem+3*nnodes, 3*nelem+ 3*nnodes) )
        q = matrix( [ [matrix(sy[0],(1,nelem))], [matrix(0.0,(1,2*nelem))], [ matrix( -p_vec_zero_reactions ).trans() ] ] ).trans()
    else:
        print(" done.")
        P = spmatrix( vals, Is, Js, (nelem+3*nnodes, nelem+ 3*nnodes) )
        q = matrix( [ [matrix(0.0,(1,nelem))], [ matrix( -p_vec_zero_reactions ).trans() ] ] ).trans()

    return P, q

def _assemble_P_and_q_dual(nodes, connec, materials, areas, n_dofs_d, nelem, d_vec, isbil ):

    youngs = np.zeros( (len(materials)))
    youngs_p = np.zeros( (len(materials)))
    sy = np.zeros( (len(materials)))

    if isbil:
        for ind in range(len(materials)):
            print(materials[ind])
            youngs[ind] = materials[ind].E
            youngs_p[ind] = materials[ind].Ep
            sy[ind] = materials[ind].sy
    else:
        for ind in range(len(materials)):
            print(materials[ind])
            youngs[ind] = materials[ind].E

    youngs = youngs[ connec[:,0]]
    if isbil:
        youngs_p = youngs_p[ connec[:,0]]
    sy = sy[ connec[:,0]]
    
    areas  = areas [ connec[:,1]]
    connec = connec[:,2:4]

    vals = []
    Is   = []
    Js   = []
    lengths = np.empty((nelem))

    print("assembling matrix P and q ...",end='')
    for ele in range( nelem ):
        # q
        Is.extend(   [ ele ] )
        Js.extend(   [ ele ] )
        lengths[ele] = np.linalg.norm( nodes[ connec[ele,1],:] - nodes[ connec[ele,0],:])
        k = youngs[ele]*areas[ele]/lengths[ele]
        vals.extend( [1/k] )
    print(" done.")
    P = spmatrix( vals, Is, Js, ((1+3)*nelem+n_dofs_d, (1+3)*nelem+n_dofs_d) )
    q = matrix( [ [matrix(lengths).trans()], [matrix(0.0,(1,3*nelem))], [ matrix( -np.array(d_vec) ).trans() ] ] ).trans()
    return P, q



def _assemble_A_and_b_primal(n_dofs_d, dofs_d, d_vec, nelem, nnodes):

    vals = []
    Is   = []
    Js   = []
    for ind in range(n_dofs_d):
        Is.extend( [ ind ]         )
        Js.extend( [ dofs_d[ind]+nelem ] )
        vals.extend( [1.0] )
    cvxA = spmatrix( vals, Is, Js, (n_dofs_d, nelem+3*nnodes) )
    cvxb = matrix( [ matrix( np.array(d_vec) ) ] )
    
    return cvxA, cvxb