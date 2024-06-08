
from cvxopt import matrix, spmatrix, spdiag, solvers

import numpy             as np
import matplotlib.colors as colo
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# input nparray matrices
def _assemble_B(nodes, connec, I ):
    nnodes = np.size( nodes, 0 )
    nelem  = np.size( connec, 0 )
    vals = []
    Is   = []
    Js   = []
    print("assembling matrix B \n")
    for ele in range( nelem ):
        ini_node, end_node = connec[ele, :]
        print( "elem :", ele, "  ini node", ini_node, " end node", end_node )
        Is.extend( range( ele*3     +0, (ele+1)*3      ) )
        Js.extend( range( ini_node*3+0, (ini_node+1)*3 ) )
        vals.extend( (-1,-1,-1) )
        
        Is.extend( range( ele*3     +0, (ele+1)*3      ) )
        Js.extend( range( end_node*3+0, (end_node+1)*3 ) )
        vals.extend( (1,1,1) )
    return spmatrix( vals, Is, Js, (3*nelem, 3*nnodes))

#
# input nparray matrix
#
def _assemble_d_or_p_vec( mat ):
    n_nodes = np.size( mat, 0 )
    print(" n_nodes ", n_nodes )
    
    vec  = []
    dofs = []

    # assemble
    for i in range( n_nodes ):
        this_node_dofs = list( range( int(mat[i,0]*3), int((mat[i,0]+1)*3) ) )
        print("CACAAA", type(this_node_dofs)) 
        dofs.extend( this_node_dofs )
        vec.extend( mat[i,1:] )
        print(" this dofs ", dofs )
        # print(" vec ", vec[ this_node_dofs,:] )
        # print(" vec ", vec )
        print(" mat ", mat[ i,1:] )
        # vec[ (i*3+0):(i+1)*3,: ] = mat[i,1:-1]
    return vec, dofs

def _assemble_G(n_dofs_d, nnodes, nelem ):
    vals = []
    Is   = []
    Js   = []
    print("assembling matrix G \n")
    for ele in range( nelem ):
        # q
        Is.extend( [ (ele)*4 ] )
        Js.extend( [  ele    ] )
        vals.extend( [-1] )
        
        # v
        Is.extend( range(       ele*4+1,        ele*4+4 ) )
        Js.extend( range( nelem+ele*3  , nelem+(ele+1)*3      ) )
        vals.extend( (-1.0,-1.0,-1.0) )
    return spmatrix( vals, Is, Js, ((1+3)*nelem, (1+3)*nelem + n_dofs_d))

def _assemble_P_and_q(nodes, connec, ks_vec, n_dofs_d, nelem, d_vec ):
    vals = []
    Is   = []
    Js   = []
    lengths = np.empty((nelem))
    print("assembling matrix P and q \n")
    for ele in range( nelem ):
        # q
        Is.extend(   [ ele ] )
        Js.extend(   [ ele ] )
        vals.extend( [1/ks_vec[ele]] )
        lengths[ele] = np.linalg.norm( nodes[ connec[ele,1],:] - nodes[ connec[ele,0],:])
    print("largos ", lengths)
    print("vals", vals)
    
    print("sizes ", ( matrix(lengths).trans()).size)
    print((matrix(0,(1,3*nelem))).size )
    print(( matrix( np.array(d_vec) ).trans()  ).size)
    P = spmatrix( vals, Is, Js, ((1+3)*nelem+n_dofs_d, (1+3)*nelem+n_dofs_d) )
    q = matrix( [ [matrix(lengths).trans()], [matrix(0.0,(1,3*nelem))], [ matrix( -np.array(d_vec) ).trans() ] ] ).trans()
    return P, q



#
# variables are x: [q,v,r] and s
#
def solve( nodes, connec, ks_vec, disp_mat, fext_mat ):

    print( "\n=== Welcome to cablenets ===\n" )
    nnodes = np.size( nodes, 0 )
    nelem  = np.size( connec, 0 )

    print( "nodes", nnodes, "nelem: ", nelem,  )
    # dims |@··@|~= {'l': 0, 'q': [n+1], 's': []}

    I = spdiag( matrix(1,(1,3)))
    print("I :\n", I)
    B = _assemble_B(nodes, connec, I)    
    print("B :\n", B)

    # assemble d and p
    d_vec, dofs_d = _assemble_d_or_p_vec(disp_mat)    
    p_vec, dofs_p = _assemble_d_or_p_vec(fext_mat)    

    print(" dvec ", d_vec)
    print(" pvec ", p_vec)
    print(" dofs d ", dofs_d)
    print(" dofs p ", dofs_p)
    
    BTp = B[:,dofs_p].trans()
    BTd = B[:,dofs_d].trans()

    n_dofs_d = len( dofs_d )
    n_dofs_p = len( dofs_p )

    cvxP, cvxq = _assemble_P_and_q(nodes, connec, ks_vec, n_dofs_d, nelem, d_vec)    

    print("B size", B.size, "   BT \n", BTd, BTd.size, "ndofs p", n_dofs_p, "ndofs d", n_dofs_d)

    # primal-dual equality constraints
    cvxG = _assemble_G( n_dofs_d, nnodes, nelem )
    cvxh = matrix(0.0,(1,(1+3)*nelem)).trans()
    print("G ", cvxG )
    print("size G",cvxG.size)
    print("size h",cvxh.size)

    # primal equality constraints
    cvxA = matrix([
    [ spmatrix([],[],[], (n_dofs_p, nelem )), spmatrix([],[],[], (n_dofs_d, nelem )) ],
    [ BTp, BTd ],
    [ spmatrix( [],[],[], (n_dofs_p, n_dofs_d )), -spdiag( matrix(1.0, (1,n_dofs_d)) ) ]
    ]) 
    cvxb = matrix( [ matrix( np.array(p_vec) ) , matrix(0.0, (n_dofs_d, 1)) ] )
    print("size A ",cvxA.size)
    print("size b ",cvxb.size)

    # import sys
    #np.set_printoptions(threshold=np.inf,linewidth=10000)
    print(" ================ ")
    print("P",cvxP)
    print("q",cvxq)
    print("size A ",cvxA.size)
    print(" A ",cvxA[:,0:2])
    print(" A ",cvxA[:,2:8])
    print(" A ",cvxA[:,8:14])
    print(" b ",cvxb)

    print(" G ",cvxG[:,0:2])
    print(" G ",cvxG[:,2:8])
    print(" G ",cvxG[:,8:14])
    print(" h ",cvxh)
    # cone set

    cvxdims = {'l': 0, 'q': [4]*nelem , 's': []}
    print("dims", cvxdims['q'])

    solu = solvers.coneqp( cvxP, cvxq, cvxG, cvxh, cvxdims, cvxA, cvxb )   

    y = solu['y']
    print("solu z", solu['z'])
    x = solu['x']
    qs = x[0:nelem]
    vs = x[(nelem):(nelem+3*nelem)]
    print("xs: ", x)
    print("y: ", y)
    print("qs: ", qs)
    print("vs: ", vs)
    nodes_def = np.zeros( np.size( y ))
    print("U: ", nodes_def)
    print("U: ", len(dofs_p))
    nodes_def[dofs_p] = -y[0:len(dofs_p)]
    nodes_def[dofs_d] = -y[len(dofs_p):(len(y))]
    print("U: ", nodes_def)
    nodes_def = np.reshape(nodes_def, (nnodes,3))
    print("U: ", nodes_def)
    print("size: ", np.shape(nodes_def) )
    normal_forces = (np.array(qs))
    normal_forces = np.reshape(normal_forces, (nelem,1))
    return nodes_def, normal_forces

# 
# 
#
def plot(nodes, connec, nodes_def, normal_forces ):
    nnodes = np.size( nodes, 0 )
    nelem  = np.size( connec, 0 )

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
        if ele==0:
            legR='reference'
            legD='deformed'
        else:
            legR = ''
            legD = ''

        ax.plot(    [nodes[ini_node,0], nodes[end_node,0]],
                    [nodes[ini_node,1], nodes[end_node,1]],
                 zs=[nodes[ini_node,2], nodes[end_node,2]], label=legR, c='gray')

        ax.plot(    [nodes_def[ini_node,0], nodes_def[end_node,0]],
                    [nodes_def[ini_node,1], nodes_def[end_node,1]],
                 zs=[nodes_def[ini_node,2], nodes_def[end_node,2]], label=legD, c=m.to_rgba(normal_forces[ele]))
    ax.legend()
    ax.axis('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    fig.colorbar(m, ax=ax)
    plt.show()

