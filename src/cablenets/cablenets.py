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

# import  cvx
from cvxopt import matrix, spmatrix, spdiag, solvers, sparse

# import dependencies
import numpy             as np
import matplotlib.colors as colo
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from cablenets._assemblers import _assemble_B, _assemble_d_or_p_vec
from cablenets._assemblers import _assemble_P_and_q_primal, _assemble_P_and_q_dual
from cablenets._assemblers import _assemble_G_dual, _assemble_G_and_h_primal
from cablenets._assemblers import _assemble_A_and_b_primal

from cablenets._classes import Material, Model, AnalySettings

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
def solve( model, analy_sett, **kwargs ):

    nodes = model.nodes
    connec = model.connec
    materials = model.materials
    areas = model.areas
    youngs = np.zeros( (len(materials)))
    youngs_p = np.zeros( (len(materials)))
    cys = np.zeros( (len(materials)))
    qys = np.zeros( (len(materials)))
    for ind in range(len(materials)):
        youngs[ind] = materials[ind].E
        if materials[ind].model_type == 'bilinear':
            youngs_p[ind] = materials[ind].Ep
            cys[ind] = materials[ind].epsy
            qys[ind] = materials[ind].sy
            isbil = True
        else:
            isbil = False
            
    def_coord_mat = model.def_coord
    fext_mat = model.fext

    primal_dual_flag = analy_sett.primal_dual

    if "A" in kwargs:
        print(kwargs["A"])
    print( "\n=== Welcome to cablenets ===\n" )
    nnodes = np.size( nodes, 0 )
    nelem  = np.size( connec, 0 )
    print( "nodes", nnodes, "nelem: ", nelem,  )

    I = spdiag( matrix(1,(1,3)))
    B = _assemble_B(nodes, connec, I)    

    # --------------------------------------------   
    # assemble d and p  and compute dofs
    d_vec, dofs_d = _assemble_d_or_p_vec(def_coord_mat)   
    x_vec = np.reshape(nodes,(nnodes*3))[dofs_d] # reference coodinates
    p_vec_aux, dofs_p_aux = _assemble_d_or_p_vec(fext_mat)    

    p_vec  = np.zeros( 3*nnodes )
    dofs_p = np.arange( 3*nnodes )
    p_vec[dofs_p_aux] = p_vec[dofs_p_aux] + p_vec_aux

    p_vec_zero_reactions = p_vec

    p_vec  = np.delete(p_vec, dofs_d)
    dofs_p = np.delete(dofs_p, dofs_d) 
    
    p_vec, dofs_p = _remove_loads_in_fixed_nodes(p_vec, dofs_p, d_vec, dofs_d)
    
    dofs_d = dofs_d.tolist()
    dofs_p = dofs_p.tolist()

    # --------------------------------------------   

    if primal_dual_flag == "primal":

        n_dofs_d = len( dofs_d )
        n_dofs_p = len( dofs_p )

        cvxP, cvxq = _assemble_P_and_q_primal(nodes, connec, materials, areas, nnodes, n_dofs_d, nelem, p_vec_zero_reactions, isbil)    

        # primary-slack equality constraints
        cvxG, cvxh = _assemble_G_and_h_primal( nnodes, nelem, B, nodes, connec, isbil )

        cvxA, cvxb = _assemble_A_and_b_primal(n_dofs_d, dofs_d, d_vec, nelem, nnodes)

        if "A" in kwargs:
            cvxA = sparse( [ [cvxA, matrix( kwargs["A"])] ])
            cvxb = matrix( [ cvxb, matrix( kwargs["b"]) ] )

    elif primal_dual_flag=="dual":

        BTp = B[:,dofs_p].trans()
        BTd = B[:,dofs_d].trans()

        n_dofs_d = len( dofs_d )
        n_dofs_p = len( dofs_p )

        cvxP, cvxq = _assemble_P_and_q_dual(nodes, connec, materials, areas, n_dofs_d, nelem, d_vec, isbil)    

        # primary-slack equality constraints
        cvxG = _assemble_G_dual( n_dofs_d, nnodes, nelem )
        cvxh = matrix(0.0,(1,(1+3)*nelem)).trans()

        # primary vars equality constraints
        cvxA = matrix([
        [ spmatrix([],[],[], (n_dofs_p, nelem )), spmatrix([],[],[], (n_dofs_d, nelem )) ],
        [ BTp, BTd ],
        [ spmatrix( [],[],[], (n_dofs_p, n_dofs_d )), -spdiag( matrix(1.0, (1,n_dofs_d)) ) ]
        ]) 
        cvxb = matrix( [ matrix( np.array(p_vec) ) , matrix(0.0, (n_dofs_d, 1)) ] )

    # cone set
    cvxdims = {'l': 0, 'q': [4]*nelem , 's': []}

    print("call to CVXOPT:")

    solvers.options['abstol']=1e-12
    solvers.options['reltol']=1e-12


    if "initsol" in kwargs:
        init_x = kwargs["initsol"]["x"]

        # apply the current dirichlet conditions
        for ind,val in enumerate(dofs_d):
            # print("ind", ind, " val ", val)
            init_x[nelem+val] = d_vec[ind] 

        # compute new elongations
        vecaux = B*init_x[nelem:(nelem+3*nnodes)]
        for ele in range( nelem ):
            ele_length = np.linalg.norm( nodes[ connec[ele,3],:] - nodes[ connec[ele,2],:])
            init_x[ele] = np.linalg.norm( vecaux[ele*3:(ele+1)*3] ) - ele_length + 1.0e0

        # update primal slack variables 
        init_s = cvxh - cvxG*init_x 
        init_dict = {'x': init_x, 's': init_s, 'y':kwargs["initsol"]['y'], 'z':kwargs["initsol"]['z'] }
        
        solu = solvers.coneqp( cvxP, cvxq, cvxG, cvxh, cvxdims, cvxA, cvxb, initvals=init_dict )   

    else:
        solu = solvers.coneqp( cvxP, cvxq, cvxG, cvxh, cvxdims, cvxA, cvxb )   

    x = solu['x']
    y = solu['y']
    s = solu['s']
    z = solu['z']

    #
    if primal_dual_flag == "primal":
        cs = x[0:nelem]
        nodes_def = np.reshape( x[(nelem):(nelem+3*nnodes)], (nnodes,3))
        normal_forces = 0.0
        reactions = -y
        normal_forces = np.zeros( (nelem, 1))
        for j in range(nelem):
            normal_forces[j] = z[j*4]
    #
    elif primal_dual_flag == "dual":
        qs = x[0:nelem]
        vs = x[(nelem):(nelem+3*nelem)]
        nodes_def = np.zeros( (nnodes*3,1 ))
        nodes_def[dofs_p] = -y[0:len(dofs_p)]
        nodes_def[dofs_d] = -y[len(dofs_p):(len(y))]
        nodes_def = np.reshape(nodes_def, (nnodes,3))
        normal_forces = (np.array(qs))
        normal_forces = np.reshape(normal_forces, (nelem,1))
        reactions = 0.0
    
    return nodes_def, normal_forces, reactions, solu

# 
# 
#
def plot(nodes, connec, nodes_def, normal_forces, bool_show = True, fig_filename = "" ):

    # colormap = 'YlGnBu'
    colormap = 'rainbow'
    # colormap = 'Set1'
    
    nnodes = np.size( nodes, 0 )
    nelem  = np.size( connec, 0 )
    connec = connec[:,2:4]

    print("nodes def", np.shape(nodes_def))
    print("type ", type(normal_forces))
    print("shape ", np.shape(normal_forces))

    max_normal_force = max( normal_forces.max(), 0)
    min_normal_force = 0 # normal_forces.min()

    normali = colo.Normalize(vmin=min_normal_force, vmax=max_normal_force)

    m = cm.ScalarMappable(norm=normali, cmap=colormap) # color

    # fig = plt.figure()
    # ax = fig.add_subplot( projection='3d')
    # for ele in range( nelem ):
    #     ini_node, end_node = connec[ele, :]
    #     ax.plot(    [nodes[ini_node,0], nodes[end_node,0]],
    #                 [nodes[ini_node,1], nodes[end_node,1]],
    #              zs=[nodes[ini_node,2], nodes[end_node,2]], c='lightgray',linestyle='--')

    #     ax.plot(    [nodes_def[ini_node,0], nodes_def[end_node,0]],
    #                 [nodes_def[ini_node,1], nodes_def[end_node,1]],
    #              zs=[nodes_def[ini_node,2], nodes_def[end_node,2]], c=m.to_rgba(normal_forces[ele]))
    # ax.axis('equal')
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    # ax.set_zlabel('z')
    # ax.set_title('normal forces')
    # fig.colorbar(m, ax=ax)

    # fig = plt.figure()
    # ax = fig.add_subplot()

    fig = plt.figure()
    ax = fig.add_subplot( projection='3d')

    for ele in range( nelem ):
        ini_node, end_node = connec[ele, :]
        ax.plot(    [nodes[ini_node,0], nodes[end_node,0]],
                    [nodes[ini_node,1], nodes[end_node,1]],
                 zs=[nodes[ini_node,2], nodes[end_node,2]], c='lightgray',linestyle='--')
        
        ax.plot(    [nodes_def[ini_node,0], nodes_def[end_node,0]],
                    [nodes_def[ini_node,1], nodes_def[end_node,1]],
                 zs=[nodes_def[ini_node,2], nodes_def[end_node,2]], c=m.to_rgba(normal_forces[ele]))

        # ax.plot(    [nodes_def[ini_node,0], nodes_def[end_node,0]],
        #             [nodes_def[ini_node,1], nodes_def[end_node,1]], c=m.to_rgba(normal_forces[ele]))
    ax.axis('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('normal forces')
    fig.colorbar(m, ax=ax)

    if len(fig_filename):
        plt.savefig( fig_filename)

    if bool_show:
        plt.show()
