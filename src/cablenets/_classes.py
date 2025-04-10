
class Material:
    def __init__(self, num_label, model_type, Eval, sy=[], epsy=[], Ep=[]):
        self.label = num_label
        if model_type=='linear':
            self.E = Eval

class Model:
    def __init__(self, nodes_coord_mat, elem_connec_mat, materials_list, areas_vec, def_coord_mat, fext_mat ):
        self.nodes = nodes_coord_mat
        self.connec = elem_connec_mat 
        self.materials = materials_list
        self.areas = areas_vec
        self.def_coord = def_coord_mat
        self.fext = fext_mat

class AnalySettings:
    def __init__(self, primal_dual_flag="primal", toler=[]):
        self.primal_dual = primal_dual_flag
    
    def set_pd_flag(self, primal_dual_flag):
        self.primal_dual = primal_dual_flag 