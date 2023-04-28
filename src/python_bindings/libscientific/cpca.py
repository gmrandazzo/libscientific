"""
cpca libscientific python binding

Copyright (C) <2022>  Giuseppe Marco Randazzo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import ctypes
import libscientific.vectlist as vlst
import libscientific.matrix as mx
import libscientific.tensor as tns
import libscientific.vector as vect
from libscientific.loadlibrary import load_libscientific_library

lsci = load_libscientific_library()


class CPCAMODEL(ctypes.Structure):
    """
    CPCA model data structure
    """
    _fields_ = [
        ("block_scores", ctypes.POINTER(tns.TENSOR)),
        ("block_loadings", ctypes.POINTER(tns.TENSOR)),
        ("super_scores", ctypes.POINTER(mx.MATRIX)),
        ("super_weights", ctypes.POINTER(mx.MATRIX)),
        ("scaling_factor", ctypes.POINTER(vect.DVECTOR)),
        ("total_expvar", ctypes.POINTER(vect.DVECTOR)),
        ("block_expvar", ctypes.POINTER(vlst.DVECTLIST)),
        ("colaverage", ctypes.POINTER(vlst.DVECTLIST)),
        ("colscaling", ctypes.POINTER(vlst.DVECTLIST))]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__



lsci.NewCPCAModel.argtypes = [ctypes.POINTER(ctypes.POINTER(CPCAMODEL))]
lsci.NewCPCAModel.restype = None


def new_cpca_model():
    """
    new_cpca_model: Allocate in memory an empty libscientific CPCA model
    """
    mpca = ctypes.POINTER(CPCAMODEL)()
    lsci.NewCPCAModel(ctypes.pointer(mpca))
    return mpca


lsci.DelCPCAModel.argtypes = [ctypes.POINTER(ctypes.POINTER(CPCAMODEL))]
lsci.DelCPCAModel.restype = None


def del_cpca_model(mpca):
    """
    del_cpca_model: Delete an allocated libscientific CPCA model
    """
    lsci.DelCPCAModel(ctypes.pointer(mpca))


lsci.CPCA.argtypes = [ctypes.POINTER(tns.TENSOR),
                     ctypes.c_size_t,
                     ctypes.c_size_t,
                     ctypes.POINTER(CPCAMODEL),
                     ctypes.POINTER(ctypes.c_int)]
lsci.CPCA.restype = None


def cpca_algorithm(t_input, scaling, npc, mpca):
    """
    CPCA: Calculate the CPCA model for a matrix m using the NIPALS algorithm
    """
    ssignal = ctypes.c_int(0)
    lsci.CPCA(t_input, scaling, npc, mpca, ssignal)


lsci.CPCAScorePredictor.argtypes = [ctypes.POINTER(tns.TENSOR),
                                   ctypes.POINTER(CPCAMODEL),
                                   ctypes.c_size_t,
                                   ctypes.POINTER(mx.MATRIX),
                                   ctypes.POINTER(tns.TENSOR)]
lsci.CPCAScorePredictor.restype = None


def cpca_score_predictor(t_input, model, npc, p_super_scores, p_block_scores):
    """
    cpca_score_predictor: Predict scores for a matrix m in the computed CPCA modes
    """
    lsci.CPCAScorePredictor(t_input,
                            model,
                            npc,
                            p_super_scores,
                            p_block_scores)


lsci.PrintCPCA.argtypes = [ctypes.POINTER(CPCAMODEL)]
lsci.PrintCPCA.restype = None

def print_cpca(mpca):
    """
    PrintCPCA: Print to video the CPCA Model
    """
    lsci.PrintCPCA(mpca)


class CPCA():
    """
    CPCA model class
    """
    def __init__(self, scaling, npc):
        self.model = new_cpca_model()
        self.scaling = scaling
        self.npc = npc

    def __del__(self):
        if self.model is not None:
            del_cpca_model(self.model)
            del self.model
        self.model = None

    def fit(self, t_input):
        """
        fit a cpca model using an input tensor
        """
        if "Tensor" in str(type(t_input)):
            cpca_algorithm(t_input.tns, self.scaling, self.npc, self.model)
        else:
            t_input_ = tns.new_tensor(t_input)
            cpca_algorithm(t_input_, self.scaling, self.npc, self.model)
            tns.del_tensor(t_input_)
            del t_input_

    def get_super_scores(self):
        """
        get the cpca super scores
        """
        return mx.matrix_to_list(self.model[0].super_scores)

    def get_super_weights(self):
        """
        get the cpca super weights
        """
        return mx.matrix_to_list(self.model[0].super_weights)

    def get_block_scores(self):
        """
        get the cpca block scores
        """
        return tns.tensor_tolist(self.model[0].block_scores)

    def get_block_loadings(self):
        """
        get the cpca block loadings
        """
        return tns.tensor_tolist(self.model[0].block_loadings)

    def get_block_expvar(self):
        """
        get the cpca block variance explained
        """
        return vlst.dvector_list_tolist(self.model[0].block_expvar)

    def get_total_exp_variance(self):
        """
        get the cpca total variance explained
        """
        return vect.dvector_tolist(self.model[0].total_expvar)

    def predict(self, t_input):
        """
        get the cpca block variance explained
        """
        t_input_ = tns.new_tensor(t_input)
        p_super_scores_ = mx.init_matrix()
        p_block_scores_ = tns.init_tensor()

        cpca_score_predictor(t_input_,
                            self.model,
                            self.npc,
                            p_super_scores_,
                            p_block_scores_)
        p_super_scores = mx.matrix_to_list(p_super_scores_)
        p_block_scores = tns.tensor_tolist(p_block_scores_)
        mx.del_matrix(p_super_scores_)
        del p_super_scores_

        tns.del_tensor(p_block_scores_)
        del p_block_scores_

        tns.del_tensor(t_input_)
        del t_input_
        return p_super_scores, p_block_scores


if __name__ == '__main__':
    def mx_to_video(m_input, decimals=5):
        """
        print a matrix to video
        """
        for row in m_input:
            print("\t".join([str(round(x, decimals)) for x in row]))

    def t_to_video(t_input):
        """
        print a tensor to video
        """
        i = 1
        for mx_input in t_input:
            print(f"Block: {i}")
            mx_to_video(mx_input, 3)
            i+=1

    import random
    random.seed(123456)
    a = [[[random.random() for j in range(4)] for i in range(10)] for k in range(4)]
    t = tns.Tensor(a)
    print("Original Tensor")
    t_to_video(a)
    print("Computing CPCA ...")
    cpca = CPCA(1, 2)
    cpca.fit(a)
    cpca2 = CPCA(1, 2)
    cpca2.fit(t)
    print("Showing the CPCA super scores")
    sscores = cpca.get_super_scores()
    mx_to_video(sscores, 3)
    print("Showing the CPCA super weights")
    sweights = cpca.get_super_weights()
    mx_to_video(sweights, 3)
    print("Showing the CPCA block scores")
    block_scores = cpca.get_block_scores()
    t_to_video(block_scores)
    print("Showing the CPCA block loadings")
    block_loadings = cpca.get_block_loadings()
    t_to_video(block_loadings)

    print(cpca.get_total_exp_variance())
    p_ss, p_bs = cpca.predict(a)
    mx_to_video(p_ss, 3)
