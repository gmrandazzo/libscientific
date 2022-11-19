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
from libscientific.loadlibrary import LoadLibrary

lsci = LoadLibrary()


class CPCAMODEL(ctypes.Structure):
    _fields_ = [
        ("block_scores", ctypes.POINTER(tns.tensor)),
        ("block_loadings", ctypes.POINTER(tns.tensor)),
        ("super_scores", ctypes.POINTER(mx.matrix)),
        ("super_weights", ctypes.POINTER(mx.matrix)),
        ("scaling_factor", ctypes.POINTER(vect.dvector)),
        ("total_expvar", ctypes.POINTER(vect.dvector)),
        ("block_expvar", ctypes.POINTER(vlst.dvectlist)),
        ("colaverage", ctypes.POINTER(vlst.dvectlist)),
        ("colscaling", ctypes.POINTER(vlst.dvectlist))]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__



lsci.NewCPCAModel.argtypes = [ctypes.POINTER(ctypes.POINTER(CPCAMODEL))]
lsci.NewCPCAModel.restype = None


def NewCPCAModel():
    """
    NewCPCAModel: Allocate in memory an empty libscientific CPCA model
    """
    mpca = ctypes.POINTER(CPCAMODEL)()
    lsci.NewCPCAModel(ctypes.pointer(mpca))
    return mpca


lsci.DelCPCAModel.argtypes = [ctypes.POINTER(ctypes.POINTER(CPCAMODEL))]
lsci.DelCPCAModel.restype = None


def DelCPCAModel(mpca):
    """
    DelCPCAModel: Delete an allocated libscientific CPCA model
    """
    lsci.DelCPCAModel(ctypes.pointer(mpca))


lsci.CPCA.argtypes = [ctypes.POINTER(tns.tensor),
                     ctypes.c_size_t,
                     ctypes.c_size_t,
                     ctypes.POINTER(CPCAMODEL),
                     ctypes.POINTER(ctypes.c_int)]
lsci.CPCA.restype = None


def CPCA_(m, scaling, npc, mpca):
    """
    CPCA: Calculate the CPCA model for a matrix m using the NIPALS algorithm
    """
    ssignal = ctypes.c_int(0)
    lsci.CPCA(m, scaling, npc, mpca, ssignal)


lsci.CPCAScorePredictor.argtypes = [ctypes.POINTER(tns.tensor),
                                   ctypes.POINTER(CPCAMODEL),
                                   ctypes.c_size_t,
                                   ctypes.POINTER(mx.matrix),
                                   ctypes.POINTER(tns.tensor)]
lsci.CPCAScorePredictor.restype = None


def CPCAScorePredictor(t, model, npc, p_super_scores, p_block_scores):
    """
    CPCAScorePredictor: Predict scores for a matrix m in the computed CPCA modes
    """
    lsci.CPCAScorePredictor(t,
                            model,
                            npc,
                            p_super_scores,
                            p_block_scores)


lsci.PrintCPCA.argtypes = [ctypes.POINTER(CPCAMODEL)]
lsci.PrintCPCA.restype = None

def PrintCPCA(mpca):
    """
    PrintCPCA: Print to video the CPCA Model
    """
    lsci.PrintCPCA(mpca)


class CPCA(object):
    def __init__(self, scaling, npc):
        self.model = NewCPCAModel()
        self.scaling = scaling
        self.npc = npc

    def __del__(self):
        if self.model is not None:
            DelCPCAModel(self.model)
            del self.model
        self.model = None

    def fit(self, t_):
        if "Tensor" in str(type(t_)):
            CPCA_(t_, self.scaling, self.npc, self.model)
        else:
            t = tns.NewTensor(t_)
            CPCA_(t, self.scaling, self.npc, self.model)
            tns.DelTensor(t)
            del t

    def get_super_scores(self):
        return mx.MatrixToList(self.model[0].super_scores)

    def get_super_weights(self):
        return mx.MatrixToList(self.model[0].super_weights)

    def get_block_scores(self):
        return tns.TensorToList(self.model[0].block_scores)
    
    def get_block_loadings(self):
        return tns.TensorToList(self.model[0].block_loadings)

    def get_block_expvar(self):
        return vlst.DVectorListToList(self.model[0].block_expvar)
    
    def get_total_exp_variance(self):
        return vect.DVectorToList(self.model[0].total_expvar)

    def predict(self, t_):
        t = tns.NewTensor(t_)
        p_super_scores_ = mx.initMatrix()
        p_block_scores_ = tns.initTensor()
        
        CPCAScorePredictor(t,
                           self.model,
                           self.npc,
                           p_super_scores_,
                           p_block_scores_)
        p_super_scores = mx.MatrixToList(p_super_scores_)
        p_block_scores = tns.TensorToList(p_block_scores_)
        mx.DelMatrix(p_super_scores_)
        del p_super_scores_
        
        tns.DelTensor(p_block_scores_)
        del p_block_scores_
        
        tns.DelTensor(t)
        del t
        return p_super_scores, p_block_scores


if __name__ == '__main__':
    def mx_to_video(m, decimals=5):
        for row in m:
            print("\t".join([str(round(x, decimals)) for x in row]))

    def t_to_video(t):
        i = 1
        for m in t:
            print("Block: %d" % (i))
            mx_to_video(m, 3)
            i+=1
    
    import random
    random.seed(123456)
    a = [[[random.random() for j in range(4)] for i in range(10)] for k in range(4)]
    t = tns.NewTensor(a)
    print("Original Tensor")
    t_to_video(a)
    print("Computing CPCA ...")
    model = CPCA(1, 2)
    model.fit(a)
    print("Showing the CPCA super scores")
    sscores = model.get_super_scores()
    mx_to_video(sscores, 3)
    print("Showing the CPCA super weights")
    sweights = model.get_super_weights()
    mx_to_video(sweights, 3)
    print("Showing the CPCA block scores")
    block_scores = model.get_block_scores()
    t_to_video(block_scores)
    print("Showing the CPCA block loadings")
    block_loadings = model.get_block_loadings()
    t_to_video(block_loadings)

    print(model.get_total_exp_variance())
    p_ss, p_bs = model.predict(a)
    mx_to_video(p_ss, 3)
