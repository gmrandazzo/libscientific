"""
pca libscientific python binding

Copyright (C) <2019>  Giuseppe Marco Randazzo

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
import libscientific.matrix as mx
import libscientific.vector as vect
from libscientific.loadlibrary import load_libscientific_library

lsci = load_libscientific_library()


class PCAMODEL(ctypes.Structure):
    """
    PCA model data structure
    """
    _fields_ = [
        ("scores", ctypes.POINTER(mx.MATRIX)),
        ("loadings", ctypes.POINTER(mx.MATRIX)),
        ("dmodx", ctypes.POINTER(mx.MATRIX)),
        ("varexp", ctypes.POINTER(vect.DVECTOR)),
        ("colaverage", ctypes.POINTER(vect.DVECTOR)),
        ("colscaling", ctypes.POINTER(vect.DVECTOR))]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__



lsci.NewPCAModel.argtypes = [ctypes.POINTER(ctypes.POINTER(PCAMODEL))]
lsci.NewPCAModel.restype = None


def new_pca_model():
    """
    new_pca_model: Allocate in memory an empty libscientific PCA model
    """
    mpca = ctypes.POINTER(PCAMODEL)()
    lsci.NewPCAModel(ctypes.pointer(mpca))
    return mpca


lsci.DelPCAModel.argtypes = [ctypes.POINTER(ctypes.POINTER(PCAMODEL))]
lsci.DelPCAModel.restype = None


def del_pca_model(mpca):
    """
    del_pca_model: Delete an allocated libscientific PCA model
    """
    lsci.DelPCAModel(ctypes.pointer(mpca))


lsci.PCA.argtypes = [ctypes.POINTER(mx.MATRIX),
                     ctypes.c_size_t,
                     ctypes.c_size_t,
                     ctypes.POINTER(PCAMODEL),
                     ctypes.POINTER(ctypes.c_int)]
lsci.PCA.restype = None


def pca_algorithm(m_input, scaling, npc, mpca):
    """
    PCA: Calculate the PCA model for a matrix m using the NIPALS algorithm
    """
    ssignal = ctypes.c_int(0)
    lsci.PCA(m_input, scaling, npc, mpca, ssignal)


lsci.PCAScorePredictor.argtypes = [ctypes.POINTER(mx.MATRIX),
                                   ctypes.POINTER(PCAMODEL),
                                   ctypes.c_size_t,
                                   ctypes.POINTER(mx.MATRIX)]
lsci.PCAScorePredictor.restype = None


def pca_score_predictor(m_input, mpca, npc, pscores):
    """
    pca_score_predictor: Predict scores for a matrix m in the computed PCA modes
    """
    lsci.PCAScorePredictor(m_input,
                           mpca,
                           npc,
                           pscores)


lsci.PCAIndVarPredictor.argtypes = [ctypes.POINTER(mx.MATRIX),
                                    ctypes.POINTER(mx.MATRIX),
                                    ctypes.POINTER(vect.DVECTOR),
                                    ctypes.POINTER(vect.DVECTOR),
                                    ctypes.c_size_t,
                                    ctypes.POINTER(mx.MATRIX)]
lsci.PCAIndVarPredictor.restype = None


def pca_ind_var_predictor(scores_input,
                          loadings_input,
                          colaverage,
                          colscaling,
                          npc,
                          ind_vars_out):
    """
    pca_ind_var_predictor: Reconstruct the original matrix
                           from PCA model using scores and loadings
    """
    lsci.PCAIndVarPredictor(scores_input,
                            loadings_input,
                            colaverage,
                            colscaling,
                            npc,
                            ind_vars_out)


lsci.PrintPCA.argtypes = [ctypes.POINTER(PCAMODEL)]
lsci.PrintPCA.restype = None


def print_pca(mpca):
    """
    PrintPCA: Print to video the PCA Model
    """
    lsci.PrintPCA(mpca)


class PCA():
    """
    PCA model class
    """
    def __init__(self, scaling, npc):
        self.mpca = new_pca_model()
        self.scaling = scaling
        self.npc = npc

    def __del__(self):
        if self.mpca is not None:
            del_pca_model(self.mpca)
            del self.mpca
        self.mpca = None

    def fit(self, m_input):
        """
        fit a PCA model using an input matrix
        """
        if "Matrix" in str(type(m_input)):
            pca_algorithm(m_input.mtx, self.scaling, self.npc, self.mpca)
        else:
            m_input_ = mx.new_matrix(m_input)
            pca_algorithm(m_input_, self.scaling, self.npc, self.mpca)
            mx.del_matrix(m_input_)
            del m_input_

    def get_scores(self):
        """
        get the PCA scores
        """
        return mx.matrix_to_list(self.mpca[0].scores)

    def get_loadings(self):
        """
        get the PCA loadings
        """
        return mx.matrix_to_list(self.mpca[0].loadings)

    def get_exp_variance(self):
        """
        get the PCA explained variance
        """
        return vect.dvector_tolist(self.mpca[0].varexp)

    def predict(self, m_input):
        """
        Project an input matrix into the PCA model
        """
        pscores_ = mx.init_matrix()
        if "Matrix" in str(type(m_input)):
            m_input_ = m_input
            pca_score_predictor(m_input.mtx, self.mpca, self.npc, pscores_)
        else:
            m_input_ = mx.new_matrix(m_input)
            pca_score_predictor(m_input_, self.mpca, self.npc, pscores_)
            mx.del_matrix(m_input_)
            del m_input_
        pscores = mx.matrix_to_list(pscores_)
        mx.del_matrix(pscores_)
        del pscores_
        return pscores

    def reconstruct_original_matrix(self,
                                    npc_input=None,
                                    scores_input=None):
        """
        Reconstruct the original input matrix giving a
        number of principal components to be used from scores and loadings
        """
        npc_ = None
        if npc_input is None:
            npc_ = self.npc
        else:
            npc_ = npc_input

        scores_ = None
        if scores_input is None:
            scores_ = self.mpca[0].scores
        else:
            scores_ = scores_input

        ind_vars = mx.init_matrix()
        pca_ind_var_predictor(scores_,
                              self.mpca[0].loadings,
                              self.mpca[0].colaverage,
                              self.mpca[0].colscaling,
                              npc_,
                              ind_vars)
        omx = mx.matrix_to_list(ind_vars)
        mx.del_matrix(ind_vars)
        del ind_vars
        return omx


if __name__ == '__main__':
    def mx_to_video(m_input, decimals=5):
        """
        print a matrix to video
        """
        for row in m_input:
            print("\t".join([str(round(x, decimals)) for x in row]))
    import random
    random.seed(123456)
    a = [[random.random() for j in range(4)] for i in range(10)]
    m = mx.Matrix(a)
    print("Original Matrix")
    mx_to_video(a)
    print("Computing PCA ...")
    model = PCA(1, 2)
    model.fit(a)
    model2 = PCA(1, 2)
    model2.fit(m)
    print("Showing the PCA scores")
    scores = model.get_scores()
    mx_to_video(scores, 3)
    print("Showing the PCA scores")
    scores2 = model2.get_scores()
    mx_to_video(scores2, 3)
    print("Showing the PCA loadings")
    loadings = model.get_loadings()
    mx_to_video(loadings, 3)
    print(model.get_exp_variance())
    print("Reconstruct the original PCA matrix using the PCA Model")
    ra = model.reconstruct_original_matrix()
    mx_to_video(ra)
    pred_scores = model.predict(a)
    for i, row1 in enumerate(pred_scores):
        row2 = scores[i]
        print(row1, row2)
    
