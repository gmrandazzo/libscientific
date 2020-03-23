# pls libscientific python binding
#
# Copyright (C) <2019>  Giuseppe Marco Randazzo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import ctypes
import scientific.matrix as mx
import scientific.vector as vect
import scientific.tensor as t
from scientific.loadlibrary import LoadLibrary

lsci = LoadLibrary()


class PLSMODEL(ctypes.Structure):
    _fields_ = [
        ("xscores", ctypes.POINTER(mx.matrix)),
        ("xloadings", ctypes.POINTER(mx.matrix)),
        ("xweights", ctypes.POINTER(mx.matrix)),
        ("yscores", ctypes.POINTER(mx.matrix)),
        ("yloadings", ctypes.POINTER(mx.matrix)),
        ("cweights", ctypes.POINTER(mx.matrix)),
        ("b", ctypes.POINTER(vect.dvector)),
        ("xvarexp", ctypes.POINTER(vect.dvector)),
        ("xcolaverage", ctypes.POINTER(vect.dvector)),
        ("xcolscaling", ctypes.POINTER(vect.dvector)),
        ("ycolaverage", ctypes.POINTER(vect.dvector)),
        ("ycolscaling", ctypes.POINTER(vect.dvector)),
        ("recalculated_y", ctypes.POINTER(mx.matrix)),
        ("recalc_residuals", ctypes.POINTER(mx.matrix)),
        ("predicted_y", ctypes.POINTER(mx.matrix)),
        ("pred_residuals", ctypes.POINTER(mx.matrix)),
        ("recalc_residuals", ctypes.POINTER(mx.matrix)),
        ("r2y_recalculated", ctypes.POINTER(mx.matrix)),
        ("r2y_validation", ctypes.POINTER(mx.matrix)),
        ("q2y", ctypes.POINTER(mx.matrix)),
        ("sdep", ctypes.POINTER(mx.matrix)),
        ("sdec", ctypes.POINTER(mx.matrix)),
        ("bias", ctypes.POINTER(mx.matrix)),
        ("roc_recalculated", ctypes.POINTER(t.tensor)),
        ("roc_validation", ctypes.POINTER(t.tensor)),
        ("roc_auc_recalculated", ctypes.POINTER(mx.matrix)),
        ("roc_auc_validation", ctypes.POINTER(mx.matrix)),
        ("roc_auc_validation", ctypes.POINTER(mx.matrix)),
        ("precision_recall_recalculated", ctypes.POINTER(t.tensor)),
        ("precision_recall_validation", ctypes.POINTER(t.tensor)),
        ("precision_recall_ap_recalculated", ctypes.POINTER(mx.matrix)),
        ("precision_recall_ap_validation", ctypes.POINTER(mx.matrix)),
        ("yscrambling", ctypes.POINTER(mx.matrix))]


lsci.NewPLSModel.argtypes = [ctypes.POINTER(ctypes.POINTER(PLSMODEL))]
lsci.NewPLSModel.restype = None


def NewPLSModel():
    """
    NewPLSModel: Allocate in memory an empty libscientific PLS model
    """
    mpls = ctypes.POINTER(PLSMODEL)()
    lsci.NewPLSModel(ctypes.pointer(mpls))
    return mpls


lsci.DelPLSModel.argtypes = [ctypes.POINTER(ctypes.POINTER(PLSMODEL))]
lsci.DelPLSModel.restype = None


def DelPLSModel(mpls):
    """
    DelPLSModel: Delete an allocated libscientific PLS model
    """
    lsci.DelPLSModel(ctypes.pointer(mpls))


lsci.PLS.argtypes = [ctypes.POINTER(mx.matrix),
                     ctypes.c_size_t,
                     ctypes.c_size_t,
                     ctypes.POINTER(PCAMODEL),
                     ctypes.c_int]
lsci.PCA.restype = None


def PCA_(m, scaling, npc, mpca):
    """
    PCA: Calculate the PCA model for a matrix m using the NIPALS algorithm
    """
    lsci.PCA(m, scaling, npc, mpca, 0)


lsci.PCAScorePredictor.argtypes = [ctypes.POINTER(mx.matrix),
                                   ctypes.POINTER(PCAMODEL),
                                   ctypes.c_size_t,
                                   ctypes.POINTER(ctypes.POINTER(mx.matrix))]
lsci.PCAScorePredictor.restype = None


def PCAScorePredictor(m, mpca, npc, pscores):
    """
    PCAScorePredictor: Predict scores for a matrix m in the computed PCA modes
    """
    lsci.PCAScorePredictor(m,
                           mpca,
                           npc,
                           ctypes.pointer(ctypes.pointer(m)))


lsci.PCAIndVarPredictor.argtypes = [ctypes.POINTER(mx.matrix),
                                    ctypes.POINTER(mx.matrix),
                                    ctypes.POINTER(vect.dvector),
                                    ctypes.POINTER(vect.dvector),
                                    ctypes.c_size_t,
                                    ctypes.POINTER(ctypes.POINTER(mx.matrix))]
lsci.PCAIndVarPredictor.restype = None


def PCAIndVarPredictor(t, p, colaverage, colscaling, npc, iv):
    """
    PCAIndVarPredictor: Reconstruct the original matrix from PCA model,
                        scores and loadings
    """
    lsci.PCAIndVarPredictor(t, p, colaverage, colscaling, npc, ctypes.pointer(iv))


lsci.PrintPCA.argtypes = [ctypes.POINTER(PCAMODEL)]
lsci.PrintPCA.restype = None


def PrintPCA(mpca):
    """
    PrintPCA: Print to video the PCA Model
    """
    lsci.PrintPCA(mpca)


class PCA(object):
    def __init__(self, scaling, npc):
        self.mpca = NewPCAModel()
        self.scaling = scaling
        self.npc = npc

    def __del__(self):
        if self.mpca is not None:
            DelPCAModel(self.mpca)
            del self.mpca
        self.mpca = None

    def fit(self, m_):
        if "Matrix" in str(type(m_)):
            PCA_(m_, self.scaling, self.npc, self.mpca)
        else:
            m = mx.NewMatrix(m_)
            PCA_(m, self.scaling, self.npc, self.mpca)
            mx.DelMatrix(m)
            del m

    def get_scores(self):
        return mx.MatrixToList(self.mpca[0].scores)

    def get_loadings(self):
        return mx.MatrixToList(self.mpca[0].loadings)

    def get_exp_variance(self):
        return vect.DVectorToList(self.mpca[0].varexp)

    def predict(self, m_):
        m = mx.NewMatrix(m_)
        pscores_ = mx.initMatrix()
        PCAScorePredictor(m, self.mpca, self.npc, ctypes.pointer(pscores_))
        pscores = pscores_.tolist()
        mx.DelMatrix(m)
        del m
        mx.DelMatrix(pscores_)
        del pscores_
        return pscores

    def reconstruct_original_matrix(self,
                                    npc_=None,
                                    scores_=None):
        npc = None
        if npc_ is None:
            npc = self.npc
        else:
            npc = npc_

        scores = None
        if scores_ is None:
            scores = self.mpca[0].scores
        else:
            scores = scores_

        ivals = mx.initMatrix()
        PCAIndVarPredictor(scores,
                           self.mpca[0].loadings,
                           self.mpca[0].colaverage,
                           self.mpca[0].colscaling,
                           npc,
                           ivals)
        omx = mx.MatrixToList(ivals)
        mx.DelMatrix(ivals)
        del ivals
        return omx


if __name__ == '__main__':
    def mx_to_video(m, decimals=5):
        for row in m:
            print("\t".join([str(round(x, decimals)) for x in row]))
    import random
    random.seed(123456)
    a = [[random.random() for j in range(4)] for i in range(10)]
    print("Original Matrix")
    mx_to_video(a)
    print("Computing PCA ...")
    model = PCA(1, 2)
    model.fit(a)
    print("Showing the PCA scores")
    scores = model.get_scores()
    mx_to_video(scores, 3)
    print("Showing the PCA loadings")
    loadings = model.get_loadings()
    mx_to_video(loadings, 3)
    print(model.get_exp_variance())
    print("Reconstruct the original PCA matrix using the PCA Model")
    ra = model.reconstruct_original_matrix()
    mx_to_video(ra)
