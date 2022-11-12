"""
pls libscientific python binding

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
import libscientific.tensor as t
from libscientific.loadlibrary import LoadLibrary

lsci = LoadLibrary()


class PLSMODEL(ctypes.Structure):
    """
    PLSMODEL data structure
    """
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

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


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
                     ctypes.POINTER(mx.matrix),
                     ctypes.c_size_t,
                     ctypes.c_size_t,
                     ctypes.c_size_t,
                     ctypes.POINTER(PLSMODEL),
                     ctypes.POINTER(ctypes.c_int)]
lsci.PLS.restype = None


def PLS_(x, y, nlv, xscaling, yscaling, mpls):
    """
    PLS: Calculate the PLS model using a matrix x and a matrix y
         according to the NIPALS algorithm

    PLS(matrix *mx,
        matrix *my,
        size_t nlv,
        size_t xautoscaling,
        size_t yautoscaling,
        PLSMODEL* model,
        ssignal *s)
    """
    print(nlv)
    ssignal = ctypes.c_int(0)
    lsci.PLS(x,
             y,
             nlv,
             xscaling,
             yscaling,
             mpls,
             ctypes.pointer(ssignal))


lsci.PLSBetasCoeff.argtypes = [ctypes.POINTER(PLSMODEL),
                               ctypes.c_size_t,
                               ctypes.POINTER(vect.dvector)]
lsci.PLSBetasCoeff.restype = None


def PLSBetasCoeff(mpls, nlv, bcoeff):
    """
    PLSBetasCoeff calculation
    """
    lsci.PLSBetasCoeff(mpls, nlv, ctypes.pointer(bcoeff))



lsci.PLSScorePredictor.argtypes = [ctypes.POINTER(mx.matrix),
                                   ctypes.POINTER(PLSMODEL),
                                   ctypes.c_size_t,
                                   ctypes.POINTER(mx.matrix)]
lsci.PLSScorePredictor.restype = None


def PLSScorePredictor(x, mpls, nlv, pscores):
    """
    PLSScorePredictor: Predict scores for a matrix m in the computed PLS modes
    """
    lsci.PLSScorePredictor(x,
                           mpls,
                           nlv,
                           pscores)


# void PLSYPredictor(matrix *tscore, PLSMODEL *model, size_t nlv, matrix **y);

lsci.PLSYPredictor.argtypes = [ctypes.POINTER(mx.matrix),
                               ctypes.POINTER(PLSMODEL),
                               ctypes.c_size_t,
                               ctypes.POINTER(mx.matrix)]
lsci.PLSYPredictor.restype = None


def PLSYPredictor(xscores, mpls, nlv, predicted_y):
    """
    PLSYPredictor: Predict the Y according the predicted scores and the calculated pls model.
                   This function is dependent on PLSScorePredictor.
    """
    lsci.PLSYPredictor(xscores, mpls, nlv, ctypes.pointer(ctypes.pointer(predicted_y)))


# void PLSYPredictorAllLV(matrix *mx, PLSMODEL *model, matrix **tscores, matrix **y);

lsci.PLSYPredictorAllLV.argtypes = [ctypes.POINTER(mx.matrix),
                                    ctypes.POINTER(PLSMODEL),
                                    ctypes.POINTER(mx.matrix),
                                    ctypes.POINTER(mx.matrix)]
lsci.PLSYPredictorAllLV.restype = None


def PLSYPredictorAllLV(x, mpls, predicted_scores, predicted_y):
    """
    PLSYPredictorAllLV: Predict the Y according the original
                        feature matrix and the calculated pls model.
                        This function is NOT dependent on PLSScorePredictor.
    """
    lsci.PLSYPredictorAllLV(x,
                            mpls,
                            predicted_scores,
                            predicted_y)



lsci.PrintPLSModel.argtypes = [ctypes.POINTER(PLSMODEL)]
lsci.PrintPLSModel.restype = None


def PrintPLS(mpls):
    """
    PrintPLS: Print to video the PLS Model
    """
    lsci.PrintPLSModel(mpls)


class PLS(object):
    """
    Partial least squares
        Arguments:
            nlv (int) : number of latent variable
            xscaling (int) : scaling type for x matrix
            yscaling (int) : scaling type for y matrix
            N.B.: xscaling and yscaling can be:
                0 -> No scaling
                1 -> Standard deviation scaling
                2 -> Root mean squared scaling
                3 -> Pareto scaling
                4 -> Range scaling
                5 -> Level scaling
    """
    def __init__(self, nlv, xscaling=1, yscaling=0):
        self.mpls = NewPLSModel()
        self.nlv = nlv
        self.xscaling = xscaling
        self.yscaling = yscaling

    def __del__(self):
        if self.mpls is not None:
            DelPLSModel(self.mpls)
            del self.mpls
        self.mpls = None

    def fit(self, x_, y_, cross_validation=None):
        """
        Fit a pls model giving x matrix and y matrix
        """
        x = None
        xalloc = False
        if "Matrix" not in str(type(x_)):
            x = mx.NewMatrix(x_)
            xalloc = True
        else:
            x = x_

        y = None
        yalloc = False
        if "Matrix" not in str(type(y_)):
            y = mx.NewMatrix(y_)
            yalloc = True
        else:
            y = y_
        PLS_(x,
             y,
             self.nlv,
             self.xscaling,
             self.yscaling,
             self.mpls)

        if xalloc is True:
            mx.DelMatrix(x)
            del x

        if yalloc is True:
            mx.DelMatrix(y)
            del y

    def get_tscores(self):
        """
        Get the T-Scores
        """
        return mx.MatrixToList(self.mpls[0].xscores)

    def get_uscores(self):
        """
        Get the U-Scores
        """
        return mx.MatrixToList(self.mpls[0].yscores)

    def get_ploadings(self):
        """
        Get the P-Loadings
        """
        return mx.MatrixToList(self.mpls[0].xloadings)

    def get_qloadings(self):
        """
        Get the Q-Loadings
        """
        return mx.MatrixToList(self.mpls[0].yloadings)

    def get_weights(self):
        """
        Get the W-weigths
        """
        return mx.MatrixToList(self.mpls[0].xweights)

    def get_exp_variance(self):
        """
        Get the explained variance
        """
        return vect.DVectorToList(self.mpls[0].xvarexp)

    def predict(self, x_, nlv_=None):
        """
        Predict the y giving an x_ matrix
        """
        x = mx.NewMatrix(x_)
        pscores_ = mx.initMatrix()
        py_ = mx.initMatrix()

        nlv = None
        if nlv_ is None:
            nlv = self.nlv
        else:
            nlv = self.nlv

        PLSYPredictorAllLV(x, self.mpls, pscores_, py_)
        pscores = mx.MatrixToList(pscores_)
        py = mx.MatrixToList(py_)
        mx.DelMatrix(x)
        del x
        mx.DelMatrix(pscores_)
        del pscores_
        mx.DelMatrix(py_)
        del py_
        return py, pscores


if __name__ == '__main__':
    def mx_to_video(m, decimals=5):
        for row in m:
            print("\t".join([str(round(x, decimals)) for x in row]))
    import random
    random.seed(123456)
    x = [[random.random() for j in range(4)] for i in range(10)]
    y = [[random.random() for j in range(1)] for i in range(10)]
    xp = [[random.random() for j in range(4)] for i in range(10)]

    print("Original Matrix")
    print("X")
    mx_to_video(x)
    print("Y")
    mx_to_video(y)
    print("XP")
    mx_to_video(xp)
    print("Computing PLS ...")
    model = PLS(nlv=2, xscaling=1, yscaling=0)
    model.fit(x, y)
    print("Showing the PLS T scores")
    tscores = model.get_tscores()
    mx_to_video(tscores, 3)

    print("Showing the PLS U scores")
    uscores = model.get_uscores()
    mx_to_video(uscores, 3)

    print("Showing the PLS P loadings")
    ploadings = model.get_ploadings()
    mx_to_video(ploadings, 3)

    print("Showing the X Variance")
    print(model.get_exp_variance())


    print("Predict XP")
    py, pscores = model.predict(xp)
    print("Predicted Y for all LVs")
    mx_to_video(py, 3)
    print("Predicted Scores")
    mx_to_video(pscores, 3)
