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
import libscientific.matrix as mx
import libscientific.vector as vect
import libscientific.tensor as t
from libscientific.loadlibrary import load_libscientific_library

lsci = load_libscientific_library()


class PLSMODEL(ctypes.Structure):
    """
    PLSMODEL data structure
    """
    _fields_ = [
        ("xscores", ctypes.POINTER(mx.MATRIX)),
        ("xloadings", ctypes.POINTER(mx.MATRIX)),
        ("xweights", ctypes.POINTER(mx.MATRIX)),
        ("yscores", ctypes.POINTER(mx.MATRIX)),
        ("yloadings", ctypes.POINTER(mx.MATRIX)),
        ("b", ctypes.POINTER(vect.DVECTOR)),
        ("xvarexp", ctypes.POINTER(vect.DVECTOR)),
        ("xcolaverage", ctypes.POINTER(vect.DVECTOR)),
        ("xcolscaling", ctypes.POINTER(vect.DVECTOR)),
        ("ycolaverage", ctypes.POINTER(vect.DVECTOR)),
        ("ycolscaling", ctypes.POINTER(vect.DVECTOR)),
        ("recalculated_y", ctypes.POINTER(mx.MATRIX)),
        ("recalc_residuals", ctypes.POINTER(mx.MATRIX)),
        ("predicted_y", ctypes.POINTER(mx.MATRIX)),
        ("pred_residuals", ctypes.POINTER(mx.MATRIX)),
        ("r2y_recalculated", ctypes.POINTER(mx.MATRIX)),
        ("r2y_validation", ctypes.POINTER(mx.MATRIX)),
        ("q2y", ctypes.POINTER(mx.MATRIX)),
        ("sdep", ctypes.POINTER(mx.MATRIX)),
        ("sdec", ctypes.POINTER(mx.MATRIX)),
        ("bias", ctypes.POINTER(mx.MATRIX)),
        ("roc_recalculated", ctypes.POINTER(t.TENSOR)),
        ("roc_validation", ctypes.POINTER(t.TENSOR)),
        ("roc_auc_recalculated", ctypes.POINTER(mx.MATRIX)),
        ("roc_auc_validation", ctypes.POINTER(mx.MATRIX)),
        ("precision_recall_recalculated", ctypes.POINTER(t.TENSOR)),
        ("precision_recall_validation", ctypes.POINTER(t.TENSOR)),
        ("precision_recall_ap_recalculated", ctypes.POINTER(mx.MATRIX)),
        ("precision_recall_ap_validation", ctypes.POINTER(mx.MATRIX)),
        ("yscrambling", ctypes.POINTER(mx.MATRIX))]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__nam


lsci.NewPLSModel.argtypes = [ctypes.POINTER(ctypes.POINTER(PLSMODEL))]
lsci.NewPLSModel.restype = None


def new_pls_model():
    """
    new_pls_model: Allocate in memory an empty libscientific PLS model
    """
    mpls = ctypes.POINTER(PLSMODEL)()
    lsci.NewPLSModel(ctypes.pointer(mpls))
    return mpls


lsci.DelPLSModel.argtypes = [ctypes.POINTER(ctypes.POINTER(PLSMODEL))]
lsci.DelPLSModel.restype = None


def del_pls_model(mpls):
    """
    del_pls_model: Delete an allocated libscientific PLS model
    """
    lsci.DelPLSModel(ctypes.pointer(mpls))


lsci.PLS.argtypes = [ctypes.POINTER(mx.MATRIX),
                     ctypes.POINTER(mx.MATRIX),
                     ctypes.c_size_t,
                     ctypes.c_size_t,
                     ctypes.c_size_t,
                     ctypes.POINTER(PLSMODEL),
                     ctypes.POINTER(ctypes.c_int)]
lsci.PLS.restype = None


def pls_algorithm(x_input, y_input, nlv, x_scaling, y_scaling, mpls):
    """
    PLS: Calculate the PLS model using a matrix x and a matrix y according to the NIPALS algorithm
    """
    ssignal = ctypes.c_int(0)
    lsci.PLS(x_input,
             y_input,
             nlv,
             x_scaling,
             y_scaling,
             mpls,
             ctypes.pointer(ssignal))

lsci.PLSBetasCoeff.argtypes = [ctypes.POINTER(PLSMODEL),
                               ctypes.c_size_t,
                               ctypes.POINTER(vect.DVECTOR)]
lsci.PLSBetasCoeff.restype = None


def pls_beta_coefficients(mpls, nlv, bcoeff):
    """
    PLSBetasCoeff calculation
    """
    lsci.PLSBetasCoeff(mpls, nlv, bcoeff)



lsci.PLSScorePredictor.argtypes = [ctypes.POINTER(mx.MATRIX),
                                   ctypes.POINTER(PLSMODEL),
                                   ctypes.c_size_t,
                                   ctypes.POINTER(mx.MATRIX)]
lsci.PLSScorePredictor.restype = None


def pls_score_predictor(x_input, mpls, nlv, p_scores):
    """
    PLSScorePredictor: Predict scores for a matrix m in the computed PLS modes
    """
    lsci.PLSScorePredictor(x_input,
                           mpls,
                           nlv,
                           p_scores)


lsci.PLSYPredictor.argtypes = [ctypes.POINTER(mx.MATRIX),
                               ctypes.POINTER(PLSMODEL),
                               ctypes.c_size_t,
                               ctypes.POINTER(mx.MATRIX)]
lsci.PLSYPredictor.restype = None


def pls_y_predictor(x_scores, mpls, nlv, predicted_y):
    """
    PLSYPredictor: Predict the Y according the predicted scores and the calculated pls model.
    """
    lsci.PLSYPredictor(x_scores,
                       mpls,
                       nlv,
                       ctypes.pointer(ctypes.pointer(predicted_y)))


# void PLSYPredictorAllLV(matrix *mx, PLSMODEL *model, matrix **tscores, matrix **y);

lsci.PLSYPredictorAllLV.argtypes = [ctypes.POINTER(mx.MATRIX),
                                    ctypes.POINTER(PLSMODEL),
                                    ctypes.POINTER(mx.MATRIX),
                                    ctypes.POINTER(mx.MATRIX)]
lsci.PLSYPredictorAllLV.restype = None


def pls_y_predictor_all_lv(x_input, mpls, predicted_scores, predicted_y):
    """
    PLSYPredictorAllLV: Predict the Y according the original feature matrix and the calculated pls model. 
    This function is NOT dependent on PLSScorePredictor.
    """
    lsci.PLSYPredictorAllLV(x_input,
                            mpls,
                            predicted_scores,
                            predicted_y)



lsci.PrintPLSModel.argtypes = [ctypes.POINTER(PLSMODEL)]
lsci.PrintPLSModel.restype = None


def print_pls(mpls):
    """
    PrintPLS: Print to video the PLS Model
    """
    lsci.PrintPLSModel(mpls)


class PLS():
    """
    Partial least squares

    Parameters
    ----------
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
        self.mpls = new_pls_model()
        self.nlv = nlv
        self.xscaling = xscaling
        self.yscaling = yscaling

    def __del__(self):
        if self.mpls is not None:
            del_pls_model(self.mpls)
            del self.mpls
        self.mpls = None

    def fit(self, x_input, y_input):
        """
        Fit a pls model giving x matrix and y matrix
        """
        x_input_ = None
        xalloc = False
        if "Matrix" not in str(type(x_input)):
            x_input_ = mx.new_matrix(x_input)
            xalloc = True
        else:
            x_input_ = x_input.mtx

        y_input_ = None
        yalloc = False
        if "Matrix" not in str(type(y_input)):
            y_input_ = mx.new_matrix(y_input)
            yalloc = True
        else:
            y_input_ = y_input.mtx

        pls_algorithm(x_input_,
                      y_input_,
                      self.nlv,
                      self.xscaling,
                      self.yscaling,
                      self.mpls)

        if xalloc is True:
            mx.del_matrix(x_input_)
            del x_input_

        if yalloc is True:
            mx.del_matrix(y_input_)
            del y_input_

    def get_tscores(self):
        """
        Get the T-Scores
        """
        return mx.matrix_to_list(self.mpls.contents.xscores)

    def get_uscores(self):
        """
        Get the U-Scores
        """
        return mx.matrix_to_list(self.mpls.contents.yscores)

    def get_ploadings(self):
        """
        Get the P-Loadings
        """
        return mx.matrix_to_list(self.mpls.contents.xloadings)

    def get_qloadings(self):
        """
        Get the Q-Loadings
        """
        return mx.matrix_to_list(self.mpls.contents.yloadings)

    def get_weights(self):
        """
        Get the W-weigths
        """
        return mx.matrix_to_list(self.mpls.contents.xweights)

    def get_beta_coefficients(self, nlv : int=1):
        """
        Get the Beta Coefficients for fast predictions
        """
        b = vect.init_dvector()
        pls_beta_coefficients(self.mpls, nlv, b)
        r = vect.dvector_tolist(b)
        vect.del_dvector(b)
        return r
    
    def get_exp_variance(self):
        """
        Get the explained variance
        """
        return vect.dvector_tolist(self.mpls.contents.xvarexp)

    def get_x_column_scaling(self):
        """
        Get the model column scaling
        """
        return vect.dvector_tolist(self.mpls.contents.xcolscaling)

    def get_x_averages(self):
        """
        Get the feature averages
        """
        return vect.dvector_tolist(self.mpls.contents.xcolaverage)


    def get_y_column_scaling(self):
        """
        Get the model column scaling
        """
        return vect.dvector_tolist(self.mpls.contents.ycolscaling)

    def get_y_averages(self):
        """
        Get the feature averages
        """
        return vect.dvector_tolist(self.mpls.contents.ycolaverage)

    def predict(self, x_input, nlv_=None):
        """
        Predict the y giving an x matrix
        """
        x_input_ = mx.new_matrix(x_input)
        p_scores_ = mx.init_matrix()
        p_y_ = mx.init_matrix()

        pls_y_predictor_all_lv(x_input_, self.mpls, p_scores_, p_y_)
        p_scores = mx.matrix_to_list(p_scores_)

        p_y = None
        if nlv_ is None:
            p_y = mx.matrix_to_list(p_y_)
        else:
            p_y = [[row[nlv_-1]] for row in mx.matrix_to_list(p_y_)]
        mx.del_matrix(x_input_)
        del x_input_
        mx.del_matrix(p_scores_)
        del p_scores_
        mx.del_matrix(p_y_)
        del p_y_
        return p_y, p_scores


if __name__ == '__main__':
    def mx_to_video(m_input, decimals=5):
        for row in m_input:
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
    vect.print_dvector(model.mpls.contents.xcolscaling)

    xm = mx.Matrix(x)
    ym = mx.Matrix(y)
    model2 = PLS(nlv=2, xscaling=1, yscaling=0)
    model2.fit(xm, ym)

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

    print("Beta Coefficients")
    b = model.get_beta_coefficients(2)
    print(b)
    
    print("Predict XP")
    py, pscores = model.predict(xp)
    print("Predicted Y for all LVs")
    mx_to_video(py, 3)
    print("Predicted Scores")
    mx_to_video(pscores, 3)

    print("Predict using beta coefficients")
    xa = model.get_x_averages()
    xs = model.get_x_column_scaling()
    ya = model.get_y_averages()
    ys = model.get_y_column_scaling()
    print(len(xs))
    print(xa)
    vect.print_dvector(model.mpls.contents.xcolscaling)
    print(len(xp[0]))
    xps = []
    for i in range(len(xp)):
        xps.append(list())
        for j in range(len(xp[i])):
            xps[-1].append(((xp[i][j] - xa[j])/xs[j]))
    xps_ = mx.new_matrix(xps)
    pby_ = vect.new_dvector([0 for i in range(len(xps))])
    b_ = vect.new_dvector(b)
    mx.matrix_dvector_dot_product(xps_, b_, pby_)
    pby = vect.dvector_tolist(pby_)
    vect.del_dvector(pby_)
    vect.del_dvector(b_)
    mx.del_matrix(xps_)
    for i in range(len(pby)):
        pby[i] = pby[i]*ys[0]+ya[0]
        print(f'{pby[i]} == {py[i][-1]}')
    
