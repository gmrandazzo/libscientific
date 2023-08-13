"""pls.py libscientific python binding

Copyright (C) <2023>  Giuseppe Marco Randazzo

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
        return self.__class__.__name__


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


def pls_algorithm(x_input,
                  y_input,
                  mpls,
                  **kwargs):
    """
    PLS: Calculate the PLS model using a matrix x and a matrix y according to the NIPALS algorithm
    """
    nlv = kwargs["nlv"]
    x_scaling = kwargs["x_scaling"]
    y_scaling = kwargs["y_scaling"]
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


lsci.PLSYPredictorAllLV.argtypes = [ctypes.POINTER(mx.MATRIX),
                                    ctypes.POINTER(PLSMODEL),
                                    ctypes.POINTER(mx.MATRIX),
                                    ctypes.POINTER(mx.MATRIX)]
lsci.PLSYPredictorAllLV.restype = None


def pls_y_predictor_all_lv(x_input,
                           mpls,
                           predicted_scores,
                           predicted_y):
    """
    PLSYPredictorAllLV: Predict the Y according the original feature 
    matrix and the calculated pls model. 
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

def apply_scaling(x_input, x_avg, x_scal):
    """Apply scaling to a matrix"""
    x_scaled = []
    for _, x_row in enumerate(x_input):
        x_scaled.append([])
        for j, x_val in enumerate(x_row):
            x_scaled[-1].append(((x_val - x_avg[j])/x_scal[j]))
    return x_scaled

def revert_scaling(x_input, x_avg, x_scal):
    """Revert scaling to a matrix"""
    x_revert = []
    for _, x_row in enumerate(x_input):
        x_rever_row = []
        if isinstance(x_row, list):
            for j, x_val in enumerate(x_row):
                x_rever_row.append(x_val*x_scal[j]+x_avg[j])
        else:
            x_rever_row.append(x_row*x_scal[0]+x_avg[0])
        x_revert.append(x_rever_row)
    return x_revert

def pls_beta_inference(x_input,
                       **kwargs):
    """Predict using PLS beta coefficients and not the whole PLS model

    Parameters
    ----------
    x_input (list):
        input matrix to predict
    x_averages (list):
        column averages coming from the training x matrix
    x_scaling (list)
        column scaling coming from the training x matrix 
    y_averages (list):
        column averages coming from the training y matrix
    y_scaling (list)
        column scaling coming from the training y matrix 
    beta_coeff (list):
        pls beta coefficients extracted from the trained pls model

    Returns
    -------
    predicted_y (list):
        Predicted values

    Examples
    --------
    >>> model = PLS(nlv=2, xscaling=1, yscaling=0)
    >>> model.fit(x, y)
    >>> pby = pls_beta_inference(xp,
                                 x_averages=model.get_x_averages(),
                                 x_scaling=model.get_x_column_scaling(),
                                 y_averages=model.get_y_averages(),
                                 y_scaling=model.get_y_column_scaling(),
                                 beta_coeff=model.get_beta_coefficients(2))
    """

    # Predict using beta coefficients
    x_avg = kwargs["x_averages"]
    x_scal = kwargs["x_scaling"]
    y_avg = kwargs["y_averages"]
    y_scal = kwargs["y_scaling"]
    b_coeff = kwargs["beta_coeff"]
    x_ps = mx.new_matrix(apply_scaling(x_input, x_avg, x_scal))
    pby_ = vect.new_dvector([0 for i in range(len(x_input))])
    b_coeff_ = vect.new_dvector(b_coeff)
    mx.matrix_dvector_dot_product(x_ps, b_coeff_, pby_)
    pby = vect.dvector_tolist(pby_)
    vect.del_dvector(pby_)
    vect.del_dvector(b_coeff_)
    mx.del_matrix(x_ps)
    return revert_scaling(pby, y_avg, y_scal)



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
                      self.mpls,
                      nlv=self.nlv,
                      x_scaling = self.xscaling,
                      y_scaling =self.yscaling)

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
        beta_coeff = vect.init_dvector()
        pls_beta_coefficients(self.mpls, nlv, beta_coeff)
        beta_coeff_lst = vect.dvector_tolist(beta_coeff)
        vect.del_dvector(beta_coeff)
        return beta_coeff_lst

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
