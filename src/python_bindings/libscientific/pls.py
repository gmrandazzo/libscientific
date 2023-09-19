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
    """
    Predict using PLS beta coefficients and not the whole PLS model.

    Parameters
    ----------
    x_input : list
        Input matrix to predict.
    **kwargs : dict
        Additional keyword arguments:
        - x_averages (list): Column averages from the training x matrix.
        - x_scaling (list): Column scaling from the training x matrix.
        - y_averages (list): Column averages from the training y matrix.
        - y_scaling (list): Column scaling from the training y matrix.
        - beta_coeff (list): PLS beta coefficients extracted from the trained PLS model.

    Returns
    -------
    predicted_y : list
        Predicted values.

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


lsci.WritePLS.argtypes = [
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(PLSMODEL)
]
lsci.WritePLS.restype = None

def write_pls(dbpath: str, pls: PLSMODEL):
    """
    write_pls(dbpath, pls):
    
    Writes a PLS (Partial Least Squares) model to
    a SQLite database.

    Parameters:
        dbpath (str): The path to the SQLite database file where the PLS 
            model will be stored.
        pls (PLSMODEL): The libscientific data structure representing the PLS 
            model.

    Returns:
        None
    """
    encoded_string = dbpath.encode('utf-8')
    ctypes_string = ctypes.c_char_p(encoded_string)
    lsci.WritePLS(ctypes_string, pls)


lsci.ReadPLS.argtypes = [
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(PLSMODEL)
]
lsci.ReadPLS.restype = None

def read_pls(dbpath: str, pls: PLSMODEL):
    """
    read_pls(dbpath, pls):
    
    Reads a PLS (Partial Least Squares) model from a SQLite database.

    Parameters:
        dbpath (str): The path to the SQLite database file from which the PLS
            model will be read.
        pls (PLSMODEL): The libscientific data structure where the PLS
            model will be loaded.

    Returns:
        None
    """
    encoded_string = dbpath.encode('utf-8')
    ctypes_string = ctypes.c_char_p(encoded_string)
    lsci.ReadPLS(ctypes_string, pls)

class PLS():
    """
    Partial Least Squares (PLS) Model

    This class provides methods for creating and utilizing a PLS model for regression analysis.

    Parameters
    ----------
    nlv (int) : Number of latent variables.
    xscaling (int) : Scaling type for x matrix. Default is 1.
    yscaling (int) : Scaling type for y matrix. Default is 0.

    Note
    ----
    Scaling options:
    - 0: No scaling
    - 1: Standard deviation scaling
    - 2: Root mean squared scaling
    - 3: Pareto scaling
    - 4: Range scaling
    - 5: Level scaling

    Attributes
    ----------
    mpls (CDataType) : The PLS model data.
    nlv (int) : Number of latent variables.
    xscaling (int) : Scaling type for x matrix.
    yscaling (int) : Scaling type for y matrix.

    Methods
    -------
    __init__(self, nlv, xscaling=1, yscaling=0)
    __del__(self)
    fit(self, x_input, y_input)
    get_tscores(self)
    get_uscores(self)
    get_ploadings(self)
    get_qloadings(self)
    get_weights(self)
    get_beta_coefficients(self, nlv : int=1)
    get_exp_variance(self)
    get_x_column_scaling(self)
    get_x_averages(self)
    get_y_column_scaling(self)
    get_y_averages(self)
    predict(self, x_input, nlv_=None)
    """

    def __init__(self, nlv=2, xscaling=1, yscaling=0):
        """
        Initialize a PLS instance.

        Parameters
        ----------
        nlv (int) : Number of latent variables.
        xscaling (int) : Scaling type for x matrix. Default is 1.
        yscaling (int) : Scaling type for y matrix. Default is 0.
        """
        self.model = new_pls_model()
        self.nlv = nlv
        self.xscaling = xscaling
        self.yscaling = yscaling

    def __del__(self):
        """
        Clean up resources associated with the PLS instance.
        """
        if self.model is not None:
            del_pls_model(self.model)
            del self.model
        self.model = None

    def fit(self, x_input, y_input):
        """
        Fit a PLS model using x and y matrices.

        Parameters
        ----------
        x_input : Matrix or List[List]
            Input x matrix for fitting the PLS model.
        y_input : Matrix or List[List]
            Input y matrix for fitting the PLS model.
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
                      self.model,
                      nlv=self.nlv,
                      x_scaling=self.xscaling,
                      y_scaling=self.yscaling)

        if xalloc is True:
            mx.del_matrix(x_input_)
            del x_input_

        if yalloc is True:
            mx.del_matrix(y_input_)
            del y_input_

    def get_tscores(self):
        """
        Get the T-Scores.

        Returns
        -------
        List[List[float]]
            The T-Scores.
        """
        return mx.matrix_to_list(self.model.contents.xscores)

    def get_uscores(self):
        """
        Get the U-Scores.

        Returns
        -------
        List[List[float]]
            The U-Scores.
        """
        return mx.matrix_to_list(self.model.contents.yscores)

    def get_ploadings(self):
        """
        Get the P-Loadings.

        Returns
        -------
        List[List[float]]
            The P-Loadings.
        """
        return mx.matrix_to_list(self.model.contents.xloadings)

    def get_qloadings(self):
        """
        Get the Q-Loadings.

        Returns
        -------
        List[List[float]]
            The Q-Loadings.
        """
        return mx.matrix_to_list(self.model.contents.yloadings)

    def get_weights(self):
        """
        Get the W-Weights.

        Returns
        -------
        List[List[float]]
            The W-Weights.
        """
        return mx.matrix_to_list(self.model.contents.xweights)

    def get_beta_coefficients(self, nlv: int = 1):
        """
        Get the Beta Coefficients for fast predictions.

        Parameters
        ----------
        nlv : int, optional
            Number of latent variables. Default is 1.

        Returns
        -------
        List[float]
            The Beta Coefficients.
        """
        beta_coeff = vect.init_dvector()
        pls_beta_coefficients(self.model, nlv, beta_coeff)
        beta_coeff_lst = vect.dvector_tolist(beta_coeff)
        vect.del_dvector(beta_coeff)
        return beta_coeff_lst

    def get_exp_variance(self):
        """
        Get the explained variance.

        Returns
        -------
        List[float]
            The explained variance.
        """
        return vect.dvector_tolist(self.model.contents.xvarexp)

    def get_x_column_scaling(self):
        """
        Get the model column scaling for x.

        Returns
        -------
        List[float]
            The model column scaling for x.
        """
        return vect.dvector_tolist(self.model.contents.xcolscaling)

    def get_x_averages(self):
        """
        Get the feature averages for x.

        Returns
        -------
        List[float]
            The feature averages for x.
        """
        return vect.dvector_tolist(self.model.contents.xcolaverage)

    def get_y_column_scaling(self):
        """
        Get the model column scaling for y.

        Returns
        -------
        List[float]
            The model column scaling for y.
        """
        return vect.dvector_tolist(self.model.contents.ycolscaling)

    def get_y_averages(self):
        """
        Get the feature averages for y.

        Returns
        -------
        List[float]
            The feature averages for y.
        """
        return vect.dvector_tolist(self.model.contents.ycolaverage)

    def predict(self, x_input, nlv_=None):
        """
        Predict the y values given an x matrix.

        Parameters
        ----------
        x_input : Matrix or List[List]
            Input x matrix for prediction.
        nlv_ : int, optional
            Number of latent variables for prediction. Default is None.

        Returns
        -------
        List[List[float]], List[List[float]]
            The predicted y values and T-Scores.
        """
        x_input_ = mx.new_matrix(x_input)
        p_scores_ = mx.init_matrix()
        p_y_ = mx.init_matrix()

        pls_y_predictor_all_lv(x_input_, self.model, p_scores_, p_y_)
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

    def save(self, dbpath):
        """
        Save PLS model to a sqlite3 file

        Parameters
        ----------
        dbpath (str): The path to the SQLite database file where the PLS 
            model will be stored.
        """
        write_pls(dbpath, self.model)

    def load(self, dbpath):
        """
        Load PLS model from a sqlite3 file

        Parameters
        ----------
        dbpath (str): The path to the SQLite database file where the PLS 
            model will be stored.
        """
        read_pls(dbpath, self.model)
        self.nlv = self.model.contents.xscores[0].col
