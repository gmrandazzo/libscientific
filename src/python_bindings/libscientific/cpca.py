"""cpca.py libscientific python binding

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

lsci.WriteCPCA.argtypes = [
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(CPCAMODEL)
]
lsci.WriteCPCA.restype = None

def write_cpca(dbpath: str, cpca: CPCAMODEL):
    """
    write_pca(dbpath, cpca):

    Writes a CPCA (Consensus Principal Component Analysis) model to
    a SQLite database.

    Parameters:
        dbpath (str): The path to the SQLite database file where the CPCA
            model will be stored.
        cpca (CPCAMODEL): The libscientific data structure representing the CPCA
            model.

    Returns:
        None
    """
    encoded_string = dbpath.encode('utf-8')
    ctypes_string = ctypes.c_char_p(encoded_string)
    lsci.WriteCPCA(ctypes_string, cpca)


lsci.ReadCPCA.argtypes = [
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(CPCAMODEL)
]
lsci.ReadCPCA.restype = None

def read_cpca(dbpath: str, cpca: CPCAMODEL):
    """
    read_cpca(dbpath, cpca):

    Reads a CPCA (Consensus Principal Component Analysis) model from
    a SQLite database.

    Parameters:
        dbpath (str): The path to the SQLite database file from which the CPCA
            model will be read.
        cpca (CPCAMODEL): The libscientific data structure where the CPCA
            model will be loaded.

    Returns:
        None
    """
    encoded_string = dbpath.encode('utf-8')
    ctypes_string = ctypes.c_char_p(encoded_string)
    lsci.ReadCPCA(ctypes_string, cpca)

class CPCA:
    """
    CPCA Class (Consensus Principal Component Analysis)

    This class provides methods for creating and manipulating
    a Consensus Principal Component Analysis model.

    Attributes:
        model : CPCA_MODEL
            The CPCA model.
        scaling : int
            Scaling type to apply. Available scalings:
            0: No scaling. Only mean centering
            1: Mean centering and STDEV scaling
            2: Mean centering and Root-Mean-Square column scaling
            3: Mean centering and Pareto scaling
            4: Mean centering and min-max range scaling
            5: Mean centering and level scaling
        npc : int
            Number of principal components.

    Methods:
        __init__(self, scaling, npc)
        __del__(self)
        fit(self, t_input)
        get_super_scores(self)
        get_super_weights(self)
        get_block_scores(self)
        get_block_loadings(self)
        get_block_expvar(self)
        get_total_exp_variance(self)
        predict(self, t_input)
    """

    def __init__(self, scaling=1, npc=2):
        """
        Initialize a CPCA instance.

        Args:
            scaling : int
                Scaling type to apply.
            npc : int
                Number of principal components.
        """
        self.model = new_cpca_model()
        self.scaling = scaling
        self.npc = npc

    def __del__(self):
        """
        Clean up resources associated with the CPCA instance.
        """
        if self.model is not None:
            del_cpca_model(self.model)
            del self.model
        self.model = None

    def fit(self, t_input):
        """
        Fit the consensus PCA model.

        Args:
            t_input : Tensor
                The input tensor.

        Returns:
            CPCA : The fitted CPCA model instance.
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
        Get the CPCA super scores.

        Returns:
            List[List] : CPCA super scores.
        """
        return mx.matrix_to_list(self.model.contents.super_scores)

    def get_super_weights(self):
        """
        Get the CPCA super weights.

        Returns:
            List[List] : CPCA super weights.
        """
        return mx.matrix_to_list(self.model.contents.super_weights)

    def get_block_scores(self):
        """
        Get the CPCA block scores.

        Returns:
            List[List[List]] : CPCA block scores.
        """
        return tns.tensor_tolist(self.model.contents.block_scores)

    def get_block_loadings(self):
        """
        Get the CPCA block loadings.

        Returns:
            List[List[List]] : CPCA block loadings.
        """
        return tns.tensor_tolist(self.model.contents.block_loadings)

    def get_block_expvar(self):
        """
        Get the CPCA block explained variance.

        Returns:
            List[List[List]] : CPCA block explained variance.
        """
        return vlst.dvector_list_tolist(self.model.contents.block_expvar)

    def get_total_exp_variance(self):
        """
        Get the CPCA total explained variance.

        Returns:
            List[List[List]] : CPCA total explained variance.
        """
        return vect.dvector_tolist(self.model.contents.total_expvar)

    def predict(self, t_input):
        """
        Predict super scores and block scores using the CPCA model.

        Args:
            t_input : Tensor
                The input tensor.

        Returns:
            p_super_scores : List[List]
                CPCA predicted super scores.
            p_block_scores : List[List[List]]
                CPCA predicted block scores.
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

    def save(self, dbpath):
        """
        Save CPCA model to a sqlite3 file

        Parameters
        ----------
        dbpath (str): The path to the SQLite database file where the CPCA
            model will be stored.
        """
        write_cpca(dbpath, self.model)

    def load(self, dbpath):
        """
        Load CPCA model from a sqlite3 file

        Parameters
        ----------
        dbpath (str): The path to the SQLite database file where the CPCA
            model will be stored.
        """
        read_cpca(dbpath, self.model)
        self.npc = self.model.contents.super_scores[0].col
