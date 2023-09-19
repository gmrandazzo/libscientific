"""pca.py libscientific python binding

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
                          pca_model,
                          npc,
                          ind_vars_out):
    """
    pca_ind_var_predictor: Reconstruct the original matrix
                           from PCA model using scores and loadings
    """
    lsci.PCAIndVarPredictor(scores_input,
                            pca_model.contents.loadings,
                            pca_model.contents.colaverage,
                            pca_model.contents.colscaling,
                            npc,
                            ind_vars_out)

lsci.PrintPCA.argtypes = [ctypes.POINTER(PCAMODEL)]
lsci.PrintPCA.restype = None


def print_pca(mpca):
    """
    PrintPCA: Print to video the PCA Model
    """
    lsci.PrintPCA(mpca)

lsci.WritePCA.argtypes = [
    ctypes.c_char_p,
    ctypes.POINTER(PCAMODEL)
]
lsci.WritePCA.restype = None


def write_pca(dbpath: str, pca: PCAMODEL):
    """
    write_pca(dbpath, pca):
    
    Writes a PCA (Principal Component Analysis) model to a SQLite database.

    Parameters:
        dbpath (str): The path to the SQLite database file where the PCA 
            model will be stored.
        pca (PCAMODEL): The libscientific data structure representing the PCA 
            model.

    Returns:
        None
    """
    encoded_string = dbpath.encode('utf-8')
    ctypes_string = ctypes.c_char_p(encoded_string)
    lsci.WritePCA(ctypes_string, pca)


lsci.ReadPCA.argtypes = [
    ctypes.c_char_p,
    ctypes.POINTER(PCAMODEL)
]
lsci.ReadPCA.restype = None

def read_pca(dbpath: str, pca: PCAMODEL):
    """
    read_pca(dbpath, pca):
    
    Reads a PCA (Principal Component Analysis) model from a SQLite database.

    Parameters:
        dbpath (str): The path to the SQLite database file from which the PCA
            model will be read.
        pca (PCAMODEL): The libscientific data structure where the PCA
            model will be loaded.

    Returns:
        None
    """
    encoded_string = dbpath.encode('utf-8')
    ctypes_string = ctypes.c_char_p(encoded_string)
    lsci.ReadPCA(ctypes_string, pca)

class PCA():
    """
    Principal Component Analysis (PCA) Model

    This class provides methods for creating and utilizing a PCA model for dimensionality reduction
    and feature extraction.

    Attributes:
        mpca (CDataType): The PCA model data.
        scaling (int): Scaling option for PCA.
        npc (int): Number of principal components.

    Methods:
        __init__(self, scaling, npc)
        __del__(self)
        fit(self, m_input)
        get_scores(self)
        get_loadings(self)
        get_exp_variance(self)
        predict(self, m_input)
        reconstruct_original_matrix(self, npc_input=None, scores_input=None)
    """

    def __init__(self, scaling=1, npc=2):
        """
        Initialize a PCA instance.

        Args:
            scaling (int): Scaling option for PCA.
            npc (int): Number of principal components.
        """
        self.model = new_pca_model()
        self.scaling = scaling
        self.npc = npc

    def __del__(self):
        """
        Clean up resources associated with the PCA instance.
        """
        if self.model is not None:
            del_pca_model(self.model)
            del self.model
        self.model = None

    def fit(self, m_input):
        """
        Fit the PCA model using an input matrix.

        Args:
            m_input (Matrix or List[List]): Input matrix for fitting the PCA model.
        """
        if "Matrix" in str(type(m_input)):
            pca_algorithm(m_input.mtx, self.scaling, self.npc, self.model)
        else:
            m_input_ = mx.new_matrix(m_input)
            pca_algorithm(m_input_, self.scaling, self.npc, self.model)
            mx.del_matrix(m_input_)
            del m_input_

    def get_scores(self):
        """
        Get the PCA scores.

        Returns:
            List[List[float]]: The PCA scores.
        """
        return mx.matrix_to_list(self.model.contents.scores)

    def get_loadings(self):
        """
        Get the PCA loadings.

        Returns:
            List[List[float]]: The PCA loadings.
        """
        return mx.matrix_to_list(self.model.contents.loadings)

    def get_exp_variance(self):
        """
        Get the explained variance of each principal component.

        Returns:
            List[float]: The explained variance of each principal component.
        """
        return vect.dvector_tolist(self.model.contents.varexp)

    def predict(self, m_input):
        """
        Project an input matrix into the PCA model.

        Args:
            m_input (Matrix or List[List]): Input matrix for projecting into the PCA model.

        Returns:
            List[List[float]]: Projected scores using the PCA model.
        """
        pscores_ = mx.init_matrix()
        if "Matrix" in str(type(m_input)):
            m_input_ = m_input
            pca_score_predictor(m_input.mtx, self.model, self.npc, pscores_)
        else:
            m_input_ = mx.new_matrix(m_input)
            pca_score_predictor(m_input_, self.model, self.npc, pscores_)
            mx.del_matrix(m_input_)
            del m_input_
        pscores = mx.matrix_to_list(pscores_)
        mx.del_matrix(pscores_)
        del pscores_
        return pscores

    def reconstruct_original_matrix(self, npc_input=None, scores_input=None):
        """
        Reconstruct the original input matrix using a specified number of principal components.

        Args:
            npc_input (int): Number of principal components to be used. Default is None.
            scores_input (CDataType): Scores used for reconstruction. Default is None.

        Returns:
            List[List[float]]: Reconstructed matrix using the specified principal components.
        """
        npc_ = None
        if npc_input is None:
            npc_ = self.npc
        else:
            npc_ = npc_input

        scores_ = None
        if scores_input is None:
            scores_ = self.model.contents.scores
        else:
            scores_ = scores_input

        ind_vars = mx.init_matrix()
        pca_ind_var_predictor(scores_,
                              self.model,
                              npc_,
                              ind_vars)
        omx = mx.matrix_to_list(ind_vars)
        mx.del_matrix(ind_vars)
        del ind_vars
        return omx

    def save(self, dbpath):
        """
        Save PCA model to a sqlite3 file

        Parameters
        ----------
        dbpath (str): The path to the SQLite database file where the PCA 
            model will be stored.
        """
        write_pca(dbpath, self.model)

    def load(self, dbpath):
        """
        Load PCA model from a sqlite3 file

        Parameters
        ----------
        dbpath (str): The path to the SQLite database file where the PCA 
            model will be stored.
        """
        read_pca(dbpath, self.model)
        self.npc = self.model.contents.scores[0].col
