# upls.py libscientific python binding
#
# Copyright (C) <2026>  Giuseppe Marco Randazzo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""Multi-way (unfolded) Partial Least Squares regression (UPLS).

UPLS relates three-way predictors to (two- or three-way) responses. Both
are given as lists of layers with the same objects::

    X = [x_layer_0, ..., x_layer_{Kx-1}]   # each n_objects x n_features
    Y = [y_layer_0, ..., y_layer_{Ky-1}]   # each n_objects x n_y

A single response matrix is a list with one layer: ``Y = [y]``. Each latent
variable has one X score per object and a weight/loading **matrix**
``n_features x Kx`` (Wold, Geladi and Esbensen, J. Chemometrics 1, 41-56,
1987).

    >>> from libscientific.upls import UPLS
    >>> X = [[[float(i + j + k) + 0.1 * ((i * 7 + j * 3 + k) % 5) for j in range(3)]
    ...       for i in range(12)] for k in range(2)]
    >>> Y = [[[sum(X[k][i][j] for k in range(2) for j in range(3))] for i in range(12)]]
    >>> model = UPLS(nlv=2)
    >>> model.fit(X, Y)
    >>> y_pred = model.predict(X)
    >>> len(y_pred), len(y_pred[0]), len(y_pred[0][0])   # 1 Y layer, 12 objects, 1 response
    (1, 12, 1)

The functions after the class are the low-level ctypes wrappers of ``upls.h``.
"""

from __future__ import annotations

import ctypes
from typing import Any

import libscientific.matrix as mx
import libscientific.tensor as tns
import libscientific.vectlist as vlst
import libscientific.vector as vect
from libscientific.loadlibrary import load_libscientific_library
from libscientific.upca import as_c_tensor

lsci = load_libscientific_library()


class UPLSMODEL(ctypes.Structure):
    """
    C structure holding a fitted UPLS model (low-level API).

    Use :class:`UPLS` instead of accessing this structure directly. With
    ``n`` objects, ``f`` features in ``Kx`` X layers, ``n_y`` responses in
    ``Ky`` Y layers and ``nlv`` latent variables:

    Attributes
    ----------
    xscores, yscores : POINTER(MATRIX)
        X and Y scores, ``n x nlv``.
    xloadings, xweights : POINTER(TENSOR)
        One ``f x Kx`` matrix per latent variable.
    yloadings : POINTER(TENSOR)
        One ``n_y x Ky`` matrix per latent variable.
    b : POINTER(DVECTOR)
        Inner regression coefficient of each latent variable.
    xvarexp : POINTER(DVECTOR)
        Percentage of the X variance explained by each latent variable.
    xcolaverage, xcolscaling, ycolaverage, ycolscaling : POINTER(DVECTLIST)
        Column means and scalings of each X and Y layer.
    recalculated_y, recalc_residuals : POINTER(TENSOR)
        Fitted Y of the training objects and residuals, one matrix per Y layer.

    The validation fields (``q2y``, ``sdep``, ``predicted_y``, ...) are filled
    by the cross-validation functions of the C library.
    """

    _fields_ = [
        ("xscores", ctypes.POINTER(mx.MATRIX)),
        ("xloadings", ctypes.POINTER(tns.TENSOR)),
        ("xweights", ctypes.POINTER(tns.TENSOR)),
        ("xvarexp", ctypes.POINTER(vect.DVECTOR)),
        ("yscores", ctypes.POINTER(mx.MATRIX)),
        ("yloadings", ctypes.POINTER(tns.TENSOR)),
        ("b", ctypes.POINTER(vect.DVECTOR)),
        ("xcolaverage", ctypes.POINTER(vlst.DVECTLIST)),
        ("xcolscaling", ctypes.POINTER(vlst.DVECTLIST)),
        ("ycolaverage", ctypes.POINTER(vlst.DVECTLIST)),
        ("ycolscaling", ctypes.POINTER(vlst.DVECTLIST)),
        ("r2x_model", ctypes.POINTER(vect.DVECTOR)),
        ("r2x_validation", ctypes.POINTER(vect.DVECTOR)),
        ("r2y_model", ctypes.POINTER(tns.TENSOR)),
        ("r2y_validation", ctypes.POINTER(tns.TENSOR)),
        ("q2y", ctypes.POINTER(tns.TENSOR)),
        ("sdep", ctypes.POINTER(tns.TENSOR)),
        ("sdec", ctypes.POINTER(tns.TENSOR)),
        ("recalculated_y", ctypes.POINTER(tns.TENSOR)),
        ("predicted_y", ctypes.POINTER(tns.TENSOR)),
        ("recalc_residuals", ctypes.POINTER(tns.TENSOR)),
        ("pred_residuals", ctypes.POINTER(tns.TENSOR)),
        ("q2y_yscrambling", ctypes.POINTER(tns.TENSOR)),
        ("sdep_yscrambling", ctypes.POINTER(tns.TENSOR)),
    ]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


class UPLS:
    """
    Multi-way Partial Least Squares regression model.

    Parameters
    ----------
    nlv : int, default 2
        Number of latent variables; capped to ``n_features``. After
        :meth:`fit`, ``nlv`` holds the number actually computed.
    xscaling, yscaling : int, default 1 and 0
        Column preprocessing of every X and Y layer: all options mean-center
        each column, then 0 = no scaling, 1 = unit variance (autoscaling),
        2 = root mean square, 3 = Pareto, 4 = range, 5 = level.

    Attributes
    ----------
    nlv, xscaling, yscaling : int
        See Parameters.
    model : ctypes.POINTER(UPLSMODEL)
        The underlying C model (low-level access, see :class:`UPLSMODEL`).

    Notes
    -----
    - The NIPALS iterations of UPLS stop at a relative change of ``1e-3``,
      a looser tolerance than PCA/PLS: results agree with other
      implementations to about 3 significant digits.
    - When ``Y`` is already fully explained, the remaining latent variables
      are skipped (a message is printed on standard error) and :attr:`nlv`
      is lowered accordingly.
    - For a single response, the inner coefficients ``b`` are 1: the
      regression coefficient is carried by the Y loadings.
    """

    def __init__(self, nlv: int = 2, xscaling: int = 1, yscaling: int = 0) -> None:
        self.model = new_upls_model()
        self.nlv = nlv
        self.xscaling = xscaling
        self.yscaling = yscaling

    def __repr__(self) -> str:
        return f"UPLS(nlv={self.nlv}, xscaling={self.xscaling}, yscaling={self.yscaling})"

    def __del__(self) -> None:
        if getattr(self, "model", None) is not None:
            del_upls_model(self.model)
            self.model = None  # type: ignore[assignment]  # already released

    def fit(self, x_input: Any, y_input: Any) -> None:
        """
        Fit the model on three-way training data.

        Calling ``fit`` again replaces the previous model.

        Parameters
        ----------
        x_input : list of matrices or 3D numpy array
            ``Kx`` X layers of shape ``(n_objects, n_features)``.
        y_input : list of matrices or 3D numpy array
            ``Ky`` Y layers of shape ``(n_objects, n_y)``; ``[y]`` for a
            single response matrix.

        Returns
        -------
        None
            The model is updated in place.
        """
        x_input_, x_owned = as_c_tensor(x_input)
        try:
            y_input_, y_owned = as_c_tensor(y_input)
            try:
                upls_algorithm(
                    x_input_, y_input_, self.nlv, self.xscaling, self.yscaling, self.model
                )
            finally:
                if y_owned:
                    tns.del_tensor(y_input_)
        finally:
            if x_owned:
                tns.del_tensor(x_input_)
        self.nlv = int(self.model.contents.xscores[0].col)

    def get_tscores(self) -> list[list[float]]:
        """
        Return the X scores of the training objects.

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_objects, nlv)``.
        """
        return mx.matrix_to_list(self.model.contents.xscores)

    def get_uscores(self) -> list[list[float]]:
        """
        Return the Y scores of the training objects.

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_objects, nlv)``.
        """
        return mx.matrix_to_list(self.model.contents.yscores)

    def get_xloadings(self) -> list[list[list[float]]]:
        """
        Return the X loading matrix of each latent variable.

        Returns
        -------
        list[list[list[float]]]
            ``nlv`` matrices of shape ``(n_features, Kx)``.
        """
        return tns.tensor_tolist(self.model.contents.xloadings)

    def get_xweights(self) -> list[list[list[float]]]:
        """
        Return the X weight matrix of each latent variable.

        Returns
        -------
        list[list[list[float]]]
            ``nlv`` matrices of shape ``(n_features, Kx)``.
        """
        return tns.tensor_tolist(self.model.contents.xweights)

    def get_yloadings(self) -> list[list[list[float]]]:
        """
        Return the Y loading matrix of each latent variable.

        Returns
        -------
        list[list[list[float]]]
            ``nlv`` matrices of shape ``(n_y, Ky)``.
        """
        return tns.tensor_tolist(self.model.contents.yloadings)

    def get_b(self) -> list[float]:
        """
        Return the inner regression coefficient of each latent variable.

        Returns
        -------
        list[float]
            ``nlv`` values.
        """
        return vect.dvector_tolist(self.model.contents.b)

    def get_exp_variance(self) -> list[float]:
        """
        Return the percentage of the X variance explained by each latent variable.

        Returns
        -------
        list[float]
            ``nlv`` values, not cumulative.
        """
        return vect.dvector_tolist(self.model.contents.xvarexp)

    def predict(self, x_input: Any, nlv_: int | None = None) -> list[list[list[float]]]:
        """
        Predict the responses of new objects.

        Parameters
        ----------
        x_input : list of matrices or 3D numpy array
            X layers of the new objects: same number of layers and features as
            the training data, any number of objects.
        nlv\\_ : int, optional
            Number of latent variables to use (default: all).

        Returns
        -------
        list[list[list[float]]]
            ``Ky`` predicted Y layers of shape ``(n_new_objects, n_y)``, in the
            original units of ``Y``.

        Raises
        ------
        ValueError
            If ``nlv_`` is not between 1 and :attr:`nlv`.
        """
        nlv = self.nlv if nlv_ is None else nlv_
        if not 1 <= nlv <= self.nlv:
            raise ValueError(f"nlv_ must be between 1 and {self.nlv}, got {nlv_}")
        x_input_, owned = as_c_tensor(x_input)
        tscores = mx.init_matrix()
        y_pred = tns.init_tensor()
        try:
            upls_score_predictor(x_input_, self.model, nlv, tscores)
            upls_y_predictor(tscores, self.model, nlv, y_pred)
            return tns.tensor_tolist(y_pred)
        finally:
            mx.del_matrix(tscores)
            tns.del_tensor(y_pred)
            if owned:
                tns.del_tensor(x_input_)

    def save(self, filepath: str) -> None:
        """
        Save the model to a file.

        Parameters
        ----------
        filepath : str
            Destination path. ``.sqlite3`` or ``.db`` use the legacy SQLite
            format; any other extension (``.lsm`` recommended) uses the binary
            libscientific model container.

        Raises
        ------
        OSError
            If the file cannot be written.
        """
        if filepath.endswith(".sqlite3") or filepath.endswith(".db"):
            write_upls(filepath, self.model)
        else:
            st = lsci.SaveUPLSModel(filepath.encode("utf-8"), self.model)
            if st != 0:
                raise OSError(f"SaveUPLSModel failed with error code {st}")

    def load(self, filepath: str) -> None:
        """
        Load a model saved with :meth:`save`, replacing the current one.

        After loading, :attr:`nlv` is set to the number of stored latent
        variables.

        Parameters
        ----------
        filepath : str
            Path of the model file; the format is detected automatically.

        Raises
        ------
        OSError
            If a binary model file cannot be read.
        """
        from libscientific.mlr import _is_binary_model

        if _is_binary_model(filepath):
            st = lsci.ImportUPLSModel(filepath.encode("utf-8"), self.model)
            if st != 0:
                raise OSError(f"ImportUPLSModel failed with error code {st}")
        else:
            read_upls(filepath, self.model)
        self.nlv = int(self.model.contents.xscores[0].col)


# ---------------------------------------------------------------------------
# Low-level API: thin ctypes wrappers around the C functions in upls.h.
# ---------------------------------------------------------------------------

lsci.NewUPLSModel.argtypes = [ctypes.POINTER(ctypes.POINTER(UPLSMODEL))]
lsci.NewUPLSModel.restype = None


def new_upls_model() -> ctypes._Pointer[UPLSMODEL]:
    """
    Allocate an empty C UPLS model (low-level API).

    Returns
    -------
    ctypes.POINTER(UPLSMODEL)
        A new, empty model; release it with :func:`del_upls_model`.
    """
    mupls = ctypes.POINTER(UPLSMODEL)()
    lsci.NewUPLSModel(ctypes.pointer(mupls))
    return mupls


lsci.DelUPLSModel.argtypes = [ctypes.POINTER(ctypes.POINTER(UPLSMODEL))]
lsci.DelUPLSModel.restype = None


def del_upls_model(mupls: ctypes._Pointer[UPLSMODEL]) -> None:
    """
    Free a C UPLS model created with :func:`new_upls_model` (low-level API).

    Parameters
    ----------
    mupls : ctypes.POINTER(UPLSMODEL)
        The model to release. It must not be used afterwards.
    """
    lsci.DelUPLSModel(ctypes.pointer(mupls))


lsci.UPLS.argtypes = [
    ctypes.POINTER(tns.TENSOR),
    ctypes.POINTER(tns.TENSOR),
    ctypes.c_size_t,
    ctypes.c_size_t,
    ctypes.c_size_t,
    ctypes.POINTER(UPLSMODEL),
    ctypes.POINTER(ctypes.c_int),
]
lsci.UPLS.restype = None


def upls_algorithm(
    x_input: ctypes._Pointer[tns.TENSOR],
    y_input: ctypes._Pointer[tns.TENSOR],
    nlv: int,
    xscaling: int,
    yscaling: int,
    mupls: ctypes._Pointer[UPLSMODEL],
) -> None:
    """
    Fit a UPLS model (low-level API, see :meth:`UPLS.fit`).

    Parameters
    ----------
    x_input, y_input : ctypes.POINTER(TENSOR)
        X and Y layers, all with the same number of rows.
    nlv : int
        Number of latent variables.
    xscaling, yscaling : int
        Column preprocessing of X and Y, 0-5.
    mupls : ctypes.POINTER(UPLSMODEL)
        Model created with :func:`new_upls_model`; filled with the results.
    """
    lsci.UPLS(x_input, y_input, nlv, xscaling, yscaling, mupls, None)


lsci.UPLSScorePredictor.argtypes = [
    ctypes.POINTER(tns.TENSOR),
    ctypes.POINTER(UPLSMODEL),
    ctypes.c_size_t,
    ctypes.POINTER(mx.MATRIX),
]
lsci.UPLSScorePredictor.restype = None


def upls_score_predictor(
    x_input: ctypes._Pointer[tns.TENSOR],
    mupls: ctypes._Pointer[UPLSMODEL],
    nlv: int,
    tscores: ctypes._Pointer[mx.MATRIX],
) -> None:
    """
    Compute the X scores of new objects (low-level API).

    Parameters
    ----------
    x_input : ctypes.POINTER(TENSOR)
        X layers of the new objects.
    mupls : ctypes.POINTER(UPLSMODEL)
        A fitted model.
    nlv : int
        Number of latent variables to use.
    tscores : ctypes.POINTER(MATRIX)
        Output matrix created with :func:`libscientific.matrix.init_matrix`;
        filled with the scores, ``m x nlv``.
    """
    lsci.UPLSScorePredictor(x_input, mupls, nlv, tscores)


lsci.UPLSYPredictor.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(UPLSMODEL),
    ctypes.c_size_t,
    ctypes.POINTER(tns.TENSOR),
]
lsci.UPLSYPredictor.restype = None


def upls_y_predictor(
    tscores: ctypes._Pointer[mx.MATRIX],
    mupls: ctypes._Pointer[UPLSMODEL],
    nlv: int,
    y_pred: ctypes._Pointer[tns.TENSOR],
) -> None:
    """
    Compute the responses from X scores (low-level API).

    Parameters
    ----------
    tscores : ctypes.POINTER(MATRIX)
        X scores from :func:`upls_score_predictor`.
    mupls : ctypes.POINTER(UPLSMODEL)
        A fitted model.
    nlv : int
        Number of latent variables to use.
    y_pred : ctypes.POINTER(TENSOR)
        Output tensor created with :func:`libscientific.tensor.init_tensor`;
        receives one ``m x n_y`` matrix per Y layer, in original units.
    """
    lsci.UPLSYPredictor(tscores, mupls, nlv, y_pred)


lsci.WriteUPLS.argtypes = [ctypes.c_char_p, ctypes.POINTER(UPLSMODEL)]
lsci.WriteUPLS.restype = None


def write_upls(dbpath: str, mupls: ctypes._Pointer[UPLSMODEL]) -> None:
    """
    Write a model to a SQLite database (legacy format, low-level API).

    Parameters
    ----------
    dbpath : str
        Path of the SQLite file to create.
    mupls : ctypes.POINTER(UPLSMODEL)
        The model to store.
    """
    lsci.WriteUPLS(dbpath.encode("utf-8"), mupls)


lsci.ReadUPLS.argtypes = [ctypes.c_char_p, ctypes.POINTER(UPLSMODEL)]
lsci.ReadUPLS.restype = None


def read_upls(dbpath: str, mupls: ctypes._Pointer[UPLSMODEL]) -> None:
    """
    Read a model from a SQLite database (legacy format, low-level API).

    Parameters
    ----------
    dbpath : str
        Path of the SQLite file written by :func:`write_upls`.
    mupls : ctypes.POINTER(UPLSMODEL)
        Model created with :func:`new_upls_model`; filled with the stored data.
    """
    lsci.ReadUPLS(dbpath.encode("utf-8"), mupls)


lsci.SaveUPLSModel.argtypes = [ctypes.c_char_p, ctypes.POINTER(UPLSMODEL)]
lsci.SaveUPLSModel.restype = ctypes.c_int
lsci.ImportUPLSModel.argtypes = [ctypes.c_char_p, ctypes.POINTER(UPLSMODEL)]
lsci.ImportUPLSModel.restype = ctypes.c_int
