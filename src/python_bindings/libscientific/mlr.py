# mlr.py libscientific python binding
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
"""Multiple Linear Regression (MLR) by ordinary least squares.

MLR fits, for every response ``y_k``, the linear equation::

    y_k = b_0k + b_1k * x_1 + b_2k * x_2 + ... + b_fk * x_f

where ``b_0k`` is the intercept. Each response is fitted independently.

    >>> from libscientific.mlr import MLR
    >>> X = [[1.0, 0.0], [2.0, 1.0], [3.0, 4.0], [4.0, 2.0], [5.0, 3.0]]
    >>> y = [[1.0 + 2.0 * a - 0.5 * b] for a, b in X]
    >>> model = MLR()
    >>> model.fit(X, y)
    >>> [round(row[0], 6) for row in model.get_coefficients()]   # intercept, x1, x2
    [1.0, 2.0, -0.5]
    >>> round(model.predict([[10.0, 2.0]])[0][0], 6)
    20.0

The functions after the class are the low-level ctypes wrappers of ``mlr.h``.
"""

from __future__ import annotations

import ctypes

import libscientific.matrix as mx
import libscientific.vector as vect
from libscientific.loadlibrary import load_libscientific_library
from libscientific.matrix import MatrixLike

lsci = load_libscientific_library()


class MLRMODEL(ctypes.Structure):
    """
    C structure holding a fitted MLR model (low-level API).

    Use :class:`MLR` instead of accessing this structure directly. With
    ``n`` objects, ``f`` features and ``n_y`` responses:

    Attributes
    ----------
    b : POINTER(MATRIX)
        Coefficients, ``(f + 1) x n_y``: row 0 is the intercept, row ``j + 1``
        the coefficient of feature ``j``.
    recalculated_y, recalc_residuals : POINTER(MATRIX)
        Fitted responses of the training objects and residuals
        (fitted - true), ``n x n_y``.
    ymean : POINTER(DVECTOR)
        Mean of each training response.
    r2y_model, sdec : POINTER(DVECTOR)
        R² and standard deviation of the errors of the training fit, one
        value per response.
    predicted_y, pred_residuals, q2y, sdep, bias, r2q2scrambling
        Filled by the validation functions of the C library
        (cross-validation, y-scrambling); empty after a plain fit.
    """

    _fields_ = [
        ("b", ctypes.POINTER(mx.MATRIX)),
        ("recalculated_y", ctypes.POINTER(mx.MATRIX)),
        ("predicted_y", ctypes.POINTER(mx.MATRIX)),
        ("recalc_residuals", ctypes.POINTER(mx.MATRIX)),
        ("pred_residuals", ctypes.POINTER(mx.MATRIX)),
        ("ymean", ctypes.POINTER(vect.DVECTOR)),
        ("r2y_model", ctypes.POINTER(vect.DVECTOR)),
        ("q2y", ctypes.POINTER(vect.DVECTOR)),
        ("sdep", ctypes.POINTER(vect.DVECTOR)),
        ("sdec", ctypes.POINTER(vect.DVECTOR)),
        ("bias", ctypes.POINTER(vect.DVECTOR)),
        ("r2q2scrambling", ctypes.POINTER(mx.MATRIX)),
    ]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


class MLR:
    """
    Multiple Linear Regression model (ordinary least squares).

    Attributes
    ----------
    model : ctypes.POINTER(MLRMODEL)
        The underlying C model (low-level access, see :class:`MLRMODEL`).

    Notes
    -----
    - ``Y`` is always 2D: one column per response (``[[y1], [y2], ...]`` for
      a single response). Every response gets its own coefficients.
    - The data are not scaled: coefficients are in the units of the data.
    - OLS needs more objects than features, and features that are not
      linear combinations of each other; otherwise prefer
      :class:`~libscientific.pls.PLS`.

    Examples
    --------
    Two responses at once:

    >>> X = [[1.0, 0.0], [2.0, 1.0], [3.0, 4.0], [4.0, 2.0], [5.0, 3.0]]
    >>> Y = [[1.0 + 2.0 * a - 0.5 * b, 3.0 - a + b] for a, b in X]
    >>> model = MLR()
    >>> model.fit(X, Y)
    >>> [[round(v, 6) for v in row] for row in model.get_coefficients()]
    [[1.0, 3.0], [2.0, -1.0], [-0.5, 1.0]]
    >>> [round(v, 6) for v in model.get_r2()]       # perfect fit
    [1.0, 1.0]
    """

    def __init__(self) -> None:
        self.model = new_mlr_model()

    def __repr__(self) -> str:
        return "MLR()"

    def __del__(self) -> None:
        if getattr(self, "model", None) is not None:
            del_mlr_model(self.model)
            self.model = None  # type: ignore[assignment]  # already released

    def fit(self, x_input: MatrixLike, y_input: MatrixLike) -> None:
        """
        Fit the coefficients by ordinary least squares.

        Calling ``fit`` again replaces the previous model.

        Parameters
        ----------
        x_input : MatrixLike
            Predictors, shape ``(n_objects, n_features)``.
        y_input : MatrixLike
            Responses, shape ``(n_objects, n_y)``.

        Returns
        -------
        None
            The model is updated in place.
        """
        x_input_, x_owned = mx.as_c_matrix(x_input)
        try:
            y_input_, y_owned = mx.as_c_matrix(y_input)
            try:
                mlr_algorithm(x_input_, y_input_, self.model)
            finally:
                if y_owned:
                    mx.del_matrix(y_input_)
        finally:
            if x_owned:
                mx.del_matrix(x_input_)

    def get_coefficients(self) -> list[list[float]]:
        """
        Return the regression coefficients.

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_features + 1, n_y)``: row 0 holds the
            intercepts, row ``j + 1`` the coefficients of feature ``j``;
            column ``k`` refers to response ``k``.
        """
        return mx.matrix_to_list(self.model.contents.b)

    def get_recalculated_y(self) -> list[list[float]]:
        """
        Return the fitted responses of the training objects.

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_objects, n_y)``.
        """
        return mx.matrix_to_list(self.model.contents.recalculated_y)

    def get_residuals(self) -> list[list[float]]:
        """
        Return the training residuals (fitted minus true response).

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_objects, n_y)``.
        """
        return mx.matrix_to_list(self.model.contents.recalc_residuals)

    def get_r2(self) -> list[float]:
        """
        Return the coefficient of determination R² of the training fit.

        Returns
        -------
        list[float]
            One value per response: ``1 - SS(residuals) / SS(y - mean(y))``.
        """
        return vect.dvector_tolist(self.model.contents.r2y_model)

    def get_sdec(self) -> list[float]:
        """
        Return the standard deviation of the training errors (SDEC).

        Returns
        -------
        list[float]
            One value per response: ``sqrt(SS(residuals) / n_objects)``, in
            the units of the response.
        """
        return vect.dvector_tolist(self.model.contents.sdec)

    def predict(self, x_input: MatrixLike) -> list[list[float]]:
        """
        Predict the responses of new objects.

        Parameters
        ----------
        x_input : MatrixLike
            Objects, shape ``(n_new_objects, n_features)``, with the same
            features (same order) as the training data.

        Returns
        -------
        list[list[float]]
            Predicted responses, shape ``(n_new_objects, n_y)``.

        Raises
        ------
        ValueError
            If the number of features differs from the training data.

        Examples
        --------
        >>> model = MLR()
        >>> model.fit([[0.0], [1.0], [2.0]], [[1.0], [3.0], [5.0]])   # y = 1 + 2x
        >>> [[round(v, 6) for v in row] for row in model.predict([[10.0], [-1.0]])]
        [[21.0], [-1.0]]
        """
        x_input_, owned = mx.as_c_matrix(x_input)
        n_coef = int(self.model.contents.b[0].row)
        try:
            if int(x_input_[0].col) != n_coef - 1:
                raise ValueError(
                    f"the model has {n_coef - 1} features, the input has {x_input_[0].col}"
                )
            predicted = mx.init_matrix()
            try:
                mlr_predict_y(x_input_, self.model, predicted)
                return mx.matrix_to_list(predicted)
            finally:
                mx.del_matrix(predicted)
        finally:
            if owned:
                mx.del_matrix(x_input_)

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
            write_mlr(filepath, self.model)
        else:
            st = lsci.SaveMLRModel(filepath.encode("utf-8"), self.model)
            if st != 0:
                raise OSError(f"SaveMLRModel failed with error code {st}")

    def load(self, filepath: str) -> None:
        """
        Load a model saved with :meth:`save`, replacing the current one.

        Parameters
        ----------
        filepath : str
            Path of the model file; the format is detected automatically.

        Raises
        ------
        OSError
            If a binary model file cannot be read.

        Examples
        --------
        >>> import os, tempfile
        >>> model = MLR()
        >>> model.fit([[0.0], [1.0], [2.0]], [[1.0], [3.0], [5.0]])
        >>> path = os.path.join(tempfile.mkdtemp(), "mlr.lsm")
        >>> model.save(path)
        >>> loaded = MLR()
        >>> loaded.load(path)
        >>> loaded.predict([[4.0]]) == model.predict([[4.0]])
        True
        """
        if _is_binary_model(filepath):
            st = lsci.ImportMLRModel(filepath.encode("utf-8"), self.model)
            if st != 0:
                raise OSError(f"ImportMLRModel failed with error code {st}")
        else:
            read_mlr(filepath, self.model)


def _is_binary_model(filepath: str) -> bool:
    """True if the file is a binary libscientific model container."""
    if filepath.endswith(".lsm") or filepath.endswith(".lsb"):
        return True
    try:
        with open(filepath, "rb") as f:
            return f.read(4) == b"LSCM"
    except OSError:
        return False


# ---------------------------------------------------------------------------
# Low-level API: thin ctypes wrappers around the C functions in mlr.h.
# ---------------------------------------------------------------------------

lsci.NewMLRModel.argtypes = [ctypes.POINTER(ctypes.POINTER(MLRMODEL))]
lsci.NewMLRModel.restype = None


def new_mlr_model() -> ctypes._Pointer[MLRMODEL]:
    """
    Allocate an empty C MLR model (low-level API).

    Returns
    -------
    ctypes.POINTER(MLRMODEL)
        A new, empty model; release it with :func:`del_mlr_model`.
    """
    mmlr = ctypes.POINTER(MLRMODEL)()
    lsci.NewMLRModel(ctypes.pointer(mmlr))
    return mmlr


lsci.DelMLRModel.argtypes = [ctypes.POINTER(ctypes.POINTER(MLRMODEL))]
lsci.DelMLRModel.restype = None


def del_mlr_model(mmlr: ctypes._Pointer[MLRMODEL]) -> None:
    """
    Free a C MLR model created with :func:`new_mlr_model` (low-level API).

    Parameters
    ----------
    mmlr : ctypes.POINTER(MLRMODEL)
        The model to release. It must not be used afterwards.
    """
    lsci.DelMLRModel(ctypes.pointer(mmlr))


lsci.MLR.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(MLRMODEL),
    ctypes.POINTER(ctypes.c_int),
]
lsci.MLR.restype = None


def mlr_algorithm(
    x_input: ctypes._Pointer[mx.MATRIX],
    y_input: ctypes._Pointer[mx.MATRIX],
    mmlr: ctypes._Pointer[MLRMODEL],
) -> None:
    """
    Fit an MLR model (low-level API, see :meth:`MLR.fit`).

    Parameters
    ----------
    x_input : ctypes.POINTER(MATRIX)
        Predictors, ``n x f``.
    y_input : ctypes.POINTER(MATRIX)
        Responses, ``n x n_y``.
    mmlr : ctypes.POINTER(MLRMODEL)
        Model created with :func:`new_mlr_model`; filled with the results.
    """
    lsci.MLR(x_input, y_input, mmlr, None)


lsci.MLRPredictY.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(MLRMODEL),
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(vect.DVECTOR),
    ctypes.POINTER(vect.DVECTOR),
]
lsci.MLRPredictY.restype = None


def mlr_predict_y(
    x_input: ctypes._Pointer[mx.MATRIX],
    mmlr: ctypes._Pointer[MLRMODEL],
    predicted_y: ctypes._Pointer[mx.MATRIX],
) -> None:
    """
    Predict the responses of new objects (low-level API, see :meth:`MLR.predict`).

    Parameters
    ----------
    x_input : ctypes.POINTER(MATRIX)
        Objects, ``m x f``.
    mmlr : ctypes.POINTER(MLRMODEL)
        A fitted model.
    predicted_y : ctypes.POINTER(MATRIX)
        Output matrix created with :func:`libscientific.matrix.init_matrix`;
        resized to ``m x n_y``.
    """
    lsci.MLRPredictY(x_input, None, mmlr, predicted_y, None, None, None)


lsci.WriteMLR.argtypes = [ctypes.c_char_p, ctypes.POINTER(MLRMODEL)]
lsci.WriteMLR.restype = None


def write_mlr(dbpath: str, mmlr: ctypes._Pointer[MLRMODEL]) -> None:
    """
    Write a model to a SQLite database (legacy format, low-level API).

    Parameters
    ----------
    dbpath : str
        Path of the SQLite file to create.
    mmlr : ctypes.POINTER(MLRMODEL)
        The model to store.
    """
    lsci.WriteMLR(dbpath.encode("utf-8"), mmlr)


lsci.ReadMLR.argtypes = [ctypes.c_char_p, ctypes.POINTER(MLRMODEL)]
lsci.ReadMLR.restype = None


def read_mlr(dbpath: str, mmlr: ctypes._Pointer[MLRMODEL]) -> None:
    """
    Read a model from a SQLite database (legacy format, low-level API).

    Parameters
    ----------
    dbpath : str
        Path of the SQLite file written by :func:`write_mlr`.
    mmlr : ctypes.POINTER(MLRMODEL)
        Model created with :func:`new_mlr_model`; filled with the stored data.
    """
    lsci.ReadMLR(dbpath.encode("utf-8"), mmlr)


lsci.SaveMLRModel.argtypes = [ctypes.c_char_p, ctypes.POINTER(MLRMODEL)]
lsci.SaveMLRModel.restype = ctypes.c_int
lsci.ImportMLRModel.argtypes = [ctypes.c_char_p, ctypes.POINTER(MLRMODEL)]
lsci.ImportMLRModel.restype = ctypes.c_int
