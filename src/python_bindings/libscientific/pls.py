# pls.py libscientific python binding
#
# Copyright (C) <2023>  Giuseppe Marco Randazzo
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
"""Partial Least Squares regression (PLS).

PLS models the relationship between a matrix of predictors ``X``
(``n_objects x n_features``) and a matrix of responses ``Y``
(``n_objects x n_y``) through ``nlv`` latent variables::

    X = T @ P.T + E        Y = T @ diag(b) @ Q.T + F

``T`` are the X **scores** (``n_objects x nlv``), ``P`` the X **loadings**
(``n_features x nlv``), ``W`` the X **weights** used to compute the scores
of new objects, ``U`` the Y scores, ``Q`` the Y loadings (``n_y x nlv``)
and ``b`` the inner regression coefficients. Latent variables are computed
one at a time with the NIPALS algorithm.

Most users only need the :class:`PLS` class::

    >>> from libscientific.pls import PLS
    >>> X = [[1.0, 2.0, 0.5], [2.0, 1.0, 1.5], [3.0, 4.0, 1.0],
    ...      [4.0, 3.0, 2.5], [5.0, 5.0, 2.0], [6.0, 4.5, 3.5]]
    >>> y = [[3.1], [3.9], [7.2], [7.8], [10.1], [10.9]]
    >>> model = PLS(nlv=2, xscaling=1, yscaling=0)
    >>> model.fit(X, y)
    >>> y_pred, scores = model.predict([[3.5, 3.0, 2.0]], nlv_=2)
    >>> round(y_pred[0][0], 2)
    7.02

The functions after the class (``new_pls_model``, ``pls_algorithm``, ...)
are the low-level ctypes wrappers of the C API.
"""

from __future__ import annotations

import ctypes
from typing import Any

import libscientific.matrix as mx
import libscientific.tensor as t
import libscientific.vector as vect
from libscientific.loadlibrary import load_libscientific_library
from libscientific.matrix import MatrixLike

lsci = load_libscientific_library()


class PLSMODEL(ctypes.Structure):
    """
    C structure holding a fitted PLS model (low-level API).

    Use :class:`PLS` instead of accessing this structure directly. All fields
    are pointers to C matrices/vectors/tensors; convert them with
    :func:`libscientific.matrix.matrix_to_list` or
    :func:`libscientific.vector.dvector_tolist`.

    With ``n`` objects, ``f`` features, ``n_y`` responses and ``nlv`` latent
    variables, the fields filled by :func:`pls_algorithm` are:

    Attributes
    ----------
    xscores, xloadings, xweights : POINTER(MATRIX)
        ``T`` (``n x nlv``), ``P`` (``f x nlv``) and ``W`` (``f x nlv``).
    yscores, yloadings : POINTER(MATRIX)
        ``U`` (``n x nlv``) and ``Q`` (``n_y x nlv``).
    b : POINTER(DVECTOR)
        Inner regression coefficients, one per latent variable.
    xvarexp : POINTER(DVECTOR)
        Percentage of the X variance explained by each latent variable.
    xcolaverage, xcolscaling, ycolaverage, ycolscaling : POINTER(DVECTOR)
        Column means and scaling factors of X (size ``f``) and Y (size ``n_y``).
    recalculated_y, recalc_residuals : POINTER(MATRIX)
        Fitted Y of the training objects and their residuals,
        ``n x (n_y * nlv)``; column ``n_y * lv + j`` refers to response ``j``
        with ``lv + 1`` latent variables.

    The remaining fields (``predicted_y``, ``q2y``, ``sdep``, ``roc_*``,
    ``yscrambling``, ...) are filled by the validation functions of the C
    library (cross-validation, y-scrambling) and are empty after a plain fit.
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
        ("yscrambling", ctypes.POINTER(mx.MATRIX)),
    ]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


class PLS:
    """
    Partial Least Squares regression model (NIPALS algorithm).

    Parameters
    ----------
    nlv : int, default 2
        Number of latent variables to compute. It is capped to the number of
        features: after :meth:`fit`, ``nlv`` holds the number actually computed.
    xscaling : int, default 1
        Column preprocessing of ``X``. Every option mean-centers the columns,
        then: 0 = no scaling, 1 = unit variance (autoscaling), 2 = root mean
        square, 3 = Pareto, 4 = range, 5 = level (see
        :class:`~libscientific.pca.PCA` for details).
    yscaling : int, default 0
        Column preprocessing of ``Y``, same options as ``xscaling``.
        Predictions are always returned in the original units of ``Y``.

    Attributes
    ----------
    nlv : int
        Number of latent variables of the model.
    xscaling, yscaling : int
        The preprocessing options given at construction.
    model : ctypes.POINTER(PLSMODEL)
        The underlying C model (low-level access, see :class:`PLSMODEL`).

    Notes
    -----
    - ``Y`` is always 2D: one column per response. A single response is
      ``[[y1], [y2], ...]``.
    - The sign of each latent variable is arbitrary: scores, loadings and
      weights of a latent variable can change sign between platforms without
      changing the predictions.
    - Inputs can be lists of rows, numpy arrays or
      :class:`~libscientific.matrix.Matrix` objects (see
      :data:`~libscientific.matrix.MatrixLike`); results are Python lists.

    Examples
    --------
    >>> from libscientific.pls import PLS
    >>> X = [[1.0, 2.0, 0.5], [2.0, 1.0, 1.5], [3.0, 4.0, 1.0],
    ...      [4.0, 3.0, 2.5], [5.0, 5.0, 2.0], [6.0, 4.5, 3.5]]
    >>> y = [[3.1], [3.9], [7.2], [7.8], [10.1], [10.9]]
    >>> model = PLS(nlv=2, xscaling=1, yscaling=0)
    >>> model.fit(X, y)
    >>> [round(v, 2) for v in model.get_exp_variance()]   # % of X variance per LV
    [82.3, 17.47]

    Predict with 1 and 2 latent variables at once (one column per latent
    variable for a single response):

    >>> y_pred, t_scores = model.predict(X)
    >>> [round(v, 2) for v in y_pred[0]]
    [3.27, 3.42]

    Predict new objects with 2 latent variables:

    >>> y_new, _ = model.predict([[3.5, 3.0, 2.0]], nlv_=2)
    >>> round(y_new[0][0], 2)
    7.02
    """

    def __init__(self, nlv: int = 2, xscaling: int = 1, yscaling: int = 0) -> None:
        self.model = new_pls_model()
        self.nlv = nlv
        self.xscaling = xscaling
        self.yscaling = yscaling

    def __repr__(self) -> str:
        return f"PLS(nlv={self.nlv}, xscaling={self.xscaling}, yscaling={self.yscaling})"

    def __str__(self) -> str:
        return self.__repr__()

    def __del__(self) -> None:
        if getattr(self, "model", None) is not None:
            del_pls_model(self.model)
            self.model = None  # type: ignore[assignment]  # already released

    def _n_y(self) -> int:
        return int(self.model.contents.yloadings[0].row)

    def fit(self, x_input: MatrixLike, y_input: MatrixLike) -> None:
        """
        Fit the model on training data.

        ``X`` and ``Y`` are preprocessed according to :attr:`xscaling` and
        :attr:`yscaling`, then :attr:`nlv` latent variables are extracted.
        Calling ``fit`` again replaces the previous model.

        Parameters
        ----------
        x_input : MatrixLike
            Predictors, shape ``(n_objects, n_features)``.
        y_input : MatrixLike
            Responses, shape ``(n_objects, n_y)``: one column per response
            (``[[y1], [y2], ...]`` for a single response).

        Returns
        -------
        None
            The model is updated in place.

        Examples
        --------
        >>> model = PLS(nlv=5)
        >>> model.fit([[1.0, 2.0], [2.0, 1.0], [3.0, 4.0], [4.0, 3.0]],
        ...           [[1.0], [2.0], [3.0], [4.0]])
        >>> model.nlv        # capped to the number of features
        2
        """
        x_input_, x_owned = mx.as_c_matrix(x_input)
        try:
            y_input_, y_owned = mx.as_c_matrix(y_input)
            try:
                pls_algorithm(
                    x_input_,
                    y_input_,
                    self.model,
                    nlv=self.nlv,
                    x_scaling=self.xscaling,
                    y_scaling=self.yscaling,
                )
            finally:
                if y_owned:
                    mx.del_matrix(y_input_)
        finally:
            if x_owned:
                mx.del_matrix(x_input_)
        self.nlv = int(self.model.contents.xscores[0].col)

    def get_tscores(self) -> list[list[float]]:
        """
        Return the X scores ``T`` of the training objects.

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_objects, nlv)``: coordinates of each training
            object on each latent variable.

        See Also
        --------
        predict : also returns the X scores of new objects.

        Examples
        --------
        >>> model = PLS(nlv=2)
        >>> model.fit([[1.0, 2.0], [2.0, 1.0], [3.0, 4.0], [4.0, 3.0]],
        ...           [[1.0], [2.0], [3.0], [4.0]])
        >>> len(model.get_tscores()), len(model.get_tscores()[0])
        (4, 2)
        """
        return mx.matrix_to_list(self.model.contents.xscores)

    def get_uscores(self) -> list[list[float]]:
        """
        Return the Y scores ``U`` of the training objects.

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_objects, nlv)``. Plotting ``U[:, k]`` against
            ``T[:, k]`` shows the inner relation of latent variable ``k``.
        """
        return mx.matrix_to_list(self.model.contents.yscores)

    def get_ploadings(self) -> list[list[float]]:
        """
        Return the X loadings ``P``.

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_features, nlv)``: how each feature is
            reconstructed from each latent variable (``X ≈ T @ P.T``).
        """
        return mx.matrix_to_list(self.model.contents.xloadings)

    def get_qloadings(self) -> list[list[float]]:
        """
        Return the Y loadings ``Q``.

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_y, nlv)``. For a single response all the
            values are ``1.0``.
        """
        return mx.matrix_to_list(self.model.contents.yloadings)

    def get_weights(self) -> list[list[float]]:
        """
        Return the X weights ``W``.

        The weights define how the scores are computed from the
        (preprocessed) features; they are the usual way to interpret which
        features drive each latent variable.

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_features, nlv)``.
        """
        return mx.matrix_to_list(self.model.contents.xweights)

    def get_beta_coefficients(self, nlv: int = 1) -> list[float]:
        """
        Return the regression coefficients of the model with ``nlv`` latent variables.

        The coefficients ``B = W (P' W)^-1 b`` turn the PLS model into a
        linear equation on the **preprocessed** data::

            (y - y_average) / y_scaling = sum_j B[j] * (x[j] - x_average[j]) / x_scaling[j]

        so predictions can be computed without the model, e.g. with
        :func:`pls_beta_inference`.

        Parameters
        ----------
        nlv : int, default 1
            Number of latent variables, from 1 to :attr:`nlv`.

        Returns
        -------
        list[float]
            ``n_features`` coefficients, one per feature.

        Raises
        ------
        ValueError
            If the model has more than one response: the coefficients are
            computed for a single response only.

        Examples
        --------
        >>> X = [[1.0, 2.0, 0.5], [2.0, 1.0, 1.5], [3.0, 4.0, 1.0],
        ...      [4.0, 3.0, 2.5], [5.0, 5.0, 2.0], [6.0, 4.5, 3.5]]
        >>> model = PLS(nlv=2)
        >>> model.fit(X, [[3.1], [3.9], [7.2], [7.8], [10.1], [10.9]])
        >>> len(model.get_beta_coefficients(2))
        3
        """
        if self._n_y() > 1:
            raise ValueError(
                "beta coefficients are available only for models with a single response"
                f" (this model has {self._n_y()})"
            )
        beta_coeff = vect.init_dvector()
        try:
            pls_beta_coefficients(self.model, nlv, beta_coeff)
            return vect.dvector_tolist(beta_coeff)
        finally:
            vect.del_dvector(beta_coeff)

    def get_vip(self) -> list[float]:
        """
        Return the Variable Importance in Projection (VIP) of each feature.

        VIP summarizes how much each feature contributes to explaining ``Y``
        over all the latent variables of the model (S. Wold's definition).
        The average of the squared VIP values is 1, so features with
        ``VIP > 1`` are usually considered important and features with
        ``VIP < 0.8`` candidates for removal.

        Returns
        -------
        list[float]
            ``n_features`` values, one per feature.

        Examples
        --------
        >>> X = [[1.0, 2.0, 0.5], [2.0, 1.0, 1.5], [3.0, 4.0, 1.0],
        ...      [4.0, 3.0, 2.5], [5.0, 5.0, 2.0], [6.0, 4.5, 3.5]]
        >>> model = PLS(nlv=2)
        >>> model.fit(X, [[3.1], [3.9], [7.2], [7.8], [10.1], [10.9]])
        >>> vip = model.get_vip()
        >>> len(vip), round(sum(v * v for v in vip) / len(vip), 6)
        (3, 1.0)
        """
        vip = mx.init_matrix()
        try:
            pls_vip(self.model, vip)
            return [row[0] for row in mx.matrix_to_list(vip)]
        finally:
            mx.del_matrix(vip)

    def get_exp_variance(self) -> list[float]:
        """
        Return the percentage of the variance of ``X`` explained by each latent variable.

        Returns
        -------
        list[float]
            ``nlv`` values, not cumulative.
        """
        return vect.dvector_tolist(self.model.contents.xvarexp)

    def get_x_column_scaling(self) -> list[float]:
        """
        Return the scaling factor applied to each column of ``X``.

        Returns
        -------
        list[float]
            ``n_features`` values (all ``1.0`` when ``xscaling=0``).
        """
        return vect.dvector_tolist(self.model.contents.xcolscaling)

    def get_x_averages(self) -> list[float]:
        """
        Return the mean of each column of the training ``X``.

        Returns
        -------
        list[float]
            ``n_features`` values.
        """
        return vect.dvector_tolist(self.model.contents.xcolaverage)

    def get_y_column_scaling(self) -> list[float]:
        """
        Return the scaling factor applied to each column of ``Y``.

        Returns
        -------
        list[float]
            ``n_y`` values (all ``1.0`` when ``yscaling=0``).
        """
        return vect.dvector_tolist(self.model.contents.ycolscaling)

    def get_y_averages(self) -> list[float]:
        """
        Return the mean of each column of the training ``Y``.

        Returns
        -------
        list[float]
            ``n_y`` values.
        """
        return vect.dvector_tolist(self.model.contents.ycolaverage)

    def predict(
        self, x_input: MatrixLike, nlv_: int | None = None
    ) -> tuple[list[list[float]], list[list[float]]]:
        """
        Predict the responses of new objects.

        Parameters
        ----------
        x_input : MatrixLike
            Objects to predict, shape ``(n_new_objects, n_features)``, with the
            same features (same order) as the training data.
        nlv\\_ : int, optional
            Number of latent variables to use, from 1 to :attr:`nlv`. If
            omitted, the predictions of **every** number of latent variables
            are returned side by side (see Returns).

        Returns
        -------
        y_pred : list[list[float]]
            Predicted responses in the original units of ``Y``.

            - With ``nlv_``: shape ``(n_new_objects, n_y)``.
            - Without ``nlv_``: shape ``(n_new_objects, n_y * nlv)``, where
              column ``n_y * k + j`` is response ``j`` predicted with
              ``k + 1`` latent variables. For a single response, column ``k``
              is the prediction with ``k + 1`` latent variables.
        t_scores : list[list[float]]
            X scores of the new objects, shape ``(n_new_objects, nlv)``.

        Raises
        ------
        ValueError
            If ``nlv_`` is not between 1 and :attr:`nlv`.

        Examples
        --------
        Two responses, three latent variables:

        >>> X = [[1.0, 2.0, 0.5], [2.0, 1.0, 1.5], [3.0, 4.0, 1.0],
        ...      [4.0, 3.0, 2.5], [5.0, 5.0, 2.0], [6.0, 4.5, 3.5]]
        >>> Y = [[3.1, 0.2], [3.9, 0.9], [7.2, 0.8], [7.8, 1.9], [10.1, 1.7], [10.9, 2.8]]
        >>> model = PLS(nlv=3)
        >>> model.fit(X, Y)
        >>> y_all, t = model.predict(X)
        >>> len(y_all[0]), len(t[0])      # 2 responses x 3 LVs, 3 LVs
        (6, 3)
        >>> y_2lv, _ = model.predict(X, nlv_=2)
        >>> len(y_2lv[0])                 # both responses, 2 LVs
        2
        >>> y_2lv[0] == y_all[0][2:4]
        True
        """
        if nlv_ is not None and not 1 <= nlv_ <= self.nlv:
            raise ValueError(f"nlv_ must be between 1 and {self.nlv}, got {nlv_}")
        x_input_, owned = mx.as_c_matrix(x_input)
        p_scores_ = mx.init_matrix()
        p_y_ = mx.init_matrix()
        try:
            pls_y_predictor_all_lv(x_input_, self.model, p_scores_, p_y_)
            p_scores = mx.matrix_to_list(p_scores_)
            p_y = mx.matrix_to_list(p_y_)
        finally:
            mx.del_matrix(p_scores_)
            mx.del_matrix(p_y_)
            if owned:
                mx.del_matrix(x_input_)
        if nlv_ is not None:
            n_y = self._n_y()
            p_y = [row[n_y * (nlv_ - 1) : n_y * nlv_] for row in p_y]
        return p_y, p_scores

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
            write_pls(filepath, self.model)
        else:
            from libscientific.model_io import SavePLSModel

            st = SavePLSModel(filepath.encode("utf-8"), self.model)
            if st != 0:
                raise OSError(f"SavePLSModel failed with error code {st}")

    def load(self, filepath: str) -> None:
        """
        Load a model saved with :meth:`save`, replacing the current one.

        The format is detected from the file content or the extension. After
        loading, :attr:`nlv` is set to the number of latent variables stored
        in the file; :attr:`xscaling` and :attr:`yscaling` are left unchanged
        (predictions use the stored column averages and scalings).

        Parameters
        ----------
        filepath : str
            Path of the model file.

        Raises
        ------
        OSError
            If a binary model file cannot be read.

        Examples
        --------
        >>> import os, tempfile
        >>> X = [[1.0, 2.0], [2.0, 1.0], [3.0, 4.0], [4.0, 3.0]]
        >>> model = PLS(nlv=2)
        >>> model.fit(X, [[1.0], [2.0], [3.0], [4.0]])
        >>> path = os.path.join(tempfile.mkdtemp(), "pls.lsm")
        >>> model.save(path)
        >>> loaded = PLS()
        >>> loaded.load(path)
        >>> loaded.predict(X) == model.predict(X)
        True
        """
        is_lscm = False
        try:
            with open(filepath, "rb") as f:
                if f.read(4) == b"LSCM":
                    is_lscm = True
        except OSError:
            pass

        if is_lscm or filepath.endswith(".lsm") or filepath.endswith(".lsb"):
            from libscientific.model_io import ImportPLSModel

            st = ImportPLSModel(filepath.encode("utf-8"), self.model)
            if st != 0:
                raise OSError(f"ImportPLSModel failed with error code {st}")
        else:
            read_pls(filepath, self.model)
        self.nlv = int(self.model.contents.xscores[0].col)


def apply_scaling(
    x_input: list[list[float]], x_avg: list[float], x_scal: list[float]
) -> list[list[float]]:
    """
    Center and scale a matrix with given column averages and scalings.

    Computes ``(x[i][j] - x_avg[j]) / x_scal[j]`` for every element.

    Parameters
    ----------
    x_input : list[list[float]]
        Matrix, shape ``(n_objects, n_features)``.
    x_avg, x_scal : list[float]
        Column averages and scalings (size ``n_features``), e.g. from
        :meth:`PLS.get_x_averages` and :meth:`PLS.get_x_column_scaling`.

    Returns
    -------
    list[list[float]]
        The preprocessed matrix, same shape as ``x_input``.

    Examples
    --------
    >>> apply_scaling([[3.0, 10.0]], [1.0, 4.0], [2.0, 3.0])
    [[1.0, 2.0]]
    """
    x_scaled: list[list[float]] = []
    for x_row in x_input:
        x_scaled.append([(x_val - x_avg[j]) / x_scal[j] for j, x_val in enumerate(x_row)])
    return x_scaled


def revert_scaling(
    x_input: list[Any], x_avg: list[float], x_scal: list[float]
) -> list[list[float]]:
    """
    Undo :func:`apply_scaling`: ``x[i][j] * x_scal[j] + x_avg[j]``.

    Parameters
    ----------
    x_input : list[list[float]] or list[float]
        Preprocessed matrix, or a flat list of values of a single column.
    x_avg, x_scal : list[float]
        Column averages and scalings.

    Returns
    -------
    list[list[float]]
        The matrix in the original units. A flat input list gives one row
        per value (``[[v1], [v2], ...]``).

    Examples
    --------
    >>> revert_scaling([[1.0, 2.0]], [1.0, 4.0], [2.0, 3.0])
    [[3.0, 10.0]]
    >>> revert_scaling([0.5, 1.0], [10.0], [2.0])
    [[11.0], [12.0]]
    """
    x_revert: list[list[float]] = []
    for x_row in x_input:
        if isinstance(x_row, list):
            x_revert.append([x_val * x_scal[j] + x_avg[j] for j, x_val in enumerate(x_row)])
        else:
            x_revert.append([x_row * x_scal[0] + x_avg[0]])
    return x_revert


def pls_beta_inference(x_input: list[list[float]], **kwargs: list[float]) -> list[list[float]]:
    """
    Predict a single response from the beta coefficients only.

    Useful to deploy a model as a simple linear equation, without the
    library or the full model. The result is the same as
    :meth:`PLS.predict` with the same number of latent variables.

    Parameters
    ----------
    x_input : list[list[float]]
        Objects to predict, shape ``(n_objects, n_features)``.
    **kwargs
        All required:

        - ``x_averages``, ``x_scaling``: from :meth:`PLS.get_x_averages` and
          :meth:`PLS.get_x_column_scaling`.
        - ``y_averages``, ``y_scaling``: from :meth:`PLS.get_y_averages` and
          :meth:`PLS.get_y_column_scaling`.
        - ``beta_coeff``: from :meth:`PLS.get_beta_coefficients`.

    Returns
    -------
    list[list[float]]
        Predictions, shape ``(n_objects, 1)``, in the original units of ``Y``.

    Examples
    --------
    >>> X = [[1.0, 2.0, 0.5], [2.0, 1.0, 1.5], [3.0, 4.0, 1.0],
    ...      [4.0, 3.0, 2.5], [5.0, 5.0, 2.0], [6.0, 4.5, 3.5]]
    >>> model = PLS(nlv=2, xscaling=1, yscaling=0)
    >>> model.fit(X, [[3.1], [3.9], [7.2], [7.8], [10.1], [10.9]])
    >>> y_beta = pls_beta_inference([[3.5, 3.0, 2.0]],
    ...                             x_averages=model.get_x_averages(),
    ...                             x_scaling=model.get_x_column_scaling(),
    ...                             y_averages=model.get_y_averages(),
    ...                             y_scaling=model.get_y_column_scaling(),
    ...                             beta_coeff=model.get_beta_coefficients(2))
    >>> y_model, _ = model.predict([[3.5, 3.0, 2.0]], nlv_=2)
    >>> round(y_beta[0][0], 8) == round(y_model[0][0], 8)
    True
    """
    x_avg = kwargs["x_averages"]
    x_scal = kwargs["x_scaling"]
    y_avg = kwargs["y_averages"]
    y_scal = kwargs["y_scaling"]
    b_coeff = kwargs["beta_coeff"]
    x_ps = mx.new_matrix(apply_scaling(x_input, x_avg, x_scal))
    pby_ = vect.new_dvector([0 for _ in range(len(x_input))])
    b_coeff_ = vect.new_dvector(b_coeff)
    mx.matrix_dvector_dot_product(x_ps, b_coeff_, pby_)
    pby = vect.dvector_tolist(pby_)
    vect.del_dvector(pby_)
    vect.del_dvector(b_coeff_)
    mx.del_matrix(x_ps)
    return revert_scaling(pby, y_avg, y_scal)


# ---------------------------------------------------------------------------
# Low-level API: thin ctypes wrappers around the C functions in pls.h.
# They take and fill raw C structures; prefer the PLS class above.
# ---------------------------------------------------------------------------

lsci.NewPLSModel.argtypes = [ctypes.POINTER(ctypes.POINTER(PLSMODEL))]
lsci.NewPLSModel.restype = None


def new_pls_model() -> ctypes._Pointer[PLSMODEL]:
    """
    Allocate an empty C PLS model (low-level API).

    Returns
    -------
    ctypes.POINTER(PLSMODEL)
        A new, empty model. It must be released with :func:`del_pls_model`.
    """
    mpls = ctypes.POINTER(PLSMODEL)()
    lsci.NewPLSModel(ctypes.pointer(mpls))
    return mpls


lsci.DelPLSModel.argtypes = [ctypes.POINTER(ctypes.POINTER(PLSMODEL))]
lsci.DelPLSModel.restype = None


def del_pls_model(mpls: ctypes._Pointer[PLSMODEL]) -> None:
    """
    Free a C PLS model created with :func:`new_pls_model` (low-level API).

    Parameters
    ----------
    mpls : ctypes.POINTER(PLSMODEL)
        The model to release. It must not be used afterwards.
    """
    lsci.DelPLSModel(ctypes.pointer(mpls))


lsci.PLS.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(mx.MATRIX),
    ctypes.c_size_t,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(PLSMODEL),
    ctypes.POINTER(ctypes.c_int),
]
lsci.PLS.restype = None


def pls_algorithm(
    x_input: ctypes._Pointer[mx.MATRIX],
    y_input: ctypes._Pointer[mx.MATRIX],
    mpls: ctypes._Pointer[PLSMODEL],
    **kwargs: int,
) -> None:
    """
    Fit a PLS model with the NIPALS algorithm (low-level API, see :meth:`PLS.fit`).

    Parameters
    ----------
    x_input : ctypes.POINTER(MATRIX)
        Predictors, shape ``(n_objects, n_features)``.
    y_input : ctypes.POINTER(MATRIX)
        Responses, shape ``(n_objects, n_y)``.
    mpls : ctypes.POINTER(PLSMODEL)
        Model created with :func:`new_pls_model`; filled with the results.
    **kwargs
        Required: ``nlv`` (number of latent variables), ``x_scaling`` and
        ``y_scaling`` (preprocessing options 0-5, see :class:`PLS`).
    """
    nlv = kwargs["nlv"]
    x_scaling = kwargs["x_scaling"]
    y_scaling = kwargs["y_scaling"]
    scisignal = ctypes.c_int(0)
    lsci.PLS(x_input, y_input, nlv, x_scaling, y_scaling, mpls, ctypes.pointer(scisignal))


lsci.PLSBetasCoeff.argtypes = [
    ctypes.POINTER(PLSMODEL),
    ctypes.c_size_t,
    ctypes.POINTER(vect.DVECTOR),
]
lsci.PLSBetasCoeff.restype = None


def pls_beta_coefficients(
    mpls: ctypes._Pointer[PLSMODEL], nlv: int, bcoeff: ctypes._Pointer[vect.DVECTOR]
) -> None:
    """
    Compute the beta coefficients of a single-response model (low-level API).

    See :meth:`PLS.get_beta_coefficients`; the result is wrong for models
    with more than one response.

    Parameters
    ----------
    mpls : ctypes.POINTER(PLSMODEL)
        A fitted model.
    nlv : int
        Number of latent variables to use.
    bcoeff : ctypes.POINTER(DVECTOR)
        Output vector created with :func:`libscientific.vector.init_dvector`;
        filled with ``n_features`` coefficients.
    """
    lsci.PLSBetasCoeff(mpls, nlv, bcoeff)


lsci.PLSVIP.argtypes = [ctypes.POINTER(PLSMODEL), ctypes.POINTER(mx.MATRIX)]
lsci.PLSVIP.restype = None


def pls_vip(mpls: ctypes._Pointer[PLSMODEL], vip: ctypes._Pointer[mx.MATRIX]) -> None:
    """
    Compute the VIP of each feature (low-level API, see :meth:`PLS.get_vip`).

    Parameters
    ----------
    mpls : ctypes.POINTER(PLSMODEL)
        A fitted model.
    vip : ctypes.POINTER(MATRIX)
        Output matrix created with :func:`libscientific.matrix.init_matrix`;
        resized to ``(n_features, 1)``.
    """
    lsci.PLSVIP(mpls, vip)


lsci.PLSScorePredictor.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(PLSMODEL),
    ctypes.c_size_t,
    ctypes.POINTER(mx.MATRIX),
]
lsci.PLSScorePredictor.restype = None


def pls_score_predictor(
    x_input: ctypes._Pointer[mx.MATRIX],
    mpls: ctypes._Pointer[PLSMODEL],
    nlv: int,
    p_scores: ctypes._Pointer[mx.MATRIX],
) -> None:
    """
    Compute the X scores of new objects (low-level API).

    Parameters
    ----------
    x_input : ctypes.POINTER(MATRIX)
        Objects, shape ``(n_new_objects, n_features)``.
    mpls : ctypes.POINTER(PLSMODEL)
        A fitted model.
    nlv : int
        Number of latent variables to use.
    p_scores : ctypes.POINTER(MATRIX)
        Output matrix created with :func:`libscientific.matrix.init_matrix`;
        filled with the scores, shape ``(n_new_objects, nlv)``.
    """
    lsci.PLSScorePredictor(x_input, mpls, nlv, p_scores)


lsci.PLSYPredictor.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(PLSMODEL),
    ctypes.c_size_t,
    ctypes.POINTER(mx.MATRIX),
]
lsci.PLSYPredictor.restype = None


def pls_y_predictor(
    x_scores: Any, mpls: ctypes._Pointer[PLSMODEL], nlv: int, predicted_y: Any
) -> None:
    """
    Compute the responses from X scores (low-level API).

    Parameters
    ----------
    x_scores : ctypes.POINTER(MATRIX) or Matrix
        X scores, e.g. from :func:`pls_score_predictor`, shape ``(n, >= nlv)``.
    mpls : ctypes.POINTER(PLSMODEL)
        A fitted model.
    nlv : int
        Number of latent variables to use.
    predicted_y : ctypes.POINTER(MATRIX) or Matrix
        Output matrix; filled with the responses in original units, shape
        ``(n, n_y)``.
    """
    x_scores_ = x_scores.mtx if isinstance(x_scores, mx.Matrix) else x_scores
    predicted_y_ = predicted_y.mtx if isinstance(predicted_y, mx.Matrix) else predicted_y
    lsci.PLSYPredictor(x_scores_, mpls, nlv, predicted_y_)


lsci.PLSYPredictorAllLV.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(PLSMODEL),
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(mx.MATRIX),
]
lsci.PLSYPredictorAllLV.restype = None


def pls_y_predictor_all_lv(
    x_input: ctypes._Pointer[mx.MATRIX],
    mpls: ctypes._Pointer[PLSMODEL],
    predicted_scores: ctypes._Pointer[mx.MATRIX],
    predicted_y: ctypes._Pointer[mx.MATRIX],
) -> None:
    """
    Predict the responses for every number of latent variables (low-level API).

    See :meth:`PLS.predict`.

    Parameters
    ----------
    x_input : ctypes.POINTER(MATRIX)
        Objects, shape ``(n, n_features)``.
    mpls : ctypes.POINTER(PLSMODEL)
        A fitted model.
    predicted_scores : ctypes.POINTER(MATRIX)
        Output: X scores, shape ``(n, nlv)``.
    predicted_y : ctypes.POINTER(MATRIX)
        Output: responses, shape ``(n, n_y * nlv)``; column ``n_y * k + j`` is
        response ``j`` with ``k + 1`` latent variables.
    """
    lsci.PLSYPredictorAllLV(x_input, mpls, predicted_scores, predicted_y)


lsci.PrintPLSModel.argtypes = [ctypes.POINTER(PLSMODEL)]
lsci.PrintPLSModel.restype = None


def print_pls(mpls: ctypes._Pointer[PLSMODEL]) -> None:
    """
    Print the content of a model to standard output (low-level API).

    Parameters
    ----------
    mpls : ctypes.POINTER(PLSMODEL)
        The model to print (e.g. ``PLS().model``).
    """
    lsci.PrintPLSModel(mpls)


lsci.WritePLS.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.POINTER(PLSMODEL)]
lsci.WritePLS.restype = None


def write_pls(dbpath: str, pls: ctypes._Pointer[PLSMODEL]) -> None:
    """
    Write a model to a SQLite database (legacy format, low-level API).

    Prefer :meth:`PLS.save`.

    Parameters
    ----------
    dbpath : str
        Path of the SQLite file to create.
    pls : ctypes.POINTER(PLSMODEL)
        The model to store.
    """
    encoded_string = dbpath.encode("utf-8")
    ctypes_string = ctypes.c_char_p(encoded_string)
    lsci.WritePLS(ctypes_string, pls)


lsci.ReadPLS.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.POINTER(PLSMODEL)]
lsci.ReadPLS.restype = None


def read_pls(dbpath: str, pls: ctypes._Pointer[PLSMODEL]) -> None:
    """
    Read a model from a SQLite database (legacy format, low-level API).

    Prefer :meth:`PLS.load`.

    Parameters
    ----------
    dbpath : str
        Path of the SQLite file written by :func:`write_pls`.
    pls : ctypes.POINTER(PLSMODEL)
        Model created with :func:`new_pls_model`; filled with the stored data.
    """
    encoded_string = dbpath.encode("utf-8")
    ctypes_string = ctypes.c_char_p(encoded_string)
    lsci.ReadPLS(ctypes_string, pls)
