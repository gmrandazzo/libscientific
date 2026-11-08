# pca.py libscientific python binding
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
"""Principal Component Analysis (PCA).

PCA decomposes a data matrix ``X`` (``n_objects`` rows x ``n_features``
columns) into ``npc`` principal components::

    X_preprocessed = T @ P.T + E

where ``T`` are the **scores** (``n_objects x npc``, the coordinates of each
object in the new space), ``P`` are the **loadings** (``n_features x npc``,
the weight of each feature in each component) and ``E`` is the residual.
Components are computed one at a time with the NIPALS algorithm, so only the
requested components are calculated.

Most users only need the :class:`PCA` class::

    >>> from libscientific.pca import PCA
    >>> X = [[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [3.0, 4.0, 2.0],
    ...      [4.0, 3.0, 1.0], [5.0, 5.0, 0.5]]
    >>> model = PCA(scaling=1, npc=2)
    >>> model.fit(X)
    >>> [round(v, 2) for v in model.get_exp_variance()]
    [90.4, 6.67]

The functions below the class (``new_pca_model``, ``pca_algorithm``, ...) are
the low-level ctypes wrappers of the C API; they work with raw C pointers and
are only needed to interoperate with other low-level code.
"""

from __future__ import annotations

import ctypes

import libscientific.matrix as mx
import libscientific.vector as vect
from libscientific.loadlibrary import load_libscientific_library
from libscientific.matrix import MatrixLike

lsci = load_libscientific_library()


class PCAMODEL(ctypes.Structure):
    """
    C structure holding a fitted PCA model (low-level API).

    Use :class:`PCA` instead of accessing this structure directly. All fields
    are pointers to C matrices/vectors; convert them with
    :func:`libscientific.matrix.matrix_to_list` or
    :func:`libscientific.vector.dvector_tolist`.

    Attributes
    ----------
    scores : POINTER(MATRIX)
        Scores ``T``, shape ``(n_objects, npc)``.
    loadings : POINTER(MATRIX)
        Loadings ``P``, shape ``(n_features, npc)``.
    dmodx : POINTER(MATRIX)
        Distance of each object to the model (DModX), shape ``(n_objects, npc)``:
        column ``k`` is the distance after removing components ``0..k``.
    varexp : POINTER(DVECTOR)
        Percentage of variance explained by each component, size ``npc``.
    colaverage : POINTER(DVECTOR)
        Column means used for centering, size ``n_features``.
    colscaling : POINTER(DVECTOR)
        Column scaling factors, size ``n_features`` (all ``1.0`` when ``scaling=0``).
    """

    _fields_ = [
        ("scores", ctypes.POINTER(mx.MATRIX)),
        ("loadings", ctypes.POINTER(mx.MATRIX)),
        ("dmodx", ctypes.POINTER(mx.MATRIX)),
        ("varexp", ctypes.POINTER(vect.DVECTOR)),
        ("colaverage", ctypes.POINTER(vect.DVECTOR)),
        ("colscaling", ctypes.POINTER(vect.DVECTOR)),
    ]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


class PCA:
    """
    Principal Component Analysis model (NIPALS algorithm).

    Parameters
    ----------
    scaling : int, default 1
        Column preprocessing applied before the decomposition. Every option
        mean-centers the columns; the options differ in how they are scaled:

        ====  ==========================================================
        0     mean centering only (no scaling)
        1     unit variance / autoscaling (divide by the standard deviation)
        2     root mean square of the column
        3     Pareto (divide by the square root of the standard deviation)
        4     range (divide by max - min)
        5     level (divide by the column mean)
        ====  ==========================================================

        Use ``1`` when the features have different units or ranges,
        ``0`` when they are already on the same scale.
    npc : int, default 2
        Number of principal components to compute. It is capped to the number
        of features: after :meth:`fit`, ``npc`` holds the number of components
        actually computed.

    Attributes
    ----------
    scaling : int
        The preprocessing option given at construction.
    npc : int
        Number of principal components of the model.
    model : ctypes.POINTER(PCAMODEL)
        The underlying C model (low-level access, see :class:`PCAMODEL`).

    Notes
    -----
    - The sign of each component is arbitrary: scores and loadings of a
      component can be multiplied by -1 between runs or platforms without
      changing the model.
    - Missing values can be passed as ``float('nan')`` (or as the value
      returned by :func:`libscientific.misc.missing_value`): they are
      skipped during the decomposition.
    - Inputs can be lists of rows, numpy arrays or :class:`~libscientific.matrix.Matrix`
      objects (see :data:`~libscientific.matrix.MatrixLike`). Results are
      returned as Python lists; wrap them with ``numpy.array(...)`` if needed.

    Examples
    --------
    >>> from libscientific.pca import PCA
    >>> X = [[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [3.0, 4.0, 2.0],
    ...      [4.0, 3.0, 1.0], [5.0, 5.0, 0.5]]
    >>> model = PCA(scaling=1, npc=2)
    >>> model.fit(X)
    >>> scores = model.get_scores()        # 5 objects x 2 components
    >>> len(scores), len(scores[0])
    (5, 2)
    >>> loadings = model.get_loadings()    # 3 features x 2 components
    >>> len(loadings), len(loadings[0])
    (3, 2)
    >>> [round(v, 2) for v in model.get_exp_variance()]   # % per component
    [90.4, 6.67]
    >>> new_scores = model.predict([[2.5, 2.5, 2.5]])     # project new objects
    >>> len(new_scores), len(new_scores[0])
    (1, 2)

    Save the model and load it back:

    >>> import os, tempfile
    >>> path = os.path.join(tempfile.mkdtemp(), "model.lsm")
    >>> model.save(path)
    >>> model2 = PCA()
    >>> model2.load(path)
    >>> model2.npc
    2
    """

    def __init__(self, scaling: int = 1, npc: int = 2) -> None:
        self.model = new_pca_model()
        self.scaling = scaling
        self.npc = npc

    def __repr__(self) -> str:
        return f"PCA(scaling={self.scaling}, npc={self.npc})"

    def __del__(self) -> None:
        if getattr(self, "model", None) is not None:
            del_pca_model(self.model)
            self.model = None  # type: ignore[assignment]  # already released

    def fit(self, m_input: MatrixLike) -> None:
        """
        Compute the principal components of a matrix.

        The columns are preprocessed according to :attr:`scaling`, then
        :attr:`npc` components are extracted. Results are stored in the model
        and read with :meth:`get_scores`, :meth:`get_loadings`,
        :meth:`get_exp_variance` and :meth:`get_dmodx`. Calling ``fit`` again
        replaces the previous model.

        Parameters
        ----------
        m_input : MatrixLike
            Training data, shape ``(n_objects, n_features)``: one row per
            object, one column per feature.

        Returns
        -------
        None
            The model is updated in place.

        Examples
        --------
        >>> model = PCA(scaling=1, npc=5)
        >>> model.fit([[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [3.0, 4.0, 2.0]])
        >>> model.npc       # capped to the number of features
        3
        """
        m_input_, owned = mx.as_c_matrix(m_input)
        try:
            pca_algorithm(m_input_, self.scaling, self.npc, self.model)
        finally:
            if owned:
                mx.del_matrix(m_input_)
        self.npc = int(self.model.contents.scores[0].col)

    def get_scores(self) -> list[list[float]]:
        """
        Return the scores of the training objects.

        The scores are the coordinates of each training object in the space of
        the principal components.

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_objects, npc)``: row ``i`` holds the
            coordinates of training object ``i``, column ``k`` refers to
            component ``k``.

        See Also
        --------
        predict : scores of new objects.

        Examples
        --------
        >>> model = PCA(npc=2)
        >>> model.fit([[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [3.0, 4.0, 2.0], [4.0, 3.0, 1.0]])
        >>> scores = model.get_scores()
        >>> len(scores), len(scores[0])
        (4, 2)
        """
        return mx.matrix_to_list(self.model.contents.scores)

    def get_loadings(self) -> list[list[float]]:
        """
        Return the loadings of the features.

        The loadings describe how much each original feature contributes to
        each principal component. Each column is a unit-length vector.

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_features, npc)``: row ``j`` refers to feature
            ``j``, column ``k`` to component ``k``.

        Examples
        --------
        >>> model = PCA(npc=2)
        >>> model.fit([[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [3.0, 4.0, 2.0], [4.0, 3.0, 1.0]])
        >>> loadings = model.get_loadings()
        >>> len(loadings), len(loadings[0])
        (3, 2)
        >>> round(sum(row[0] ** 2 for row in loadings), 6)   # unit-length columns
        1.0
        """
        return mx.matrix_to_list(self.model.contents.loadings)

    def get_exp_variance(self) -> list[float]:
        """
        Return the variance explained by each component, in percent.

        Values are per component (not cumulative) and relative to the total
        variance of the preprocessed training matrix, so their sum is at most
        100.

        Returns
        -------
        list[float]
            ``npc`` values; element ``k`` is the percentage of variance
            explained by component ``k``.

        Examples
        --------
        >>> model = PCA(scaling=1, npc=3)
        >>> model.fit([[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [3.0, 4.0, 2.0],
        ...            [4.0, 3.0, 1.0], [5.0, 5.0, 0.5]])
        >>> [round(v, 2) for v in model.get_exp_variance()]
        [90.4, 6.67, 2.94]
        >>> round(sum(model.get_exp_variance()), 6)   # all components: 100%
        100.0
        """
        return vect.dvector_tolist(self.model.contents.varexp)

    def get_dmodx(self) -> list[list[float]]:
        """
        Return the distance of each training object to the model (DModX).

        DModX is the norm of the residual of an object, i.e. the part of the
        object that is not explained by the components. Large values point to
        objects that the model does not describe well (possible outliers).

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_objects, npc)``: element ``[i][k]`` is the
            distance of object ``i`` after removing components ``0..k``, so
            each row decreases from left to right.

        Examples
        --------
        >>> model = PCA(npc=2)
        >>> model.fit([[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [3.0, 4.0, 2.0], [4.0, 3.0, 1.0]])
        >>> dmodx = model.get_dmodx()
        >>> len(dmodx), len(dmodx[0])
        (4, 2)
        >>> all(row[0] >= row[1] for row in dmodx)
        True
        """
        return mx.matrix_to_list(self.model.contents.dmodx)

    def predict(self, m_input: MatrixLike) -> list[list[float]]:
        """
        Project new objects into the model and return their scores.

        The new data are preprocessed with the column means and scalings of
        the training data, then projected onto the loadings. Projecting the
        training matrix itself gives the same values as :meth:`get_scores`.

        Parameters
        ----------
        m_input : MatrixLike
            Objects to project, shape ``(n_new_objects, n_features)``. The
            columns must be the same features, in the same order, as the
            training data.

        Returns
        -------
        list[list[float]]
            Predicted scores, shape ``(n_new_objects, npc)``.

        See Also
        --------
        reconstruct_original_matrix : go back from scores to the original features.

        Examples
        --------
        >>> X = [[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [3.0, 4.0, 2.0], [4.0, 3.0, 1.0]]
        >>> model = PCA(npc=2)
        >>> model.fit(X)
        >>> new_scores = model.predict([[2.5, 2.5, 2.5], [1.0, 1.0, 1.0]])
        >>> len(new_scores), len(new_scores[0])
        (2, 2)
        """
        m_input_, owned = mx.as_c_matrix(m_input)
        pscores_ = mx.init_matrix()
        try:
            pca_score_predictor(m_input_, self.model, self.npc, pscores_)
            return mx.matrix_to_list(pscores_)
        finally:
            mx.del_matrix(pscores_)
            if owned:
                mx.del_matrix(m_input_)

    def reconstruct_original_matrix(
        self, npc_input: int | None = None, scores_input: MatrixLike | None = None
    ) -> list[list[float]]:
        """
        Rebuild data in the original feature space from scores.

        Computes ``T @ P.T`` with the first ``npc_input`` components, then
        reverts the scaling and centering, so the result is in the same units
        as the training data. With all the components the training matrix is
        reconstructed exactly; with fewer components the result is a denoised
        approximation.

        Parameters
        ----------
        npc_input : int, optional
            Number of components to use, from 1 to :attr:`npc`.
            Default: all the components of the model.
        scores_input : MatrixLike, optional
            Scores to transform back, shape ``(n_objects, npc)``, e.g. the
            output of :meth:`predict`. Default: the scores of the training data.

        Returns
        -------
        list[list[float]]
            Reconstructed data, shape ``(n_objects, n_features)``.

        Examples
        --------
        >>> X = [[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [3.0, 4.0, 2.0], [4.0, 3.0, 1.0]]
        >>> model = PCA(npc=3)
        >>> model.fit(X)
        >>> Xr = model.reconstruct_original_matrix()            # all components
        >>> [round(v, 6) for v in Xr[0]]
        [1.0, 2.0, 3.0]
        >>> Xr2 = model.reconstruct_original_matrix(npc_input=1)  # approximation
        >>> len(Xr2), len(Xr2[0])
        (4, 3)

        Reconstruct new objects from their predicted scores:

        >>> s = model.predict([[2.5, 2.5, 2.5]])
        >>> [round(v, 6) for v in model.reconstruct_original_matrix(scores_input=s)[0]]
        [2.5, 2.5, 2.5]
        """
        npc_ = self.npc if npc_input is None else npc_input
        if scores_input is None:
            scores_, owned = self.model.contents.scores, False
        else:
            scores_, owned = mx.as_c_matrix(scores_input)
        ind_vars = mx.init_matrix()
        try:
            pca_ind_var_predictor(scores_, self.model, npc_, ind_vars)
            return mx.matrix_to_list(ind_vars)
        finally:
            mx.del_matrix(ind_vars)
            if owned:
                mx.del_matrix(scores_)

    def save(self, filepath: str) -> None:
        """
        Save the model to a file.

        Parameters
        ----------
        filepath : str
            Destination path. The format is chosen from the extension:
            ``.sqlite3`` or ``.db`` use the legacy SQLite format; any other
            extension (``.lsm`` recommended) uses the binary libscientific
            model container, which is faster and needs no extra dependency.

        Raises
        ------
        OSError
            If the file cannot be written.

        Examples
        --------
        >>> import os, tempfile
        >>> model = PCA(npc=2)
        >>> model.fit([[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [3.0, 4.0, 2.0], [4.0, 3.0, 1.0]])
        >>> model.save(os.path.join(tempfile.mkdtemp(), "pca.lsm"))
        """
        if filepath.endswith(".sqlite3") or filepath.endswith(".db"):
            write_pca(filepath, self.model)
        else:
            from libscientific.model_io import SavePCAModel

            st = SavePCAModel(filepath.encode("utf-8"), self.model)
            if st != 0:
                raise OSError(f"SavePCAModel failed with error code {st}")

    def load(self, filepath: str) -> None:
        """
        Load a model saved with :meth:`save`, replacing the current one.

        The format is detected from the file content (binary container) or
        the extension (``.lsm``/``.lsb``; anything else is read as SQLite).
        After loading, :attr:`npc` is set to the number of components stored
        in the file. :attr:`scaling` is not stored and is left unchanged; it
        is not needed for :meth:`predict`, which uses the saved column
        averages and scalings.

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
        >>> path = os.path.join(tempfile.mkdtemp(), "pca.lsm")
        >>> X = [[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [3.0, 4.0, 2.0], [4.0, 3.0, 1.0]]
        >>> model = PCA(npc=2)
        >>> model.fit(X)
        >>> model.save(path)
        >>> loaded = PCA()
        >>> loaded.load(path)
        >>> loaded.npc
        2
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
            from libscientific.model_io import ImportPCAModel

            st = ImportPCAModel(filepath.encode("utf-8"), self.model)
            if st != 0:
                raise OSError(f"ImportPCAModel failed with error code {st}")
        else:
            read_pca(filepath, self.model)
        self.npc = int(self.model.contents.scores[0].col)


# ---------------------------------------------------------------------------
# Low-level API: thin ctypes wrappers around the C functions in pca.h.
# They take and fill raw C structures; prefer the PCA class above.
# ---------------------------------------------------------------------------

lsci.NewPCAModel.argtypes = [ctypes.POINTER(ctypes.POINTER(PCAMODEL))]
lsci.NewPCAModel.restype = None


def new_pca_model() -> ctypes._Pointer[PCAMODEL]:
    """
    Allocate an empty C PCA model (low-level API).

    Returns
    -------
    ctypes.POINTER(PCAMODEL)
        A new, empty model. It must be released with :func:`del_pca_model`.

    Examples
    --------
    >>> m = new_pca_model()
    >>> del_pca_model(m)
    """
    mpca = ctypes.POINTER(PCAMODEL)()
    lsci.NewPCAModel(ctypes.pointer(mpca))
    return mpca


lsci.DelPCAModel.argtypes = [ctypes.POINTER(ctypes.POINTER(PCAMODEL))]
lsci.DelPCAModel.restype = None


def del_pca_model(mpca: ctypes._Pointer[PCAMODEL]) -> None:
    """
    Free a C PCA model created with :func:`new_pca_model` (low-level API).

    Parameters
    ----------
    mpca : ctypes.POINTER(PCAMODEL)
        The model to release. It must not be used afterwards.
    """
    lsci.DelPCAModel(ctypes.pointer(mpca))


lsci.PCA.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.c_int,
    ctypes.c_size_t,
    ctypes.POINTER(PCAMODEL),
    ctypes.POINTER(ctypes.c_int),
]
lsci.PCA.restype = None


def pca_algorithm(
    m_input: ctypes._Pointer[mx.MATRIX], scaling: int, npc: int, mpca: ctypes._Pointer[PCAMODEL]
) -> None:
    """
    Fit a PCA model with the NIPALS algorithm (low-level API, see :meth:`PCA.fit`).

    Parameters
    ----------
    m_input : ctypes.POINTER(MATRIX)
        Training data, shape ``(n_objects, n_features)``.
    scaling : int
        Column preprocessing, 0-5 (see :class:`PCA`).
    npc : int
        Number of components to compute (capped to ``n_features``).
    mpca : ctypes.POINTER(PCAMODEL)
        Model created with :func:`new_pca_model`; it is filled with the results.

    Examples
    --------
    >>> import libscientific.matrix as mx
    >>> X = mx.new_matrix([[1.0, 2.0], [2.0, 1.0], [3.0, 4.0]])
    >>> m = new_pca_model()
    >>> pca_algorithm(X, 1, 2, m)
    >>> m[0].scores[0].row, m[0].scores[0].col
    (3, 2)
    >>> del_pca_model(m); mx.del_matrix(X)
    """
    scisignal = ctypes.c_int(0)
    lsci.PCA(m_input, scaling, npc, mpca, scisignal)


lsci.PCAScorePredictor.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(PCAMODEL),
    ctypes.c_size_t,
    ctypes.POINTER(mx.MATRIX),
]
lsci.PCAScorePredictor.restype = None


def pca_score_predictor(
    m_input: ctypes._Pointer[mx.MATRIX],
    mpca: ctypes._Pointer[PCAMODEL],
    npc: int,
    pscores: ctypes._Pointer[mx.MATRIX],
) -> None:
    """
    Project new objects into a fitted model (low-level API, see :meth:`PCA.predict`).

    Parameters
    ----------
    m_input : ctypes.POINTER(MATRIX)
        Objects to project, shape ``(n_new_objects, n_features)``.
    mpca : ctypes.POINTER(PCAMODEL)
        A fitted model.
    npc : int
        Number of components to use (capped to the components of the model).
    pscores : ctypes.POINTER(MATRIX)
        Output matrix created with :func:`libscientific.matrix.init_matrix`;
        it is resized and filled with the scores, shape ``(n_new_objects, npc)``.
    """
    lsci.PCAScorePredictor(m_input, mpca, npc, pscores)


lsci.PCAIndVarPredictor.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(vect.DVECTOR),
    ctypes.POINTER(vect.DVECTOR),
    ctypes.c_size_t,
    ctypes.POINTER(mx.MATRIX),
]
lsci.PCAIndVarPredictor.restype = None


def pca_ind_var_predictor(
    scores_input: ctypes._Pointer[mx.MATRIX],
    pca_model: ctypes._Pointer[PCAMODEL],
    npc: int,
    ind_vars_out: ctypes._Pointer[mx.MATRIX],
) -> None:
    """
    Reconstruct original-space data from scores (low-level API).

    See :meth:`PCA.reconstruct_original_matrix`.

    Parameters
    ----------
    scores_input : ctypes.POINTER(MATRIX)
        Scores, shape ``(n_objects, >= npc)``.
    pca_model : ctypes.POINTER(PCAMODEL)
        A fitted model (its loadings, column averages and scalings are used).
    npc : int
        Number of components to use.
    ind_vars_out : ctypes.POINTER(MATRIX)
        Output matrix created with :func:`libscientific.matrix.init_matrix`;
        filled with the reconstructed data, shape ``(n_objects, n_features)``.
    """
    lsci.PCAIndVarPredictor(
        scores_input,
        pca_model.contents.loadings,
        pca_model.contents.colaverage,
        pca_model.contents.colscaling,
        npc,
        ind_vars_out,
    )


lsci.PrintPCA.argtypes = [ctypes.POINTER(PCAMODEL)]
lsci.PrintPCA.restype = None


def print_pca(mpca: ctypes._Pointer[PCAMODEL]) -> None:
    """
    Print scores, loadings and explained variance of a model to stdout (low-level API).

    Parameters
    ----------
    mpca : ctypes.POINTER(PCAMODEL)
        The model to print (e.g. ``PCA().model``).
    """
    lsci.PrintPCA(mpca)


lsci.WritePCA.argtypes = [ctypes.c_char_p, ctypes.POINTER(PCAMODEL)]
lsci.WritePCA.restype = None


def write_pca(dbpath: str, pca: ctypes._Pointer[PCAMODEL]) -> None:
    """
    Write a model to a SQLite database (legacy format, low-level API).

    Prefer :meth:`PCA.save`, which also supports the binary ``.lsm`` format.

    Parameters
    ----------
    dbpath : str
        Path of the SQLite file to create.
    pca : ctypes.POINTER(PCAMODEL)
        The model to store.
    """
    encoded_string = dbpath.encode("utf-8")
    ctypes_string = ctypes.c_char_p(encoded_string)
    lsci.WritePCA(ctypes_string, pca)


lsci.ReadPCA.argtypes = [ctypes.c_char_p, ctypes.POINTER(PCAMODEL)]
lsci.ReadPCA.restype = None


def read_pca(dbpath: str, pca: ctypes._Pointer[PCAMODEL]) -> None:
    """
    Read a model from a SQLite database (legacy format, low-level API).

    Prefer :meth:`PCA.load`, which detects the file format automatically.

    Parameters
    ----------
    dbpath : str
        Path of the SQLite file written by :func:`write_pca`.
    pca : ctypes.POINTER(PCAMODEL)
        Model created with :func:`new_pca_model`; it is filled with the stored data.
    """
    encoded_string = dbpath.encode("utf-8")
    ctypes_string = ctypes.c_char_p(encoded_string)
    lsci.ReadPCA(ctypes_string, pca)
