# lda.py libscientific python binding
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
"""Linear Discriminant Analysis (Fisher LDA) for classification.

LDA finds the directions that best separate the classes (maximum ratio of
between-class to within-class scatter) and classifies an object into the
class with the highest linear discriminant score::

    score_k(x) = mu_k' S^-1 x - 1/2 mu_k' S^-1 mu_k + ln(prior_k)

where ``mu_k`` is the mean of class ``k``, ``S`` the pooled within-class
covariance and ``prior_k`` the fraction of training objects in class ``k``.

Class labels are integers numbered from 0 (``0..K-1``) or from 1 (``1..K``),
without gaps; predictions use the same numbering.

    >>> from libscientific.lda import LDA
    >>> X = [[1.0, 1.0], [1.5, 2.0], [2.0, 1.5], [6.0, 5.0], [6.5, 6.0], [7.0, 5.5]]
    >>> y = [1, 1, 1, 2, 2, 2]
    >>> model = LDA()
    >>> model.fit(X, y)
    >>> model.predict([[1.2, 1.4], [6.8, 5.9]])
    [1, 2]

The functions after the class are the low-level ctypes wrappers of ``lda.h``.
"""

from __future__ import annotations

import ctypes
from collections.abc import Sequence

import libscientific.matrix as mx
import libscientific.tensor as tns
import libscientific.vector as vect
from libscientific.loadlibrary import load_libscientific_library
from libscientific.matrix import MatrixLike

lsci = load_libscientific_library()


class LDAMODEL(ctypes.Structure):
    """
    C structure holding a fitted LDA model (low-level API).

    Use :class:`LDA` instead of accessing this structure directly. With
    ``n`` objects, ``f`` features and ``K`` classes:

    Attributes
    ----------
    inv_cov : POINTER(MATRIX)
        Inverse (or pseudo-inverse) of the pooled within-class scatter, ``f x f``.
    evect, eval : POINTER(MATRIX), POINTER(DVECTOR)
        Discriminant directions (``f x f``, one per column) and their
        eigenvalues (size ``f``).
    mu : POINTER(MATRIX)
        Class means, ``K x f``.
    pprob : POINTER(DVECTOR)
        Prior probability of each class (fraction of training objects).
    features, mnpdf : POINTER(TENSOR)
        Training objects projected on the discriminant directions, and their
        normal densities: one ``n_k x f`` matrix per class.
    fmean, fsdev : POINTER(MATRIX)
        Mean and standard deviation of the projected features per class, ``K x f``.
    classid : POINTER(UIVECTOR)
        Row indexes of the training objects, grouped by class.
    nclass, class_start : int
        Number of classes and first label (0 or 1).
    """

    _fields_ = [
        ("inv_cov", ctypes.POINTER(mx.MATRIX)),
        ("features", ctypes.POINTER(tns.TENSOR)),
        ("mnpdf", ctypes.POINTER(tns.TENSOR)),
        ("evect", ctypes.POINTER(mx.MATRIX)),
        ("mu", ctypes.POINTER(mx.MATRIX)),
        ("fmean", ctypes.POINTER(mx.MATRIX)),
        ("fsdev", ctypes.POINTER(mx.MATRIX)),
        ("yscrambling", ctypes.POINTER(mx.MATRIX)),
        ("recalculated_y", ctypes.POINTER(mx.MATRIX)),
        ("recalculated_residuals", ctypes.POINTER(mx.MATRIX)),
        ("predicted_y", ctypes.POINTER(mx.MATRIX)),
        ("predicted_residuals", ctypes.POINTER(mx.MATRIX)),
        ("eval", ctypes.POINTER(vect.DVECTOR)),
        ("pprob", ctypes.POINTER(vect.DVECTOR)),
        ("classid", ctypes.POINTER(vect.UIVECTOR)),
        ("nclass", ctypes.c_size_t),
        ("class_start", ctypes.c_size_t),
        ("roc", ctypes.POINTER(tns.TENSOR)),
        ("pr", ctypes.POINTER(tns.TENSOR)),
        ("roc_aucs", ctypes.POINTER(vect.DVECTOR)),
        ("pr_aucs", ctypes.POINTER(vect.DVECTOR)),
    ]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


def _labels_to_column(
    y_input: Sequence[int] | Sequence[Sequence[float]] | MatrixLike,
) -> MatrixLike:
    """Accept ``[0, 1, 1]`` as well as ``[[0], [1], [1]]``."""
    if isinstance(y_input, (mx.Matrix, ctypes.POINTER(mx.MATRIX))):
        return y_input
    if "numpy" in str(type(y_input)):
        y_input = y_input.tolist()  # type: ignore[union-attr]
    rows = list(y_input)  # type: ignore[arg-type]
    if rows and not hasattr(rows[0], "__len__"):
        return [[float(v)] for v in rows]
    return rows


class LDA:
    """
    Linear Discriminant Analysis classifier.

    Attributes
    ----------
    model : ctypes.POINTER(LDAMODEL)
        The underlying C model (low-level access, see :class:`LDAMODEL`).

    Notes
    -----
    - Labels must be integers ``0..K-1`` or ``1..K`` with every class
      present; predictions use the same numbering.
    - When the within-class scatter matrix is singular (e.g. more features
      than objects), its pseudo-inverse is used (Liu et al., Int. J. Patt.
      Recogn. Artif. Intell. 21, 1265, 2007).

    Examples
    --------
    >>> X = [[1.0, 1.0], [1.5, 2.0], [2.0, 1.5], [6.0, 5.0], [6.5, 6.0], [7.0, 5.5],
    ...      [1.0, 6.0], [1.5, 7.0], [2.0, 6.5]]
    >>> y = [0, 0, 0, 1, 1, 1, 2, 2, 2]
    >>> model = LDA()
    >>> model.fit(X, y)
    >>> model.nclass, model.classes
    (3, [0, 1, 2])
    >>> model.predict(X) == y
    True
    >>> [round(p, 4) for p in model.get_priors()]
    [0.3333, 0.3333, 0.3333]
    """

    def __init__(self) -> None:
        self.model = new_lda_model()

    def __repr__(self) -> str:
        return "LDA()"

    def __del__(self) -> None:
        if getattr(self, "model", None) is not None:
            del_lda_model(self.model)
            self.model = None  # type: ignore[assignment]  # already released

    @property
    def nclass(self) -> int:
        """Number of classes of the fitted model."""
        return int(self.model.contents.nclass)

    @property
    def classes(self) -> list[int]:
        """Class labels of the fitted model, e.g. ``[0, 1, 2]`` or ``[1, 2, 3]``."""
        start = int(self.model.contents.class_start)
        return list(range(start, start + self.nclass))

    def fit(self, x_input: MatrixLike, y_input: Sequence[int] | MatrixLike) -> None:
        """
        Fit the discriminant model.

        Calling ``fit`` again replaces the previous model.

        Parameters
        ----------
        x_input : MatrixLike
            Training objects, shape ``(n_objects, n_features)``.
        y_input : list[int] or MatrixLike
            Class label of each object: a flat list (``[0, 1, 1, ...]``) or a
            one-column matrix (``[[0], [1], [1], ...]``).

        Returns
        -------
        None
            The model is updated in place.
        """
        x_input_, x_owned = mx.as_c_matrix(x_input)
        try:
            y_input_, y_owned = mx.as_c_matrix(_labels_to_column(y_input))
            try:
                lda_algorithm(x_input_, y_input_, self.model)
            finally:
                if y_owned:
                    mx.del_matrix(y_input_)
        finally:
            if x_owned:
                mx.del_matrix(x_input_)

    def _predict_all(
        self, x_input: MatrixLike
    ) -> tuple[list[list[float]], list[list[float]], list[list[float]], list[int]]:
        x_input_, owned = mx.as_c_matrix(x_input)
        pfeatures = mx.init_matrix()
        scores = mx.init_matrix()
        mnpdf = mx.init_matrix()
        labels = mx.init_matrix()
        try:
            lda_prediction(x_input_, self.model, pfeatures, scores, mnpdf, labels)
            return (
                mx.matrix_to_list(pfeatures),
                mx.matrix_to_list(scores),
                mx.matrix_to_list(mnpdf),
                [int(round(row[0])) for row in mx.matrix_to_list(labels)],
            )
        finally:
            for m in (pfeatures, scores, mnpdf, labels):
                mx.del_matrix(m)
            if owned:
                mx.del_matrix(x_input_)

    def predict(self, x_input: MatrixLike) -> list[int]:
        """
        Predict the class of new objects.

        Parameters
        ----------
        x_input : MatrixLike
            Objects, shape ``(n_new_objects, n_features)``.

        Returns
        -------
        list[int]
            One label per object, in the numbering of the training labels.

        Examples
        --------
        >>> model = LDA()
        >>> model.fit([[1.0, 1.0], [1.5, 2.0], [2.0, 1.5], [6.0, 5.0], [6.5, 6.0], [7.0, 5.5]],
        ...           [1, 1, 1, 2, 2, 2])
        >>> model.predict([[0.5, 0.8], [8.0, 7.0]])
        [1, 2]
        """
        return self._predict_all(x_input)[3]

    def predict_scores(self, x_input: MatrixLike) -> list[list[float]]:
        """
        Return the linear discriminant score of each object for each class.

        The predicted class is the one with the highest score. Scores are
        log-scale values, not probabilities: only their differences matter.

        Parameters
        ----------
        x_input : MatrixLike
            Objects, shape ``(n_new_objects, n_features)``.

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_new_objects, nclass)``; column ``k`` refers
            to label ``classes[k]``.
        """
        return self._predict_all(x_input)[1]

    def transform(self, x_input: MatrixLike) -> list[list[float]]:
        """
        Project objects on the discriminant directions.

        Useful to plot the classes: the first columns (largest eigenvalues)
        separate them best.

        Parameters
        ----------
        x_input : MatrixLike
            Objects, shape ``(n_new_objects, n_features)``.

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_new_objects, n_features)``: column ``l`` is
            the projection on eigenvector ``l`` (see :meth:`get_eigenvectors`).
        """
        return self._predict_all(x_input)[0]

    def get_eigenvalues(self) -> list[float]:
        """
        Return the eigenvalues of the discriminant directions.

        Returns
        -------
        list[float]
            One value per direction (``n_features`` values); larger values
            mean better class separation. At most ``nclass - 1`` are
            significantly different from zero.
        """
        return vect.dvector_tolist(self.model.contents.eval)

    def get_eigenvectors(self) -> list[list[float]]:
        """
        Return the discriminant directions.

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_features, n_features)``: column ``l`` is the
            direction associated with eigenvalue ``l``.
        """
        return mx.matrix_to_list(self.model.contents.evect)

    def get_class_means(self) -> list[list[float]]:
        """
        Return the mean of each class in the original feature space.

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(nclass, n_features)``; row ``k`` refers to label
            ``classes[k]``.
        """
        return mx.matrix_to_list(self.model.contents.mu)

    def get_priors(self) -> list[float]:
        """
        Return the prior probability of each class.

        Returns
        -------
        list[float]
            ``nclass`` values: the fraction of training objects in each class.
        """
        return vect.dvector_tolist(self.model.contents.pprob)

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
            write_lda(filepath, self.model)
        else:
            st = lsci.SaveLDAModel(filepath.encode("utf-8"), self.model)
            if st != 0:
                raise OSError(f"SaveLDAModel failed with error code {st}")

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
        >>> X = [[1.0, 1.0], [1.5, 2.0], [2.0, 1.5], [6.0, 5.0], [6.5, 6.0], [7.0, 5.5]]
        >>> model = LDA()
        >>> model.fit(X, [0, 0, 0, 1, 1, 1])
        >>> path = os.path.join(tempfile.mkdtemp(), "lda.lsm")
        >>> model.save(path)
        >>> loaded = LDA()
        >>> loaded.load(path)
        >>> loaded.predict(X)
        [0, 0, 0, 1, 1, 1]
        """
        from libscientific.mlr import _is_binary_model

        if _is_binary_model(filepath):
            st = lsci.ImportLDAModel(filepath.encode("utf-8"), self.model)
            if st != 0:
                raise OSError(f"ImportLDAModel failed with error code {st}")
        else:
            read_lda(filepath, self.model)


# ---------------------------------------------------------------------------
# Low-level API: thin ctypes wrappers around the C functions in lda.h.
# ---------------------------------------------------------------------------

lsci.NewLDAModel.argtypes = [ctypes.POINTER(ctypes.POINTER(LDAMODEL))]
lsci.NewLDAModel.restype = None


def new_lda_model() -> ctypes._Pointer[LDAMODEL]:
    """
    Allocate an empty C LDA model (low-level API).

    Returns
    -------
    ctypes.POINTER(LDAMODEL)
        A new, empty model; release it with :func:`del_lda_model`.
    """
    mlda = ctypes.POINTER(LDAMODEL)()
    lsci.NewLDAModel(ctypes.pointer(mlda))
    return mlda


lsci.DelLDAModel.argtypes = [ctypes.POINTER(ctypes.POINTER(LDAMODEL))]
lsci.DelLDAModel.restype = None


def del_lda_model(mlda: ctypes._Pointer[LDAMODEL]) -> None:
    """
    Free a C LDA model created with :func:`new_lda_model` (low-level API).

    Parameters
    ----------
    mlda : ctypes.POINTER(LDAMODEL)
        The model to release. It must not be used afterwards.
    """
    lsci.DelLDAModel(ctypes.pointer(mlda))


lsci.LDA.argtypes = [ctypes.POINTER(mx.MATRIX), ctypes.POINTER(mx.MATRIX), ctypes.POINTER(LDAMODEL)]
lsci.LDA.restype = None


def lda_algorithm(
    x_input: ctypes._Pointer[mx.MATRIX],
    y_input: ctypes._Pointer[mx.MATRIX],
    mlda: ctypes._Pointer[LDAMODEL],
) -> None:
    """
    Fit an LDA model (low-level API, see :meth:`LDA.fit`).

    Parameters
    ----------
    x_input : ctypes.POINTER(MATRIX)
        Training objects, ``n x f``.
    y_input : ctypes.POINTER(MATRIX)
        Class labels, ``n x 1`` (``0..K-1`` or ``1..K``).
    mlda : ctypes.POINTER(LDAMODEL)
        Model created with :func:`new_lda_model`; filled with the results.
    """
    lsci.LDA(x_input, y_input, mlda)


lsci.LDAPrediction.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(LDAMODEL),
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(mx.MATRIX),
]
lsci.LDAPrediction.restype = None


def lda_prediction(
    x_input: ctypes._Pointer[mx.MATRIX],
    mlda: ctypes._Pointer[LDAMODEL],
    pfeatures: ctypes._Pointer[mx.MATRIX],
    scores: ctypes._Pointer[mx.MATRIX],
    mnpdf: ctypes._Pointer[mx.MATRIX],
    labels: ctypes._Pointer[mx.MATRIX],
) -> None:
    """
    Classify new objects (low-level API, see :meth:`LDA.predict`).

    All outputs must be empty matrices created with
    :func:`libscientific.matrix.init_matrix`.

    Parameters
    ----------
    x_input : ctypes.POINTER(MATRIX)
        Objects, ``m x f``.
    mlda : ctypes.POINTER(LDAMODEL)
        A fitted model.
    pfeatures : ctypes.POINTER(MATRIX)
        Output: projections on the discriminant directions, ``m x f``.
    scores : ctypes.POINTER(MATRIX)
        Output: discriminant score of each class, ``m x K``.
    mnpdf : ctypes.POINTER(MATRIX)
        Output: normal density of each projection in the predicted class, ``m x f``.
    labels : ctypes.POINTER(MATRIX)
        Output: predicted labels, ``m x 1``.
    """
    lsci.LDAPrediction(x_input, mlda, pfeatures, scores, mnpdf, labels)


lsci.WriteLDA.argtypes = [ctypes.c_char_p, ctypes.POINTER(LDAMODEL)]
lsci.WriteLDA.restype = None


def write_lda(dbpath: str, mlda: ctypes._Pointer[LDAMODEL]) -> None:
    """
    Write a model to a SQLite database (legacy format, low-level API).

    Parameters
    ----------
    dbpath : str
        Path of the SQLite file to create.
    mlda : ctypes.POINTER(LDAMODEL)
        The model to store.
    """
    lsci.WriteLDA(dbpath.encode("utf-8"), mlda)


lsci.ReadLDA.argtypes = [ctypes.c_char_p, ctypes.POINTER(LDAMODEL)]
lsci.ReadLDA.restype = None


def read_lda(dbpath: str, mlda: ctypes._Pointer[LDAMODEL]) -> None:
    """
    Read a model from a SQLite database (legacy format, low-level API).

    Parameters
    ----------
    dbpath : str
        Path of the SQLite file written by :func:`write_lda`.
    mlda : ctypes.POINTER(LDAMODEL)
        Model created with :func:`new_lda_model`; filled with the stored data.
    """
    lsci.ReadLDA(dbpath.encode("utf-8"), mlda)


lsci.SaveLDAModel.argtypes = [ctypes.c_char_p, ctypes.POINTER(LDAMODEL)]
lsci.SaveLDAModel.restype = ctypes.c_int
lsci.ImportLDAModel.argtypes = [ctypes.c_char_p, ctypes.POINTER(LDAMODEL)]
lsci.ImportLDAModel.restype = ctypes.c_int
