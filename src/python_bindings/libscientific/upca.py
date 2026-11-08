# upca.py libscientific python binding
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
"""Multi-way (unfolded) Principal Component Analysis (UPCA).

UPCA analyses three-way data: ``n_objects`` objects, each described by
``n_features`` variables measured in ``n_layers`` conditions (e.g. time
points, wavelengths, batches). The data are given as a list of layers::

    X = [layer_0, layer_1, ..., layer_{K-1}]      # each n_objects x n_features

and every layer must have the same shape. Each component has one score per
object and a loading **matrix** ``n_features x n_layers`` constrained to be
the outer product of a feature profile and a layer profile (Wold, Geladi and
Esbensen, *Multi-way principal components and PLS analysis*,
J. Chemometrics 1, 41-56, 1987).

    >>> import math
    >>> from libscientific.upca import UPCA
    >>> X = [[[3.0 * (j + 1) + math.sin(0.9 * (i + 1)) * (1 + 0.5 * j) * (1 + 0.7 * k)
    ...        for j in range(4)] for i in range(15)] for k in range(3)]
    >>> model = UPCA(npc=1, scaling=1)
    >>> model.fit(X)                      # 3 layers of 15 objects x 4 features
    >>> [round(v, 4) for v in model.get_exp_variance()]
    [100.0]

The functions after the class are the low-level ctypes wrappers of ``upca.h``.
"""

from __future__ import annotations

import ctypes
from typing import Any

import libscientific.matrix as mx
import libscientific.tensor as tns
import libscientific.vectlist as vlst
import libscientific.vector as vect
from libscientific.loadlibrary import load_libscientific_library

lsci = load_libscientific_library()


class UPCAMODEL(ctypes.Structure):
    """
    C structure holding a fitted UPCA model (low-level API).

    Use :class:`UPCA` instead of accessing this structure directly. With
    ``n`` objects, ``f`` features, ``K`` layers and ``npc`` components:

    Attributes
    ----------
    scores : POINTER(MATRIX)
        Scores, ``n x npc``.
    loadings : POINTER(TENSOR)
        One ``f x K`` loading matrix per component.
    varexp : POINTER(DVECTOR)
        Percentage of variance explained by each component.
    colaverage, colscaling : POINTER(DVECTLIST)
        Column means and scalings of each layer.
    """

    _fields_ = [
        ("scores", ctypes.POINTER(mx.MATRIX)),
        ("loadings", ctypes.POINTER(tns.TENSOR)),
        ("varexp", ctypes.POINTER(vect.DVECTOR)),
        ("colaverage", ctypes.POINTER(vlst.DVECTLIST)),
        ("colscaling", ctypes.POINTER(vlst.DVECTLIST)),
    ]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


def as_c_tensor(layers: Any) -> tuple[Any, bool]:
    """
    Convert three-way data into a C tensor.

    Parameters
    ----------
    layers : list of matrices, 3D numpy array, Tensor or C tensor
        ``[layer_0, layer_1, ...]`` with every layer ``n_objects x n_columns``
        (list of rows or 2D numpy array), or a numpy array of shape
        ``(n_layers, n_objects, n_columns)``.

    Returns
    -------
    tuple[ctypes.POINTER(TENSOR), bool]
        ``(pointer, owned)``; free the pointer with
        :func:`libscientific.tensor.del_tensor` when ``owned`` is ``True``.
    """
    if isinstance(layers, tns.Tensor):
        return layers.tns, False
    if isinstance(layers, ctypes.POINTER(tns.TENSOR)):
        return layers, False
    if "numpy" in str(type(layers)):
        layers = layers.tolist()
    layers_ = [m.tolist() if "numpy" in str(type(m)) else m for m in layers]
    return tns.new_tensor(layers_), True


class UPCA:
    """
    Multi-way Principal Component Analysis model.

    Parameters
    ----------
    npc : int, default 2
        Number of components; capped to ``n_features``. After :meth:`fit`,
        ``npc`` holds the number of components actually computed.
    scaling : int, default 1
        Column preprocessing applied to every layer: all options mean-center
        each column, then 0 = no scaling, 1 = unit variance (autoscaling),
        2 = root mean square, 3 = Pareto, 4 = range, 5 = level.

    Attributes
    ----------
    npc, scaling : int
        See Parameters.
    model : ctypes.POINTER(UPCAMODEL)
        The underlying C model (low-level access, see :class:`UPCAMODEL`).

    Notes
    -----
    Because the loading matrix of each component is constrained to rank 1,
    the components are computed one after another and are not exactly
    orthogonal: with several components the explained variances can add up
    to slightly more than 100%, and the reconstruction of general data is
    approximate. Data that follow the model exactly (a sum of
    ``score x feature profile x layer profile`` terms with orthogonal
    profiles) are reconstructed exactly.
    """

    def __init__(self, npc: int = 2, scaling: int = 1) -> None:
        self.model = new_upca_model()
        self.npc = npc
        self.scaling = scaling

    def __repr__(self) -> str:
        return f"UPCA(npc={self.npc}, scaling={self.scaling})"

    def __del__(self) -> None:
        if getattr(self, "model", None) is not None:
            del_upca_model(self.model)
            self.model = None  # type: ignore[assignment]  # already released

    def fit(self, x_input: Any) -> None:
        """
        Compute the multi-way components.

        Calling ``fit`` again replaces the previous model.

        Parameters
        ----------
        x_input : list of matrices or 3D numpy array
            ``n_layers`` layers of shape ``(n_objects, n_features)`` (see
            :func:`as_c_tensor`). All layers must have the same shape.

        Returns
        -------
        None
            The model is updated in place.
        """
        x_input_, owned = as_c_tensor(x_input)
        try:
            upca_algorithm(x_input_, self.npc, self.scaling, self.model)
        finally:
            if owned:
                tns.del_tensor(x_input_)
        self.npc = int(self.model.contents.scores[0].col)

    def get_scores(self) -> list[list[float]]:
        """
        Return the scores of the training objects.

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_objects, npc)``.
        """
        return mx.matrix_to_list(self.model.contents.scores)

    def get_loadings(self) -> list[list[list[float]]]:
        """
        Return the loading matrix of each component.

        Returns
        -------
        list[list[list[float]]]
            ``npc`` matrices of shape ``(n_features, n_layers)``: element
            ``[k][j][l]`` is the weight of feature ``j`` in layer ``l`` for
            component ``k``.
        """
        return tns.tensor_tolist(self.model.contents.loadings)

    def get_exp_variance(self) -> list[float]:
        """
        Return the percentage of variance explained by each component.

        Returns
        -------
        list[float]
            ``npc`` values, not cumulative.
        """
        return vect.dvector_tolist(self.model.contents.varexp)

    def predict(self, x_input: Any) -> list[list[float]]:
        """
        Project new objects into the model.

        Parameters
        ----------
        x_input : list of matrices or 3D numpy array
            Same number of layers and features as the training data, any
            number of objects.

        Returns
        -------
        list[list[float]]
            Scores of the new objects, shape ``(n_new_objects, npc)``.
            Projecting the training data gives :meth:`get_scores`.

        Examples
        --------
        >>> import math
        >>> X = [[[3.0 * (j + 1) + math.sin(0.9 * (i + 1)) * (1 + 0.5 * j) * (1 + 0.7 * k)
        ...        for j in range(4)] for i in range(15)] for k in range(3)]
        >>> model = UPCA(npc=1)
        >>> model.fit(X)
        >>> scores = model.predict([layer[:2] for layer in X])   # first 2 objects
        >>> [round(a - b, 6) for a, b in zip(scores[0], model.get_scores()[0])]
        [0.0]
        """
        x_input_, owned = as_c_tensor(x_input)
        pscores = mx.init_matrix()
        try:
            upca_score_predictor(x_input_, self.model, self.npc, pscores)
            return mx.matrix_to_list(pscores)
        finally:
            mx.del_matrix(pscores)
            if owned:
                tns.del_tensor(x_input_)

    def reconstruct_original_tensor(
        self, npc_input: int | None = None, scores_input: Any = None
    ) -> list[list[list[float]]]:
        """
        Rebuild three-way data in the original units from scores.

        Parameters
        ----------
        npc_input : int, optional
            Number of components to use (default: all).
        scores_input : MatrixLike, optional
            Scores, shape ``(n_objects, npc)``, e.g. from :meth:`predict`.
            Default: the training scores.

        Returns
        -------
        list[list[list[float]]]
            ``n_layers`` matrices of shape ``(n_objects, n_features)``.

        Examples
        --------
        Data that follow a one-component model are rebuilt exactly:

        >>> import math
        >>> X = [[[3.0 * (j + 1) + math.sin(0.9 * (i + 1)) * (1 + 0.5 * j) * (1 + 0.7 * k)
        ...        for j in range(4)] for i in range(15)] for k in range(3)]
        >>> model = UPCA(npc=1, scaling=1)
        >>> model.fit(X)
        >>> Xr = model.reconstruct_original_tensor()
        >>> max(abs(Xr[k][i][j] - X[k][i][j])
        ...     for k in range(3) for i in range(15) for j in range(4)) < 1e-9
        True
        """
        npc_ = self.npc if npc_input is None else npc_input
        if scores_input is None:
            scores_, owned = self.model.contents.scores, False
        else:
            scores_, owned = mx.as_c_matrix(scores_input)
        x_out = tns.init_tensor()
        try:
            upca_ind_var_predictor(scores_, self.model, npc_, x_out)
            return tns.tensor_tolist(x_out)
        finally:
            tns.del_tensor(x_out)
            if owned:
                mx.del_matrix(scores_)

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
            write_upca(filepath, self.model)
        else:
            st = lsci.SaveUPCAModel(filepath.encode("utf-8"), self.model)
            if st != 0:
                raise OSError(f"SaveUPCAModel failed with error code {st}")

    def load(self, filepath: str) -> None:
        """
        Load a model saved with :meth:`save`, replacing the current one.

        After loading, :attr:`npc` is set to the number of stored components.

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
            st = lsci.ImportUPCAModel(filepath.encode("utf-8"), self.model)
            if st != 0:
                raise OSError(f"ImportUPCAModel failed with error code {st}")
        else:
            read_upca(filepath, self.model)
        self.npc = int(self.model.contents.scores[0].col)


# ---------------------------------------------------------------------------
# Low-level API: thin ctypes wrappers around the C functions in upca.h.
# ---------------------------------------------------------------------------

lsci.NewUPCAModel.argtypes = [ctypes.POINTER(ctypes.POINTER(UPCAMODEL))]
lsci.NewUPCAModel.restype = None


def new_upca_model() -> ctypes._Pointer[UPCAMODEL]:
    """
    Allocate an empty C UPCA model (low-level API).

    Returns
    -------
    ctypes.POINTER(UPCAMODEL)
        A new, empty model; release it with :func:`del_upca_model`.
    """
    mupca = ctypes.POINTER(UPCAMODEL)()
    lsci.NewUPCAModel(ctypes.pointer(mupca))
    return mupca


lsci.DelUPCAModel.argtypes = [ctypes.POINTER(ctypes.POINTER(UPCAMODEL))]
lsci.DelUPCAModel.restype = None


def del_upca_model(mupca: ctypes._Pointer[UPCAMODEL]) -> None:
    """
    Free a C UPCA model created with :func:`new_upca_model` (low-level API).

    Parameters
    ----------
    mupca : ctypes.POINTER(UPCAMODEL)
        The model to release. It must not be used afterwards.
    """
    lsci.DelUPCAModel(ctypes.pointer(mupca))


lsci.UPCA.argtypes = [
    ctypes.POINTER(tns.TENSOR),
    ctypes.c_size_t,
    ctypes.c_size_t,
    ctypes.POINTER(UPCAMODEL),
    ctypes.POINTER(ctypes.c_int),
]
lsci.UPCA.restype = None


def upca_algorithm(
    x_input: ctypes._Pointer[tns.TENSOR], npc: int, scaling: int, mupca: ctypes._Pointer[UPCAMODEL]
) -> None:
    """
    Fit a UPCA model (low-level API, see :meth:`UPCA.fit`).

    Parameters
    ----------
    x_input : ctypes.POINTER(TENSOR)
        The layers, all ``n x f``.
    npc : int
        Number of components.
    scaling : int
        Column preprocessing, 0-5.
    mupca : ctypes.POINTER(UPCAMODEL)
        Model created with :func:`new_upca_model`; filled with the results.
    """
    lsci.UPCA(x_input, npc, scaling, mupca, None)


lsci.UPCAScorePredictor.argtypes = [
    ctypes.POINTER(tns.TENSOR),
    ctypes.POINTER(UPCAMODEL),
    ctypes.c_size_t,
    ctypes.POINTER(mx.MATRIX),
]
lsci.UPCAScorePredictor.restype = None


def upca_score_predictor(
    x_input: ctypes._Pointer[tns.TENSOR],
    mupca: ctypes._Pointer[UPCAMODEL],
    npc: int,
    pscores: ctypes._Pointer[mx.MATRIX],
) -> None:
    """
    Project new objects into a UPCA model (low-level API, see :meth:`UPCA.predict`).

    Parameters
    ----------
    x_input : ctypes.POINTER(TENSOR)
        The layers of the new objects.
    mupca : ctypes.POINTER(UPCAMODEL)
        A fitted model.
    npc : int
        Number of components to use.
    pscores : ctypes.POINTER(MATRIX)
        Output matrix created with :func:`libscientific.matrix.init_matrix`;
        filled with the scores, ``m x npc``.
    """
    lsci.UPCAScorePredictor(x_input, mupca, npc, pscores)


lsci.UPCAIndVarPredictor.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(tns.TENSOR),
    ctypes.POINTER(vlst.DVECTLIST),
    ctypes.POINTER(vlst.DVECTLIST),
    ctypes.c_size_t,
    ctypes.POINTER(tns.TENSOR),
]
lsci.UPCAIndVarPredictor.restype = None


def upca_ind_var_predictor(
    scores: ctypes._Pointer[mx.MATRIX],
    mupca: ctypes._Pointer[UPCAMODEL],
    npc: int,
    x_out: ctypes._Pointer[tns.TENSOR],
) -> None:
    """
    Rebuild original-space layers from scores (low-level API).

    See :meth:`UPCA.reconstruct_original_tensor`.

    Parameters
    ----------
    scores : ctypes.POINTER(MATRIX)
        Scores, ``n x >= npc``.
    mupca : ctypes.POINTER(UPCAMODEL)
        A fitted model.
    npc : int
        Number of components to use.
    x_out : ctypes.POINTER(TENSOR)
        Output tensor created with :func:`libscientific.tensor.init_tensor`;
        receives one ``n x f`` matrix per layer.
    """
    lsci.UPCAIndVarPredictor(
        scores,
        mupca.contents.loadings,
        mupca.contents.colaverage,
        mupca.contents.colscaling,
        npc,
        x_out,
    )


lsci.WriteUPCA.argtypes = [ctypes.c_char_p, ctypes.POINTER(UPCAMODEL)]
lsci.WriteUPCA.restype = None


def write_upca(dbpath: str, mupca: ctypes._Pointer[UPCAMODEL]) -> None:
    """
    Write a model to a SQLite database (legacy format, low-level API).

    Parameters
    ----------
    dbpath : str
        Path of the SQLite file to create.
    mupca : ctypes.POINTER(UPCAMODEL)
        The model to store.
    """
    lsci.WriteUPCA(dbpath.encode("utf-8"), mupca)


lsci.ReadUPCA.argtypes = [ctypes.c_char_p, ctypes.POINTER(UPCAMODEL)]
lsci.ReadUPCA.restype = None


def read_upca(dbpath: str, mupca: ctypes._Pointer[UPCAMODEL]) -> None:
    """
    Read a model from a SQLite database (legacy format, low-level API).

    Parameters
    ----------
    dbpath : str
        Path of the SQLite file written by :func:`write_upca`.
    mupca : ctypes.POINTER(UPCAMODEL)
        Model created with :func:`new_upca_model`; filled with the stored data.
    """
    lsci.ReadUPCA(dbpath.encode("utf-8"), mupca)


lsci.SaveUPCAModel.argtypes = [ctypes.c_char_p, ctypes.POINTER(UPCAMODEL)]
lsci.SaveUPCAModel.restype = ctypes.c_int
lsci.ImportUPCAModel.argtypes = [ctypes.c_char_p, ctypes.POINTER(UPCAMODEL)]
lsci.ImportUPCAModel.restype = ctypes.c_int
