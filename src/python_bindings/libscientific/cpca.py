# cpca.py libscientific python binding
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
"""Consensus Principal Component Analysis (CPCA), a multi-block PCA.

CPCA analyses several **blocks** of variables measured on the same
objects (for example spectra from different instruments, or different
groups of descriptors). Every block is a matrix with the same rows (objects)
and its own columns (features). The model finds:

- **super scores** (``n_objects x npc``): a consensus description of the
  objects common to all the blocks. They equal the PCA scores of the
  concatenated blocks, each block scaled by ``1/sqrt(n_columns)``;
- **super weights** (``n_blocks x npc``): how much each block contributes
  to each component;
- **block scores** and **block loadings**: the description of each block.

Algorithm: Westerhuis, Kourti and MacGregor, *Analysis of multiblock and
hierarchical PCA and PLS models*, J. Chemometrics 12, 301-321 (1998).

Most users only need the :class:`CPCA` class::

    >>> from libscientific.cpca import CPCA
    >>> block1 = [[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [3.0, 4.0, 2.0],
    ...           [4.0, 3.0, 1.0], [5.0, 5.0, 0.5], [6.0, 4.0, 2.5]]
    >>> block2 = [[0.5, 1.0], [1.5, 0.5], [1.0, 2.5],
    ...           [2.5, 1.5], [2.0, 3.0], [3.5, 2.0]]
    >>> model = CPCA(scaling=1, npc=2).fit([block1, block2])
    >>> [round(v, 2) for v in model.get_total_exp_variance()]
    [70.31, 24.68]

The functions after the class are the low-level ctypes wrappers of the C API.
"""

from __future__ import annotations

import ctypes
from collections.abc import Sequence
from typing import Any

import libscientific.matrix as mx
import libscientific.tensor as tns
import libscientific.vectlist as vlst
import libscientific.vector as vect
from libscientific.loadlibrary import load_libscientific_library
from libscientific.matrix import MatrixLike

lsci = load_libscientific_library()


class CPCAMODEL(ctypes.Structure):
    """
    C structure holding a fitted CPCA model (low-level API).

    Use :class:`CPCA` instead of accessing this structure directly. With
    ``n`` objects, ``B`` blocks (block ``b`` has ``f_b`` columns) and ``npc``
    components:

    Attributes
    ----------
    block_scores : POINTER(TENSOR)
        One ``n x B`` matrix per component: column ``b`` holds the block
        scores of block ``b``.
    block_loadings : POINTER(TENSOR)
        One ``f_b x npc`` matrix per block.
    super_scores : POINTER(MATRIX)
        Consensus scores, ``n x npc``.
    super_weights : POINTER(MATRIX)
        Block weights, ``B x npc`` (each column has unit length).
    scaling_factor : POINTER(DVECTOR)
        Block scaling factor ``sqrt(f_b)`` of each block.
    total_expvar : POINTER(DVECTOR)
        Percentage of the total variance explained by each component.
    block_expvar : POINTER(DVECTLIST)
        For each component, the cumulative percentage of the variance of each
        block explained by the components up to that one.
    colaverage, colscaling : POINTER(DVECTLIST)
        Column means and scalings of each block.
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
        ("colscaling", ctypes.POINTER(vlst.DVECTLIST)),
    ]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


def _as_c_tensor(blocks: Any) -> tuple[Any, bool]:
    """Convert blocks (list of matrices, Tensor or C tensor) to a C tensor.

    Returns ``(pointer, owned)``; free the pointer with
    :func:`libscientific.tensor.del_tensor` when ``owned`` is ``True``.
    """
    if isinstance(blocks, tns.Tensor):
        return blocks.tns, False
    if isinstance(blocks, ctypes.POINTER(tns.TENSOR)):
        return blocks, False
    blocks_ = [b.tolist() if "numpy" in str(type(b)) else b for b in blocks]
    return tns.new_tensor(blocks_), True


class CPCA:
    """
    Consensus Principal Component Analysis model.

    Parameters
    ----------
    scaling : int, default 1
        Column preprocessing applied to every block. All options mean-center
        the columns, then: 0 = no scaling, 1 = unit variance (autoscaling),
        2 = root mean square, 3 = Pareto, 4 = range, 5 = level (see
        :class:`~libscientific.pca.PCA`). After the column scaling each block
        is also divided by ``sqrt(n_columns)`` so that large blocks do not
        dominate.
    npc : int, default 2
        Number of components. It is capped to the number of columns of the
        smallest block: after :meth:`fit`, ``npc`` holds the number of
        components actually computed.

    Attributes
    ----------
    scaling : int
        The preprocessing option given at construction.
    npc : int
        Number of components of the model.
    model : ctypes.POINTER(CPCAMODEL)
        The underlying C model (low-level access, see :class:`CPCAMODEL`).

    Notes
    -----
    - The input is a list of blocks. Every block is a matrix (list of rows or
      2D numpy array) and all the blocks must have the same number of rows,
      in the same object order. A :class:`~libscientific.tensor.Tensor` is
      also accepted.
    - The sign of each component is arbitrary.

    Examples
    --------
    >>> block1 = [[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [3.0, 4.0, 2.0],
    ...           [4.0, 3.0, 1.0], [5.0, 5.0, 0.5], [6.0, 4.0, 2.5]]
    >>> block2 = [[0.5, 1.0], [1.5, 0.5], [1.0, 2.5],
    ...           [2.5, 1.5], [2.0, 3.0], [3.5, 2.0]]
    >>> model = CPCA(scaling=1, npc=2)
    >>> model.fit([block1, block2])
    CPCA(scaling=1, npc=2)
    >>> len(model.get_super_scores()), len(model.get_super_scores()[0])  # 6 objects x 2 PCs
    (6, 2)
    >>> [[round(v, 2) for v in pc] for pc in model.get_block_expvar()]   # cumulative % per block
    [[80.09, 60.53], [91.57, 98.4]]
    >>> super_scores, block_scores = model.predict([block1[:2], block2[:2]])  # 2 new objects
    >>> len(super_scores), len(super_scores[0])
    (2, 2)
    """

    def __init__(self, scaling: int = 1, npc: int = 2) -> None:
        self.model = new_cpca_model()
        self.scaling = scaling
        self.npc = npc

    def __repr__(self) -> str:
        return f"CPCA(scaling={self.scaling}, npc={self.npc})"

    def __str__(self) -> str:
        return self.__repr__()

    def __del__(self) -> None:
        if getattr(self, "model", None) is not None:
            del_cpca_model(self.model)
            self.model = None  # type: ignore[assignment]  # already released

    def fit(self, t_input: Sequence[MatrixLike] | Any) -> CPCA:
        """
        Compute the consensus components of a set of blocks.

        Calling ``fit`` again replaces the previous model.

        Parameters
        ----------
        t_input : list of MatrixLike, or Tensor
            The blocks: ``[block_0, block_1, ...]``, block ``b`` with shape
            ``(n_objects, n_columns_b)``. All blocks have the same rows.

        Returns
        -------
        CPCA
            The fitted model itself, so calls can be chained:
            ``CPCA(npc=2).fit(blocks)``.
        """
        t_input_, owned = _as_c_tensor(t_input)
        try:
            cpca_algorithm(t_input_, self.scaling, self.npc, self.model)
        finally:
            if owned:
                tns.del_tensor(t_input_)
        self.npc = int(self.model.contents.super_scores[0].col)
        return self

    def get_super_scores(self) -> list[list[float]]:
        """
        Return the super (consensus) scores of the training objects.

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_objects, npc)``: the position of each object
            in the space shared by all the blocks.
        """
        return mx.matrix_to_list(self.model.contents.super_scores)

    def get_super_weights(self) -> list[list[float]]:
        """
        Return the super weights: the contribution of each block to each component.

        Returns
        -------
        list[list[float]]
            Matrix of shape ``(n_blocks, npc)``. Each column has unit length;
            a large value means that the block drives that component.

        Examples
        --------
        >>> block1 = [[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [3.0, 4.0, 2.0],
        ...           [4.0, 3.0, 1.0], [5.0, 5.0, 0.5], [6.0, 4.0, 2.5]]
        >>> block2 = [[0.5, 1.0], [1.5, 0.5], [1.0, 2.5],
        ...           [2.5, 1.5], [2.0, 3.0], [3.5, 2.0]]
        >>> w = CPCA(npc=2).fit([block1, block2]).get_super_weights()
        >>> len(w), len(w[0])            # 2 blocks x 2 components
        (2, 2)
        >>> round(w[0][0] ** 2 + w[1][0] ** 2, 6)
        1.0
        """
        return mx.matrix_to_list(self.model.contents.super_weights)

    def get_block_scores(self) -> list[list[list[float]]]:
        """
        Return the block scores of the training objects, grouped by component.

        Returns
        -------
        list[list[list[float]]]
            ``npc`` matrices, one per component. Matrix ``k`` has shape
            ``(n_objects, n_blocks)``: column ``b`` holds the scores of block
            ``b`` on component ``k``. So ``get_block_scores()[k][i][b]`` is
            the score of object ``i`` in block ``b`` on component ``k``.
        """
        return tns.tensor_tolist(self.model.contents.block_scores)

    def get_block_loadings(self) -> list[list[list[float]]]:
        """
        Return the loadings of each block.

        Returns
        -------
        list[list[list[float]]]
            ``n_blocks`` matrices, one per block. Matrix ``b`` has shape
            ``(n_columns_b, npc)``: the weight of each feature of block ``b``
            in each component.
        """
        return tns.tensor_tolist(self.model.contents.block_loadings)

    def get_block_expvar(self) -> list[list[float]]:
        """
        Return the variance of each block explained by the model, in percent.

        Returns
        -------
        list[list[float]]
            ``npc`` rows, one per component, with ``n_blocks`` values each:
            element ``[k][b]`` is the **cumulative** percentage of the variance
            of block ``b`` explained by components ``0..k``.
        """
        return vlst.dvector_list_tolist(self.model.contents.block_expvar)

    def get_total_exp_variance(self) -> list[float]:
        """
        Return the percentage of the total variance explained by each component.

        Returns
        -------
        list[float]
            ``npc`` values, one per component (not cumulative).
        """
        return vect.dvector_tolist(self.model.contents.total_expvar)

    def predict(
        self, t_input: Sequence[MatrixLike] | Any
    ) -> tuple[list[list[float]], list[list[list[float]]]]:
        """
        Project new objects into the model.

        Parameters
        ----------
        t_input : list of MatrixLike, or Tensor
            The blocks of the new objects, with the same number of blocks and
            columns, in the same order, as the training data. Every block has
            ``n_new_objects`` rows.

        Returns
        -------
        super_scores : list[list[float]]
            Shape ``(n_new_objects, npc)``.
        block_scores : list[list[list[float]]]
            ``npc`` matrices of shape ``(n_new_objects, n_blocks)``, same layout
            as :meth:`get_block_scores`.

        Notes
        -----
        Projecting the training objects reproduces :meth:`get_super_scores`
        within the convergence tolerance of the algorithm (about 1e-5).
        """
        t_input_, owned = _as_c_tensor(t_input)
        p_super_scores_ = mx.init_matrix()
        p_block_scores_ = tns.init_tensor()
        try:
            cpca_score_predictor(t_input_, self.model, self.npc, p_super_scores_, p_block_scores_)
            return mx.matrix_to_list(p_super_scores_), tns.tensor_tolist(p_block_scores_)
        finally:
            mx.del_matrix(p_super_scores_)
            tns.del_tensor(p_block_scores_)
            if owned:
                tns.del_tensor(t_input_)

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
            write_cpca(filepath, self.model)
        else:
            from libscientific.model_io import SaveCPCAModel

            st = SaveCPCAModel(filepath.encode("utf-8"), self.model)
            if st != 0:
                raise OSError(f"SaveCPCAModel failed with error code {st}")

    def load(self, filepath: str) -> None:
        """
        Load a model saved with :meth:`save`, replacing the current one.

        The format is detected from the file content or the extension. After
        loading, :attr:`npc` is set to the number of stored components.

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
        >>> block1 = [[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [3.0, 4.0, 2.0],
        ...           [4.0, 3.0, 1.0], [5.0, 5.0, 0.5], [6.0, 4.0, 2.5]]
        >>> block2 = [[0.5, 1.0], [1.5, 0.5], [1.0, 2.5],
        ...           [2.5, 1.5], [2.0, 3.0], [3.5, 2.0]]
        >>> model = CPCA(npc=2).fit([block1, block2])
        >>> path = os.path.join(tempfile.mkdtemp(), "cpca.lsm")
        >>> model.save(path)
        >>> loaded = CPCA()
        >>> loaded.load(path)
        >>> loaded.npc, loaded.predict([block1, block2]) == model.predict([block1, block2])
        (2, True)
        """
        is_lscm = False
        try:
            with open(filepath, "rb") as f:
                if f.read(4) == b"LSCM":
                    is_lscm = True
        except OSError:
            pass

        if is_lscm or filepath.endswith(".lsm") or filepath.endswith(".lsb"):
            from libscientific.model_io import ImportCPCAModel

            st = ImportCPCAModel(filepath.encode("utf-8"), self.model)
            if st != 0:
                raise OSError(f"ImportCPCAModel failed with error code {st}")
        else:
            read_cpca(filepath, self.model)
        self.npc = int(self.model.contents.super_scores[0].col)


# ---------------------------------------------------------------------------
# Low-level API: thin ctypes wrappers around the C functions in cpca.h.
# ---------------------------------------------------------------------------

lsci.NewCPCAModel.argtypes = [ctypes.POINTER(ctypes.POINTER(CPCAMODEL))]
lsci.NewCPCAModel.restype = None


def new_cpca_model() -> ctypes._Pointer[CPCAMODEL]:
    """
    Allocate an empty C CPCA model (low-level API).

    Returns
    -------
    ctypes.POINTER(CPCAMODEL)
        A new, empty model; release it with :func:`del_cpca_model`.
    """
    mpca = ctypes.POINTER(CPCAMODEL)()
    lsci.NewCPCAModel(ctypes.pointer(mpca))
    return mpca


lsci.DelCPCAModel.argtypes = [ctypes.POINTER(ctypes.POINTER(CPCAMODEL))]
lsci.DelCPCAModel.restype = None


def del_cpca_model(mpca: ctypes._Pointer[CPCAMODEL]) -> None:
    """
    Free a C CPCA model created with :func:`new_cpca_model` (low-level API).

    Parameters
    ----------
    mpca : ctypes.POINTER(CPCAMODEL)
        The model to release. It must not be used afterwards.
    """
    lsci.DelCPCAModel(ctypes.pointer(mpca))


lsci.CPCA.argtypes = [
    ctypes.POINTER(tns.TENSOR),
    ctypes.c_int,
    ctypes.c_size_t,
    ctypes.POINTER(CPCAMODEL),
]
lsci.CPCA.restype = None


def cpca_algorithm(
    t_input: ctypes._Pointer[tns.TENSOR], scaling: int, npc: int, mpca: ctypes._Pointer[CPCAMODEL]
) -> None:
    """
    Fit a CPCA model (low-level API, see :meth:`CPCA.fit`).

    Parameters
    ----------
    t_input : ctypes.POINTER(TENSOR)
        The blocks, one tensor matrix per block, all with the same rows.
    scaling : int
        Column preprocessing, 0-5 (see :class:`CPCA`).
    npc : int
        Number of components (capped to the smallest block width).
    mpca : ctypes.POINTER(CPCAMODEL)
        Model created with :func:`new_cpca_model`; filled with the results.
    """
    lsci.CPCA(t_input, scaling, npc, mpca)


lsci.CPCAScorePredictor.argtypes = [
    ctypes.POINTER(tns.TENSOR),
    ctypes.POINTER(CPCAMODEL),
    ctypes.c_size_t,
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(tns.TENSOR),
]
lsci.CPCAScorePredictor.restype = None


def cpca_score_predictor(
    t_input: ctypes._Pointer[tns.TENSOR],
    model: ctypes._Pointer[CPCAMODEL],
    npc: int,
    p_super_scores: ctypes._Pointer[mx.MATRIX],
    p_block_scores: ctypes._Pointer[tns.TENSOR],
) -> None:
    """
    Project new objects into a CPCA model (low-level API, see :meth:`CPCA.predict`).

    Parameters
    ----------
    t_input : ctypes.POINTER(TENSOR)
        The blocks of the new objects.
    model : ctypes.POINTER(CPCAMODEL)
        A fitted model.
    npc : int
        Number of components to use.
    p_super_scores : ctypes.POINTER(MATRIX)
        Output matrix created with :func:`libscientific.matrix.init_matrix`;
        filled with the super scores, ``n x npc``.
    p_block_scores : ctypes.POINTER(TENSOR)
        Output tensor created with :func:`libscientific.tensor.init_tensor`;
        filled with one ``n x n_blocks`` matrix per component.
    """
    lsci.CPCAScorePredictor(t_input, model, npc, p_super_scores, p_block_scores)


lsci.PrintCPCA.argtypes = [ctypes.POINTER(CPCAMODEL)]
lsci.PrintCPCA.restype = None


def print_cpca(mpca: ctypes._Pointer[CPCAMODEL]) -> None:
    """
    Print the content of a CPCA model to standard output (low-level API).

    Parameters
    ----------
    mpca : ctypes.POINTER(CPCAMODEL)
        The model to print (e.g. ``CPCA().model``).
    """
    lsci.PrintCPCA(mpca)


lsci.WriteCPCA.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.POINTER(CPCAMODEL)]
lsci.WriteCPCA.restype = None


def write_cpca(dbpath: str, cpca: ctypes._Pointer[CPCAMODEL]) -> None:
    """
    Write a model to a SQLite database (legacy format, low-level API).

    Prefer :meth:`CPCA.save`.

    Parameters
    ----------
    dbpath : str
        Path of the SQLite file to create.
    cpca : ctypes.POINTER(CPCAMODEL)
        The model to store.
    """
    encoded_string = dbpath.encode("utf-8")
    ctypes_string = ctypes.c_char_p(encoded_string)
    lsci.WriteCPCA(ctypes_string, cpca)


lsci.ReadCPCA.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.POINTER(CPCAMODEL)]
lsci.ReadCPCA.restype = None


def read_cpca(dbpath: str, cpca: ctypes._Pointer[CPCAMODEL]) -> None:
    """
    Read a model from a SQLite database (legacy format, low-level API).

    Prefer :meth:`CPCA.load`.

    Parameters
    ----------
    dbpath : str
        Path of the SQLite file written by :func:`write_cpca`.
    cpca : ctypes.POINTER(CPCAMODEL)
        Model created with :func:`new_cpca_model`; filled with the stored data.
    """
    encoded_string = dbpath.encode("utf-8")
    ctypes_string = ctypes.c_char_p(encoded_string)
    lsci.ReadCPCA(ctypes_string, cpca)
