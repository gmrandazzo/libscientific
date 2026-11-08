# tensor.py libscientific python binding
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
"""Tensors (lists of matrices) stored in the C library.

A tensor is a list of layers (matrices). It is the input of the multi-way
models (:mod:`~libscientific.upca`, :mod:`~libscientific.upls`) and of
:mod:`~libscientific.cpca` (one layer per block).

    >>> from libscientific.tensor import Tensor
    >>> t = Tensor([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]])
    >>> len(t), t[1, 0, 0]
    (2, 5.0)
"""

import ctypes
from typing import Any

from libscientific import matrix as mx
from libscientific import misc
from libscientific import vector as vect
from libscientific.loadlibrary import load_libscientific_library

lsci = load_libscientific_library()


class TENSOR(ctypes.Structure):
    """C structure of a libscientific tensor (low-level API): a list of matrices.

    Attributes
    ----------
    m : POINTER(POINTER(MATRIX))
        The matrices (layers); ``m[k]`` is layer ``k``.
    order : int
        Number of layers.

    Layers can have different shapes. The multi-way models (UPCA, UPLS) need
    every layer with the same shape; CPCA uses one layer per block.
    """

    _fields_ = [("m", ctypes.POINTER(ctypes.POINTER(mx.MATRIX))), ("order", ctypes.c_size_t)]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


lsci.initTensor.argtypes = [ctypes.POINTER(ctypes.POINTER(TENSOR))]
lsci.initTensor.restype = None


def init_tensor():
    """Allocate an empty C tensor (no layers).

    Returns
    -------
    ctypes.POINTER(TENSOR)
        A new tensor; release it with :func:`del_tensor`.
    """
    tns = ctypes.POINTER(TENSOR)()
    lsci.initTensor(ctypes.pointer(tns))
    return tns


lsci.NewTensorMatrix.argtypes = [
    ctypes.POINTER(TENSOR),
    ctypes.c_size_t,
    ctypes.c_size_t,
    ctypes.c_size_t,
]
lsci.NewTensorMatrix.restype = None


def new_tensor_matrix(tns, k_indx, num_row, num_col):
    """Allocate layer ``k_indx`` of a tensor created with ``NewTensor`` (low-level).

    Parameters
    ----------
    tns : ctypes.POINTER(TENSOR)
        The tensor.
    k_indx : int
        Index of the layer to allocate.
    num_row, num_col : int
        Shape of the layer (filled with zeros).
    """
    lsci.NewTensorMatrix(tns, k_indx, num_row, num_col)


lsci.AddTensorMatrix.argtypes = [ctypes.POINTER(TENSOR), ctypes.c_size_t, ctypes.c_size_t]
lsci.AddTensorMatrix.restype = None


def add_tensor_matrix(tns, num_row, num_col):
    """Append a new layer filled with zeros.

    Parameters
    ----------
    tns : ctypes.POINTER(TENSOR)
        The tensor; its order grows by one.
    num_row, num_col : int
        Shape of the new layer.
    """
    lsci.AddTensorMatrix(tns, num_row, num_col)


lsci.NewTensor.argtypes = [ctypes.POINTER(ctypes.POINTER(TENSOR)), ctypes.c_size_t]
lsci.NewTensor.restype = None


def new_tensor(tns_input_):
    """Create a C tensor from a list of matrices or a 3D numpy array.

    Parameters
    ----------
    tns_input_ : list of list[list[float]], or numpy.ndarray
        ``[layer_0, layer_1, ...]`` where every layer is a list of rows (or a
        2D numpy array), or a 3D array of shape ``(n_layers, n_rows, n_cols)``.
        Non-numeric values are stored as the missing value.

    Returns
    -------
    ctypes.POINTER(TENSOR)
        A new tensor; release it with :func:`del_tensor`.

    Examples
    --------
    >>> t = new_tensor([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]])
    >>> t[0].order, tensor_tolist(t)[1]
    (2, [[5.0, 6.0], [7.0, 8.0]])
    >>> del_tensor(t)
    """
    if "numpy" in str(type(tns_input_)):
        return _new_tensor_from_numpy(tns_input_)
    tns_input = tns_input_
    order = len(tns_input)
    tns = ctypes.POINTER(TENSOR)()
    lsci.NewTensor(ctypes.pointer(tns), order)
    for k in range(order):
        layer = tns_input[k]
        if "numpy" in str(type(layer)):
            layer = layer.tolist()
        nrows = len(layer)
        ncols = len(layer[0]) if nrows > 0 else 0
        lsci.NewTensorMatrix(tns, k, nrows, ncols)
        data = tns.contents.m[k].contents.data
        row_type = ctypes.c_double * ncols
        nbytes = ctypes.sizeof(row_type)
        for i in range(nrows):
            try:
                # bulk copy of a numeric row
                ctypes.memmove(data[i], row_type(*map(float, layer[i])), nbytes)
            except (ValueError, TypeError):
                # non numeric values are stored as missing
                for j in range(ncols):
                    try:
                        data[i][j] = float(layer[i][j])
                    except (ValueError, TypeError):
                        data[i][j] = misc.missing_value()
    return tns


def _new_tensor_from_numpy(arr: Any) -> Any:
    """Bulk copy of a 3D numpy array ``(n_layers, n_rows, n_cols)`` into a C tensor."""
    import numpy as np  # pylint: disable=import-outside-toplevel

    a = np.ascontiguousarray(arr, dtype=np.float64)
    if a.ndim != 3:
        raise ValueError(f"expected a 3D array (layers, rows, columns), got {a.ndim}D")
    order, nrows, ncols = a.shape
    tns = ctypes.POINTER(TENSOR)()
    lsci.NewTensor(ctypes.pointer(tns), order)
    nbytes = ncols * ctypes.sizeof(ctypes.c_double)
    for k in range(order):
        lsci.NewTensorMatrix(tns, k, nrows, ncols)
        data = tns.contents.m[k].contents.data
        base = a[k].ctypes.data
        for i in range(nrows):
            ctypes.memmove(data[i], base + i * nbytes, nbytes)
    return tns


lsci.DelTensor.argtypes = [ctypes.POINTER(ctypes.POINTER(TENSOR))]
lsci.DelTensor.restype = None


def del_tensor(tns):
    """Free a C tensor and all its layers.

    Parameters
    ----------
    tns : ctypes.POINTER(TENSOR)
        The tensor to release. It must not be used afterwards.
    """
    lsci.DelTensor(ctypes.pointer(tns))


lsci.PrintTensor.argtypes = [ctypes.POINTER(TENSOR)]
lsci.PrintTensor.restype = None


def print_tensor(tns):
    """Print a C tensor to standard output.

    Parameters
    ----------
    tns : ctypes.POINTER(TENSOR)
        The tensor to print.
    """
    lsci.PrintTensor(tns)


lsci.setTensorValue.argtypes = [ctypes.POINTER(TENSOR)]
lsci.setTensorValue.restype = None


def set_tensor_value(tns, k_indx, i_indx, j_indx, val):
    """Set one element of a C tensor.

    Parameters
    ----------
    tns : ctypes.POINTER(TENSOR)
        The tensor, modified in place.
    k_indx, i_indx, j_indx : int
        Layer, row and column index (0-based).
    val : float
        The new value.
    """
    lsci.setTensorValue(tns, k_indx, i_indx, j_indx, val)


lsci.TensorAppendColumn.argtypes = [
    ctypes.POINTER(TENSOR),
    ctypes.c_size_t,
    ctypes.POINTER(vect.DVECTOR),
]
lsci.TensorAppendColumn.restype = None


def tensor_append_column(tns, k_indx, dvect):
    """Append a column to one layer of a tensor.

    Parameters
    ----------
    tns : ctypes.POINTER(TENSOR)
        The tensor.
    k_indx : int
        Index of the layer.
    dvect : POINTER(DVECTOR)
        The column, with one value per row of the layer.
    """
    lsci.TensorAppendColumn(tns, k_indx, dvect)


lsci.TensorAppendRow.argtypes = [
    ctypes.POINTER(TENSOR),
    ctypes.c_size_t,
    ctypes.POINTER(vect.DVECTOR),
]
lsci.TensorAppendRow.restype = None


def tensor_append_row(tns, k_indx, dvect):
    """Append a row to one layer of a tensor.

    Parameters
    ----------
    tns : ctypes.POINTER(TENSOR)
        The tensor.
    k_indx : int
        Index of the layer.
    dvect : POINTER(DVECTOR)
        The row, with one value per column of the layer.
    """
    lsci.TensorAppendRow(tns, k_indx, dvect)


lsci.TensorSet.argtypes = [ctypes.POINTER(TENSOR), ctypes.c_double]
lsci.TensorSet.restype = None


def tensor_set(tns, val):
    """Set every element of every layer to the same value.

    Parameters
    ----------
    tns : ctypes.POINTER(TENSOR)
        The tensor, modified in place.
    val : float
        The value to assign.
    """
    lsci.TensorSet(tns, val)


lsci.TensorCopy.argtypes = [ctypes.POINTER(TENSOR), ctypes.POINTER(ctypes.POINTER(TENSOR))]
lsci.TensorCopy.restype = None


def tensor_copy(src_tns, dst_tns):
    """Copy a C tensor into another one.

    Parameters
    ----------
    src_tns : ctypes.POINTER(TENSOR)
        Source tensor.
    dst_tns : ctypes.POINTER(TENSOR)
        Destination tensor (e.g. from :func:`init_tensor`); it receives a copy of
        every layer.
    """
    lsci.TensorCopy(src_tns, ctypes.pointer(dst_tns))


def tensor_tolist(tns: Any) -> list[list[list[float]]]:
    """Convert a C tensor into nested Python lists.

    Parameters
    ----------
    tns : ctypes.POINTER(TENSOR)
        The tensor.

    Returns
    -------
    list[list[list[float]]]
        One list of rows per layer: ``result[k][i][j]``.
    """
    tns_lst: list[list[list[float]]] = []
    for k in range(tns.contents.order):
        tns_lst.append(mx.matrix_to_list(tns.contents.m[k]))
    return tns_lst


class Tensor:
    """A tensor (list of matrices) stored in the C library.

    Parameters
    ----------
    tns_ : list of matrices, 3D numpy array, Tensor, C tensor or None
        Initial content (copied); see :func:`new_tensor`. ``None`` gives an
        empty tensor.

    Attributes
    ----------
    tns : ctypes.POINTER(TENSOR)
        The underlying C tensor, owned by this object.

    Notes
    -----
    ``t[k, i, j]`` reads/writes one element and ``t[k]`` returns layer ``k`` as a
    list of rows; out-of-range indexes raise ``IndexError``. ``len(t)`` is the
    number of layers.

    Examples
    --------
    >>> t = Tensor([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]])
    >>> len(t), t[1, 0, 1], t[0]
    (2, 6.0, [[1.0, 2.0], [3.0, 4.0]])
    """

    def __init__(self, tns_=None):
        """Create the tensor (see the class documentation)."""
        if tns_ is None:
            self.tns = init_tensor()
        elif isinstance(tns_, Tensor):
            self.tns = init_tensor()
            tensor_copy(tns_.tns, self.tns)
        elif isinstance(tns_, ctypes.POINTER(TENSOR)):
            self.tns = init_tensor()
            tensor_copy(tns_, self.tns)
        else:
            self.tns = new_tensor(tns_)

    def __del__(self):
        """Release the C tensor."""
        if getattr(self, "tns", None) is not None:
            del_tensor(self.tns)
            del self.tns
            self.tns = None

    def __repr__(self):
        return f"Tensor({self.tolist()})"

    def __str__(self):
        return f"Tensor({self.tolist()})"

    def __len__(self):
        """Return the number of layers."""
        return self.order()

    def __iter__(self):
        """Iterate over the layers, each as a list of rows."""
        return iter(self.tolist())

    def __getitem__(self, keys):
        """Return an element (``t[k, i, j]``) or a layer (``t[k]``).

        Raises
        ------
        IndexError
            If an index is out of range.
        """
        if isinstance(keys, tuple):
            k, i, j = self._check_indexes(*keys)
            return self.data_ptr().m[k].contents.data[i][j]
        if isinstance(keys, int):
            k, _, _ = self._check_indexes(keys, None, None)
            return mx.matrix_to_list(self.data_ptr().m[k])
        raise TypeError(f"Invalid index type: {type(keys)}")

    def _check_indexes(self, k: int, i: int | None, j: int | None) -> tuple[int, int, int]:
        """Validate (layer, row, column) indexes; negative values count from the end."""

        def norm(idx: int, size: int, what: str) -> int:
            n = idx + size if idx < 0 else idx
            if not 0 <= n < size:
                raise IndexError(f"{what} index {idx} out of range for size {size}")
            return n

        kk = norm(k, self.order(), "layer")
        if i is None or j is None:
            return kk, 0, 0
        return kk, norm(i, self.nrow(kk), "row"), norm(j, self.ncol(kk), "column")

    def __setitem__(self, keys, value):
        """Set an element: ``t[k, i, j] = value``.

        Raises
        ------
        IndexError
            If an index is out of range.
        """
        k, i, j = self._check_indexes(*keys)
        self.data_ptr().m[k].contents.data[i][j] = value

    def order(self):
        """Return the number of layers.

        Returns
        -------
        int
        """
        return self.data_ptr().order

    def nrow(self, k):
        """Return the number of rows of a layer.

        Parameters
        ----------
        k : int
            Index of the layer.

        Returns
        -------
        int
        """
        return self.data_ptr().m[k].contents.row

    def ncol(self, k):
        """Return the number of columns of a layer.

        Parameters
        ----------
        k : int
            Index of the layer.

        Returns
        -------
        int
        """
        return self.data_ptr().m[k].contents.col

    def data_ptr(self):
        """Return the C tensor structure (low-level access).

        Returns
        -------
        TENSOR
        """
        return self.tns.contents

    def tolist(self):
        """Return the content as nested lists ``[layer][row][column]``.

        Returns
        -------
        list[list[list[float]]]
        """
        return tensor_tolist(self.tns)

    def fromlist(self, tns_):
        """Replace the content with a list of matrices.

        Parameters
        ----------
        tns_ : list of list[list[float]]
            The new layers.
        """
        del_tensor(self.tns)
        del self.tns
        self.tns = new_tensor(tns_)

    def fromnumpy(self, npt):
        """Replace the content with a 3D numpy array.

        Parameters
        ----------
        npt : numpy.ndarray
            Array of shape ``(n_layers, n_rows, n_cols)``.
        """
        del_tensor(self.tns)
        self.tns = new_tensor(npt)

    def debug(self):
        """Print the tensor to standard output."""
        print_tensor(self.tns)
