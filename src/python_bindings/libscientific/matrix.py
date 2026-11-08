# matrix.py libscientific python binding
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
"""Matrices stored in the C library.

- :class:`Matrix`: Pythonic matrix (indexing, shape, transpose, SVD, ...).
- :data:`MatrixLike` and :func:`as_c_matrix`: the inputs accepted by all the
  models (lists of rows, numpy arrays, :class:`Matrix`).
- Low-level functions working on C matrix pointers (``new_matrix``,
  ``matrix_to_list``, linear algebra, ...), mainly for interoperability.

    >>> from libscientific.matrix import Matrix
    >>> m = Matrix([[1.0, 2.0], [3.0, 4.0]])
    >>> m.shape, m[1, 0]
    ((2, 2), 3.0)
"""

from __future__ import annotations

import csv
import ctypes
from typing import TYPE_CHECKING, Any, TypeAlias

from libscientific import misc, vector
from libscientific.loadlibrary import load_libscientific_library

if TYPE_CHECKING:  # only for type annotations; numpy is optional at runtime
    from collections.abc import Sequence

    import numpy as np

lsci = load_libscientific_library()

MatrixLike: TypeAlias = "Matrix | Sequence[Sequence[float]] | np.ndarray | ctypes._Pointer[MATRIX]"
"""Any 2D input accepted by the libscientific models and matrix functions.

One of:

- a list of rows (``list[list[float]]``), e.g. ``[[1.0, 2.0], [3.0, 4.0]]``;
- a 2D :class:`numpy.ndarray` (any dtype convertible to ``float64``);
- a :class:`Matrix`;
- a raw C matrix pointer (``ctypes.POINTER(MATRIX)``) from the low-level API.

Rows are objects (samples) and columns are features (variables).
"""


class MATRIX(ctypes.Structure):
    """C structure of a libscientific matrix (low-level API).

    Attributes
    ----------
    data : POINTER(POINTER(c_double))
        Row pointers: ``data[i][j]`` is the element in row ``i``, column ``j``.
    row, col : int
        Number of rows and columns.

    Use :class:`Matrix` or the functions of this module instead of accessing it
    directly.
    """

    _fields_ = [
        ("data", ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
        ("row", ctypes.c_size_t),
        ("col", ctypes.c_size_t),
    ]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


lsci.initMatrix.argtypes = [ctypes.POINTER(ctypes.POINTER(MATRIX))]
lsci.initMatrix.restype = None


def init_matrix():
    """Allocate an empty (0 x 0) C matrix.

    Returns
    -------
    ctypes.POINTER(MATRIX)
        A new matrix; release it with :func:`del_matrix`. Functions that fill an
        output matrix (e.g. :func:`matrix_inversion`) resize it as needed.

    Examples
    --------
    >>> m = init_matrix()
    >>> (m[0].row, m[0].col)
    (0, 0)
    >>> del_matrix(m)
    """
    mtx = ctypes.POINTER(MATRIX)()
    # m = matrix()
    lsci.initMatrix(ctypes.pointer(mtx))
    return mtx


lsci.NewMatrix.argtypes = [ctypes.POINTER(ctypes.POINTER(MATRIX)), ctypes.c_size_t, ctypes.c_size_t]
lsci.NewMatrix.restype = None


def _new_matrix_from_list_slow(mx_lst, nrows, ncols):
    """Element-by-element import (ragged rows, non numeric values as MISSING)."""
    mtx = ctypes.POINTER(MATRIX)()
    lsci.NewMatrix(ctypes.pointer(mtx), nrows, ncols)

    for i in range(nrows):
        for j in range(ncols):
            val = None
            try:
                val = float(mx_lst[i][j])
            except TypeError:
                val = float(mx_lst[i])
            except ValueError:
                val = None

            if val is None:
                set_missing_matrix_value(mtx, i, j)
            else:
                lsci.setMatrixValue(mtx, i, j, val)
    return mtx


def new_matrix_from_list(mx_lst):
    """Create a C matrix from a list of rows.

    Non-numeric values (e.g. ``"NA"``) are stored as the missing value (see
    :func:`libscientific.misc.missing_value`). A flat list becomes a single
    column.

    Parameters
    ----------
    mx_lst : list[list[float]] or list[float]
        The rows of the matrix.

    Returns
    -------
    ctypes.POINTER(MATRIX)
        A new matrix; release it with :func:`del_matrix`.

    Examples
    --------
    >>> m = new_matrix_from_list([[1, 2], [3, "NA"]])
    >>> matrix_to_list(m)[0], matrix_to_list(m)[1][1] > 1e7
    ([1.0, 2.0], True)
    >>> del_matrix(m)
    """
    nrows = 0
    ncols = 0
    try:
        nrows = len(mx_lst)
        try:
            ncols = len(mx_lst[0])
        except TypeError:
            ncols = 1
        except IndexError:
            ncols = 0
    except (IndexError, TypeError):
        nrows = 0

    if nrows > 0 and ncols > 0:
        # Fast path 1: let numpy parse the whole list in C, then bulk-copy.
        try:
            import numpy as np  # pylint: disable=import-outside-toplevel

            arr = np.asarray(mx_lst, dtype=np.float64)
            if arr.ndim == 1:
                arr = arr.reshape(-1, 1)
            if arr.shape == (nrows, ncols):
                return matrix_from_numpy(arr)
        except (ImportError, ValueError, TypeError):
            pass

        # Fast path 2 (no numpy): convert one row at a time and copy it with
        # a single memmove instead of one ctypes call per element.
        try:
            mtx = ctypes.POINTER(MATRIX)()
            lsci.NewMatrix(ctypes.pointer(mtx), nrows, ncols)
            row_type = ctypes.c_double * ncols
            nbytes = ctypes.sizeof(row_type)
            data = mtx[0].data
            try:
                for i, row in enumerate(mx_lst):
                    if ncols == 1 and not hasattr(row, "__len__"):
                        row = (row,)
                    if len(row) != ncols:
                        raise ValueError("ragged row")
                    ctypes.memmove(data[i], row_type(*map(float, row)), nbytes)
                return mtx
            except (ValueError, TypeError):
                del_matrix(mtx)
        except (ValueError, TypeError):
            pass

    return _new_matrix_from_list_slow(mx_lst, nrows, ncols)


def new_matrix_from_csv(fcsv):
    """Create a C matrix from a CSV file (comma separated, one row per line).

    Parameters
    ----------
    fcsv : str
        Path of the CSV file. Every row must have the same number of values.

    Returns
    -------
    ctypes.POINTER(MATRIX)
        A new matrix; release it with :func:`del_matrix`.
    """
    mtx = init_matrix()
    with open(fcsv, encoding="utf-8") as fcsv_:
        reader = csv.reader(fcsv_, delimiter=",", quotechar='"')
        for csvrow in reader:
            row_lst = vector.DVector(csvrow)
            matrix_append_row(mtx, row_lst)
            del row_lst
    return mtx


def new_matrix(mtx_input_):
    """Create a C matrix from a list of rows or from a CSV file.

    Parameters
    ----------
    mtx_input_ : list[list[float]] or str
        The rows of the matrix, or the path of a CSV file.

    Returns
    -------
    ctypes.POINTER(MATRIX)
        A new matrix; release it with :func:`del_matrix`. To accept numpy arrays
        and :class:`Matrix` objects too, use :func:`as_c_matrix`.

    Examples
    --------
    >>> m = new_matrix([[1.0, 2.0], [3.0, 4.0]])
    >>> matrix_to_list(m)
    [[1.0, 2.0], [3.0, 4.0]]
    >>> del_matrix(m)
    """
    if "str" in str(type(mtx_input_)):
        return new_matrix_from_csv(mtx_input_)
    return new_matrix_from_list(mtx_input_)


def as_c_matrix(m_input: MatrixLike) -> tuple[ctypes._Pointer[MATRIX], bool]:
    """
    Convert any :data:`MatrixLike` input into a C matrix pointer.

    Parameters
    ----------
    m_input : MatrixLike
        A list of rows, a 2D numpy array, a :class:`Matrix` or a C matrix pointer.

    Returns
    -------
    tuple[ctypes.POINTER(MATRIX), bool]
        ``(pointer, owned)``. When ``owned`` is ``True`` the matrix was
        allocated by this call and the caller must release it with
        :func:`del_matrix`; when ``False`` the pointer belongs to ``m_input``
        and must not be freed.

    Examples
    --------
    >>> ptr, owned = as_c_matrix([[1.0, 2.0], [3.0, 4.0]])
    >>> (ptr[0].row, ptr[0].col, owned)
    (2, 2, True)
    >>> del_matrix(ptr)
    """
    if isinstance(m_input, Matrix):
        return m_input.mtx, False
    if isinstance(m_input, ctypes.POINTER(MATRIX)):
        return m_input, False
    if "numpy" in str(type(m_input)):
        return matrix_from_numpy(m_input), True
    return new_matrix(m_input), True


lsci.ResizeMatrix.argtypes = [ctypes.POINTER(MATRIX), ctypes.c_size_t, ctypes.c_size_t]
lsci.ResizeMatrix.restype = None


def resize_matrix(mtx, nrows, ncols):
    """Resize a C matrix; all its values are set to zero.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        The matrix to resize.
    nrows, ncols : int
        New shape.
    """
    lsci.ResizeMatrix(mtx, nrows, ncols)


lsci.DelMatrix.argtypes = [ctypes.POINTER(ctypes.POINTER(MATRIX))]
lsci.DelMatrix.restype = None


def del_matrix(mtx):
    """Free a C matrix created by this module.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        The matrix to release. It must not be used afterwards.
    """
    lsci.DelMatrix(ctypes.pointer(mtx))


lsci.MatrixCheck.argtypes = [ctypes.POINTER(MATRIX)]
lsci.MatrixCheck.restype = None


def matrix_check(mtx):
    """Replace NaN and infinite values of a C matrix with the missing value.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        The matrix, modified in place.

    Examples
    --------
    >>> m = new_matrix([[1.0, float("nan")]])
    >>> matrix_check(m)
    >>> matrix_to_list(m)[0][1] == misc.missing_value()
    True
    >>> del_matrix(m)
    """
    lsci.MatrixCheck(mtx)


lsci.PrintMatrix.argtypes = [ctypes.POINTER(MATRIX)]
lsci.PrintMatrix.restype = None


def print_matrix(mtx):
    """Print a C matrix to standard output.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        The matrix to print.
    """
    lsci.PrintMatrix(mtx)


lsci.ValInMatrix.argtypes = [ctypes.POINTER(MATRIX), ctypes.c_double]
lsci.ValInMatrix.restype = ctypes.c_int


def val_in_matrix(mtx, val):
    """Check whether a C matrix contains a value (exact comparison).

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        The matrix to search.
    val : float
        The value to look for.

    Returns
    -------
    int
        1 if the value is present, 0 otherwise.
    """
    return lsci.ValInMatrix(mtx, val)


lsci.MatrixSet.argtypes = [ctypes.POINTER(MATRIX), ctypes.c_double]
lsci.MatrixSet.restype = None


def matrix_set(mtx, val):
    """Set every element of a C matrix to the same value.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        The matrix, modified in place.
    val : float
        The value to assign.
    """
    lsci.MatrixSet(mtx, val)


lsci.MatrixCopy.argtypes = [ctypes.POINTER(MATRIX), ctypes.POINTER(ctypes.POINTER(MATRIX))]
lsci.MatrixCopy.restype = None


def matrix_copy(msrc, mdst):
    """Copy a C matrix into another one.

    Parameters
    ----------
    msrc : ctypes.POINTER(MATRIX)
        Source matrix.
    mdst : ctypes.POINTER(MATRIX)
        Destination matrix (e.g. from :func:`init_matrix`); it is resized to the
        shape of ``msrc``.
    """
    lsci.MatrixCopy(msrc, ctypes.pointer(mdst))


lsci.setMatrixValue.argtypes = [
    ctypes.POINTER(MATRIX),
    ctypes.c_size_t,
    ctypes.c_size_t,
    ctypes.c_double,
]
lsci.setMatrixValue.restype = None


def set_matrix_value(mtx, irow, jcol, value):
    """Set one element of a C matrix.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        The matrix, modified in place.
    irow, jcol : int
        Row and column index (0-based, must be inside the matrix).
    value : float
        The new value.
    """
    lsci.setMatrixValue(mtx, irow, jcol, value)


lsci.getMatrixValue.argtypes = [ctypes.POINTER(MATRIX), ctypes.c_size_t, ctypes.c_size_t]
lsci.getMatrixValue.restype = ctypes.c_double


def get_matrix_value(mtx, irow, jcol):
    """Return one element of a C matrix.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        The matrix.
    irow, jcol : int
        Row and column index (0-based).

    Returns
    -------
    float
    """
    return lsci.getMatrixValue(mtx, irow, jcol)


def get_matrix_row(mtx, irow):
    """Return a row of a C matrix.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        The matrix.
    irow : int
        Row index (0-based).

    Returns
    -------
    list[float]
        The ``col`` values of the row.
    """
    row_lst = []
    for j in range(mtx[0].col):
        row_lst.append(get_matrix_value(mtx, irow, j))
    return row_lst


def get_matrix_column(mtx, jcol):
    """Return a column of a C matrix.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        The matrix.
    jcol : int
        Column index (0-based).

    Returns
    -------
    list[float]
        The ``row`` values of the column.
    """
    col_lst = []
    for i in range(mtx[0].row):
        col_lst.append(get_matrix_value(mtx, i, jcol))
    return col_lst


def matrix_to_list(mtx: Any) -> list[list[float]]:
    """Convert a C matrix into a list of rows.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX) or MATRIX
        The matrix.

    Returns
    -------
    list[list[float]]
        ``row`` lists of ``col`` values.
    """
    try:
        mtx_c = mtx[0]
    except TypeError:
        mtx_c = mtx
    ncols = mtx_c.col
    data = mtx_c.data
    # Slicing a ctypes double pointer converts a whole row in C.
    return [data[i][:ncols] for i in range(mtx_c.row)]


def matrix_to_numpy(mtx):
    """Copy a C matrix into a new numpy array (requires numpy).

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX) or MATRIX
        The matrix.

    Returns
    -------
    numpy.ndarray
        A ``float64`` array of shape ``(row, col)``.
    """
    import numpy as np  # pylint: disable=import-outside-toplevel

    try:
        mtx_c = mtx[0]
    except TypeError:
        mtx_c = mtx
    nrows, ncols = mtx_c.row, mtx_c.col
    out = np.empty((nrows, ncols), dtype=np.float64)
    if nrows == 0 or ncols == 0:
        return out
    nbytes = ncols * ctypes.sizeof(ctypes.c_double)
    base = out.ctypes.data
    data = mtx_c.data
    for i in range(nrows):
        ctypes.memmove(base + i * nbytes, data[i], nbytes)  # C-contiguous rows
    return out


def matrix_from_numpy(npm):
    """Create a C matrix from a numpy array.

    The data are copied row by row (no intermediate Python list). NaN/inf
    values are kept; the algorithms convert them to the missing value.

    Parameters
    ----------
    npm : numpy.ndarray
        A 2D array (a 1D array becomes a single column), any dtype convertible
        to ``float64``; non-contiguous arrays are accepted.

    Returns
    -------
    ctypes.POINTER(MATRIX)
        A new matrix; release it with :func:`del_matrix`.

    Raises
    ------
    ValueError
        If the array has more than 2 dimensions.
    """
    import numpy as np  # pylint: disable=import-outside-toplevel

    arr = np.ascontiguousarray(npm, dtype=np.float64)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    elif arr.ndim != 2:
        raise ValueError(f"matrix_from_numpy expects a 1D or 2D array, got {arr.ndim}D")
    nrows, ncols = arr.shape
    mtx = ctypes.POINTER(MATRIX)()
    lsci.NewMatrix(ctypes.pointer(mtx), nrows, ncols)
    if nrows == 0 or ncols == 0:
        return mtx
    nbytes = ncols * ctypes.sizeof(ctypes.c_double)
    base = arr.ctypes.data
    data = mtx[0].data
    for i in range(nrows):
        ctypes.memmove(data[i], base + i * nbytes, nbytes)  # C-contiguous rows
    return mtx


def set_missing_matrix_value(mtx, row_id, col_id):
    """Mark one element of a C matrix as missing.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        The matrix, modified in place.
    row_id, col_id : int
        Row and column index.
    """
    mtx[0].data[row_id][col_id] = misc.missing_value()


lsci.MatrixAppendRow.argtypes = [ctypes.POINTER(MATRIX), ctypes.POINTER(vector.DVECTOR)]
lsci.MatrixAppendRow.restype = None


def matrix_append_row(mtx, row_vect):
    """Append a row at the bottom of a C matrix.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        The matrix, modified in place.
    row_vect : list[float], DVector or POINTER(DVECTOR)
        The new row; its size must equal the number of columns (any size for an
        empty matrix).
    """
    alloc = False
    if isinstance(row_vect, vector.DVector):
        row_ptr = row_vect.dvect
    elif isinstance(row_vect, list):
        row_ptr = vector.new_dvector(row_vect)
        alloc = True
    else:
        row_ptr = row_vect
    res = lsci.MatrixAppendRow(mtx, row_ptr)
    if alloc:
        vector.del_dvector(row_ptr)
    return res


lsci.MatrixAppendCol.argtypes = [ctypes.POINTER(MATRIX), ctypes.POINTER(vector.DVECTOR)]
lsci.MatrixAppendCol.restype = None


def matrix_append_col(mtx, col_vect):
    """Append a column on the right of a C matrix.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        The matrix, modified in place.
    col_vect : list[float], DVector or POINTER(DVECTOR)
        The new column; its size must equal the number of rows (any size for an
        empty matrix).
    """
    alloc = False
    if isinstance(col_vect, vector.DVector):
        col_ptr = col_vect.dvect
    elif isinstance(col_vect, list):
        col_ptr = vector.new_dvector(col_vect)
        alloc = True
    else:
        col_ptr = col_vect
    res = lsci.MatrixAppendCol(mtx, col_ptr)
    if alloc:
        vector.del_dvector(col_ptr)
    return res


lsci.MatrixDVectorDotProduct.argtypes = [
    ctypes.POINTER(MATRIX),
    ctypes.POINTER(vector.DVECTOR),
    ctypes.POINTER(vector.DVECTOR),
]
lsci.MatrixDVectorDotProduct.restype = None


def matrix_dvector_dot_product(mtx, dvect, res):
    """Matrix-vector product ``res += mtx @ dvect``.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        Matrix, ``n x m``.
    dvect : POINTER(DVECTOR)
        Vector of size ``m``.
    res : POINTER(DVECTOR)
        Output vector of size ``n``. The product is **added** to it: pass a
        vector of zeros, e.g. ``vector.new_dvector([0.0] * n)``.

    Examples
    --------
    >>> import libscientific.vector as vector
    >>> m = new_matrix([[1.0, 2.0], [3.0, 4.0]])
    >>> v, r = vector.new_dvector([1.0, 1.0]), vector.new_dvector([0.0, 0.0])
    >>> matrix_dvector_dot_product(m, v, r)
    >>> vector.dvector_tolist(r)
    [3.0, 7.0]
    """
    return lsci.MatrixDVectorDotProduct(mtx, dvect, res)


lsci.MT_MatrixDVectorDotProduct.argtypes = [
    ctypes.POINTER(MATRIX),
    ctypes.POINTER(vector.DVECTOR),
    ctypes.POINTER(vector.DVECTOR),
]
lsci.MT_MatrixDVectorDotProduct.restype = None


def mt_matrix_dvector_dot_product(mtx, dvect, res):
    """Multi-threaded version of :func:`matrix_dvector_dot_product`.

    Threads are used only for large matrices (at least 65536 elements). Pass
    ``res`` filled with zeros: depending on the size the product is either
    written or added to it.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        Matrix, ``n x m``.
    dvect : POINTER(DVECTOR)
        Vector of size ``m``.
    res : POINTER(DVECTOR)
        Output vector of size ``n``, filled with zeros.
    """
    return lsci.MT_MatrixDVectorDotProduct(mtx, dvect, res)


lsci.DVectorMatrixDotProduct.argtypes = [
    ctypes.POINTER(MATRIX),
    ctypes.POINTER(vector.DVECTOR),
    ctypes.POINTER(vector.DVECTOR),
]
lsci.DVectorMatrixDotProduct.restype = None


def dvector_matrix_dot_product(mtx, dvect, res):
    """Vector-matrix product ``res += dvect @ mtx``.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        Matrix, ``n x m``.
    dvect : POINTER(DVECTOR)
        Vector of size ``n``.
    res : POINTER(DVECTOR)
        Output vector of size ``m``. The product is **added** to it: pass a
        vector of zeros.

    Examples
    --------
    >>> import libscientific.vector as vector
    >>> m = new_matrix([[1.0, 2.0], [3.0, 4.0]])
    >>> v, r = vector.new_dvector([1.0, 1.0]), vector.new_dvector([0.0, 0.0])
    >>> dvector_matrix_dot_product(m, v, r)
    >>> vector.dvector_tolist(r)
    [4.0, 6.0]
    """
    return lsci.DVectorMatrixDotProduct(mtx, dvect, res)


lsci.MT_DVectorMatrixDotProduct.argtypes = [
    ctypes.POINTER(MATRIX),
    ctypes.POINTER(vector.DVECTOR),
    ctypes.POINTER(vector.DVECTOR),
]
lsci.MT_DVectorMatrixDotProduct.restype = None


def mt_dvector_matrix_dot_product(mtx, dvect, res):
    """Multi-threaded version of :func:`dvector_matrix_dot_product`.

    Threads are used only for large matrices. Pass ``res`` filled with zeros.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        Matrix, ``n x m``.
    dvect : POINTER(DVECTOR)
        Vector of size ``n``.
    res : POINTER(DVECTOR)
        Output vector of size ``m``, filled with zeros.
    """
    return lsci.MT_DVectorMatrixDotProduct(mtx, dvect, res)


lsci.DVectorTrasposedDVectorDotProduct.argtypes = [
    ctypes.POINTER(vector.DVECTOR),
    ctypes.POINTER(vector.DVECTOR),
    ctypes.POINTER(MATRIX),
]
lsci.DVectorTrasposedDVectorDotProduct.restype = None


def dvector_transposed_dvector_dot_product(dvect1, dvect2, res):
    """Outer product of two vectors: ``res[i][j] = dvect1[i] * dvect2[j]``.

    Parameters
    ----------
    dvect1 : POINTER(DVECTOR)
        Vector of size ``n``.
    dvect2 : POINTER(DVECTOR)
        Vector of size ``m``.
    res : ctypes.POINTER(MATRIX)
        Output matrix, already allocated with shape ``n x m`` (see
        :func:`resize_matrix`).
    """
    return lsci.DVectorTrasposedDVectorDotProduct(dvect1, dvect2, res)


lsci.DVectorTransposedMatrixDivision.argtypes = [
    ctypes.POINTER(vector.DVECTOR),
    ctypes.POINTER(MATRIX),
    ctypes.POINTER(vector.DVECTOR),
]
lsci.DVectorTransposedMatrixDivision.restype = None


def dvector_transposed_matrix_division(dvect, mtx, res):
    """Solve ``res = inverse(mtx.T) @ dvect``.

    Parameters
    ----------
    dvect : POINTER(DVECTOR)
        Vector of size ``m``.
    mtx : ctypes.POINTER(MATRIX)
        Square matrix ``m x m`` (its transpose must be invertible).
    res : POINTER(DVECTOR)
        Output vector; resized to ``m``.
    """
    return lsci.DVectorTransposedMatrixDivision(dvect, mtx, res)


lsci.MatrixDotProduct.argtypes = [
    ctypes.POINTER(MATRIX),
    ctypes.POINTER(MATRIX),
    ctypes.POINTER(MATRIX),
]
lsci.MatrixDotProduct.restype = None


def matrix_dot_product(mtx_t, mtx, res):
    """Matrix product ``res += mtx_t @ mtx``.

    Parameters
    ----------
    mtx_t : ctypes.POINTER(MATRIX)
        Left matrix, ``n x k``.
    mtx : ctypes.POINTER(MATRIX)
        Right matrix, ``k x m``.
    res : ctypes.POINTER(MATRIX)
        Output matrix. An empty matrix (from :func:`init_matrix`) is allocated
        with shape ``n x m``; an already allocated ``n x m`` matrix gets the
        product **added** to its values.

    Examples
    --------
    >>> a, b = new_matrix([[1.0, 2.0], [3.0, 4.0]]), new_matrix([[1.0], [1.0]])
    >>> r = init_matrix()
    >>> matrix_dot_product(a, b, r)
    >>> matrix_to_list(r)
    [[3.0], [7.0]]
    """
    return lsci.MatrixDotProduct(mtx_t, mtx, res)


lsci.RowColOuterProduct.argtypes = [
    ctypes.POINTER(vector.DVECTOR),
    ctypes.POINTER(vector.DVECTOR),
    ctypes.POINTER(MATRIX),
]
lsci.RowColOuterProduct.restype = None


def row_col_outer_product(dvect_a, dvect_b, res):
    """Outer product of two vectors: ``res[i][j] = dvect_a[i] * dvect_b[j]``.

    Parameters
    ----------
    dvect_a : POINTER(DVECTOR)
        Vector of size ``n``.
    dvect_b : POINTER(DVECTOR)
        Vector of size ``m``.
    res : ctypes.POINTER(MATRIX)
        Output matrix, already allocated with shape ``n x m``.
    """
    return lsci.RowColOuterProduct(dvect_a, dvect_b, res)


lsci.MatrixTranspose.argtypes = [ctypes.POINTER(MATRIX), ctypes.POINTER(MATRIX)]
lsci.MatrixTranspose.restype = None


def matrix_transpose(mtx, mtx_t):
    """Transpose a matrix: ``mtx_t = mtx.T``.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        Matrix, ``n x m``.
    mtx_t : ctypes.POINTER(MATRIX)
        Output matrix; allocated or resized to ``m x n``.

    Examples
    --------
    >>> t = init_matrix()
    >>> matrix_transpose(new_matrix([[1.0, 2.0, 3.0]]), t)
    >>> matrix_to_list(t)
    [[1.0], [2.0], [3.0]]
    """
    return lsci.MatrixTranspose(mtx, mtx_t)


lsci.MatrixInversion.argtypes = [ctypes.POINTER(MATRIX), ctypes.POINTER(MATRIX)]
lsci.MatrixInversion.restype = None


def matrix_inversion(mtx, mtx_inv):
    """Invert a square matrix.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        Square, invertible matrix ``n x n``. For singular or rectangular
        matrices use the pseudo-inverse of the C API (``MatrixPseudoinversion``).
    mtx_inv : ctypes.POINTER(MATRIX)
        Output matrix; resized to ``n x n``.

    Examples
    --------
    >>> inv = init_matrix()
    >>> matrix_inversion(new_matrix([[2.0, 0.0], [0.0, 4.0]]), inv)
    >>> matrix_to_list(inv)
    [[0.5, 0.0], [0.0, 0.25]]
    """
    return lsci.MatrixInversion(mtx, mtx_inv)


lsci.EVectEval.argtypes = [
    ctypes.POINTER(MATRIX),
    ctypes.POINTER(vector.DVECTOR),
    ctypes.POINTER(MATRIX),
]
lsci.EVectEval.restype = None


def evect_eval(mtx, evals, evects):
    """Eigenvalues and eigenvectors of a square matrix (LAPACK).

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        Square matrix ``n x n``.
    evals : POINTER(DVECTOR)
        Output vector; receives the ``n`` eigenvalues, in the order returned by
    LAPACK (not sorted).
    evects : ctypes.POINTER(MATRIX)
        Output matrix; receives the eigenvectors as **columns**, ``n x n``:
        column ``k`` belongs to eigenvalue ``k``.
    """
    return lsci.EVectEval(mtx, evals, evects)


lsci.SVD.argtypes = [
    ctypes.POINTER(MATRIX),
    ctypes.POINTER(MATRIX),
    ctypes.POINTER(MATRIX),
    ctypes.POINTER(MATRIX),
]
lsci.SVD.restype = None


def svd(mtx, u_mtx, s_mtx, vt_mtx):
    """Thin singular value decomposition ``mtx = U @ S @ VT``.

    Same result as :func:`svd_lapack` (both use LAPACK).

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        Matrix, ``n x m``; let ``k = min(n, m)``.
    u_mtx : ctypes.POINTER(MATRIX)
        Output: left singular vectors, ``n x k``.
    s_mtx : ctypes.POINTER(MATRIX)
        Output: diagonal matrix of the singular values (largest first), ``k x k``.
    vt_mtx : ctypes.POINTER(MATRIX)
        Output: right singular vectors, transposed, ``k x m``.
    """
    lsci.SVD(mtx, u_mtx, s_mtx, vt_mtx)


lsci.SVDlapack.argtypes = [
    ctypes.POINTER(MATRIX),
    ctypes.POINTER(MATRIX),
    ctypes.POINTER(MATRIX),
    ctypes.POINTER(MATRIX),
]
lsci.SVDlapack.restype = None


def svd_lapack(mtx, u_mtx, s_mtx, vt_mtx):
    """Thin singular value decomposition ``mtx = U @ S @ VT`` with LAPACK ``dgesdd``.

    Parameters
    ----------
    mtx : ctypes.POINTER(MATRIX)
        Matrix, ``n x m``; let ``k = min(n, m)``.
    u_mtx : ctypes.POINTER(MATRIX)
        Output: left singular vectors, ``n x k``.
    s_mtx : ctypes.POINTER(MATRIX)
        Output: diagonal matrix of the singular values (largest first), ``k x k``.
    vt_mtx : ctypes.POINTER(MATRIX)
        Output: right singular vectors, transposed, ``k x m``.
    """
    lsci.SVDlapack(mtx, u_mtx, s_mtx, vt_mtx)


class Matrix:
    """A matrix stored in the C library, with a Pythonic interface.

    Parameters
    ----------
    mx_ : list[list[float]], numpy.ndarray, Matrix, C matrix pointer or None
        Initial content: a list of rows, a 2D numpy array, another
        :class:`Matrix` or C matrix (copied), or ``None`` for an empty matrix.

    Attributes
    ----------
    mtx : ctypes.POINTER(MATRIX)
        The underlying C matrix, owned by this object.

    Notes
    -----
    Indexing: ``m[i, j]`` reads/writes one element and ``m[i]`` returns row
    ``i`` as a list. ``len(m)`` is the number of rows and iterating yields the
    rows as lists.

    Examples
    --------
    >>> m = Matrix([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
    >>> m.shape, m[1], m[2, 0]
    ((3, 2), [3.0, 4.0], 5.0)
    >>> m[0, 0] = 10.0
    >>> m.tolist()[0]
    [10.0, 2.0]
    >>> m.transpose()
    >>> m.shape
    (2, 3)
    """

    def __init__(self, mx_=None):
        """Create the matrix (see the class documentation for the accepted inputs)."""
        if mx_ is None:
            self.mtx = init_matrix()
        elif isinstance(mx_, Matrix):
            self.mtx = init_matrix()
            matrix_copy(mx_.mtx, self.mtx)
        elif isinstance(mx_, ctypes.POINTER(MATRIX)):
            self.mtx = init_matrix()
            matrix_copy(mx_, self.mtx)
        elif "numpy" in str(type(mx_)):
            self.mtx = matrix_from_numpy(mx_)
        else:
            self.mtx = new_matrix(mx_)

    def __del__(self):
        """Release the C matrix."""
        if getattr(self, "mtx", None) is not None:
            del_matrix(self.mtx)
            del self.mtx
            self.mtx = None

    def __repr__(self):
        return f"Matrix({self.tolist()})"

    def __str__(self):
        return f"Matrix({self.tolist()})"

    def __len__(self):
        """Return the number of rows."""
        return self.nrow()

    def __iter__(self):
        """Iterate over the rows, as lists."""
        return iter(self.tolist())

    @property
    def shape(self):
        """``(n_rows, n_columns)`` of the matrix."""
        return (self.nrow(), self.ncol())

    def __getitem__(self, keys):
        """Return an element (``m[i, j]``) or a row (``m[i]``).

        Parameters
        ----------
        keys : tuple[int, int] or int
            ``(row, column)`` for an element, or a row index.

        Returns
        -------
        float or list[float]

        Raises
        ------
        TypeError
            For any other kind of index.
        """
        if isinstance(keys, tuple):
            i, j = self._check_indexes(*keys)
            return self.data_ptr()[i][j]
        if isinstance(keys, int):
            i, _ = self._check_indexes(keys, 0 if self.ncol() > 0 else None)
            return get_matrix_row(self.mtx, i)
        raise TypeError(f"Invalid index type: {type(keys)}")

    def _check_indexes(self, i: int, j: int | None) -> tuple[int, int]:
        """Validate (row, column) indexes; negative values count from the end."""
        nrow, ncol = self.nrow(), self.ncol()
        ii = i + nrow if i < 0 else i
        if not 0 <= ii < nrow:
            raise IndexError(f"row index {i} out of range for {nrow} rows")
        if j is None:
            return ii, 0
        jj = j + ncol if j < 0 else j
        if not 0 <= jj < ncol:
            raise IndexError(f"column index {j} out of range for {ncol} columns")
        return ii, jj

    def __setitem__(self, keys, value):
        """Set an element: ``m[i, j] = value``.

        Parameters
        ----------
        keys : tuple[int, int]
            ``(row, column)``, inside the matrix.
        value : float
            The new value.
        """
        i, j = self._check_indexes(*keys)
        set_matrix_value(self.mtx, i, j, value)

    def nrow(self):
        """Return the number of rows.

        Returns
        -------
        int
        """
        return self.mtx[0].row

    def ncol(self):
        """Return the number of columns.

        Returns
        -------
        int
        """
        return self.mtx[0].col

    def data_ptr(self):
        """Return the C row pointers of the matrix (low-level access).

        Returns
        -------
        POINTER(POINTER(c_double))
            ``data_ptr()[i][j]`` is element ``(i, j)``; valid while the matrix exists.
        """
        return self.mtx[0].data

    def tolist(self):
        """Return the content as a list of rows.

        Returns
        -------
        list[list[float]]
        """
        return matrix_to_list(self.mtx)

    def fromlist(self, mtx_):
        """Replace the content with a list of rows.

        Parameters
        ----------
        mtx_ : list[list[float]]
            The new rows.
        """
        del_matrix(self.mtx)
        del self.mtx
        self.mtx = new_matrix(mtx_)

    def fromnumpy(self, npmx):
        """Replace the content with a numpy array.

        Parameters
        ----------
        npmx : numpy.ndarray
            A 2D array.
        """
        del_matrix(self.mtx)
        del self.mtx
        self.mtx = matrix_from_numpy(npmx)

    def appendrow(self, row_lst_):
        """Append a row at the bottom.

        Parameters
        ----------
        row_lst_ : list[float] or DVector
            The new row, with one value per column.
        """
        matrix_append_row(self.mtx, row_lst_)

    def appendcol(self, col_lst_):
        """Append a column on the right.

        Parameters
        ----------
        col_lst_ : list[float] or DVector
            The new column, with one value per row.
        """
        matrix_append_col(self.mtx, col_lst_)

    def transpose(self):
        """Transpose the matrix in place (``n x m`` becomes ``m x n``)."""
        mtx_t = init_matrix()
        resize_matrix(mtx_t, self.mtx[0].col, self.mtx[0].row)
        matrix_transpose(self.mtx, mtx_t)
        matrix_copy(mtx_t, self.mtx)
        del_matrix(mtx_t)
        del mtx_t

    def get_evect_evals(self):
        """Return the eigenvalues and eigenvectors of the (square) matrix.

        Returns
        -------
        eigenvalues : list[float]
            ``n`` values, in the order returned by LAPACK (not sorted).
        eigenvectors : list[list[float]]
            ``n x n`` matrix whose **columns** are the eigenvectors: column ``k``
            belongs to ``eigenvalues[k]``.

        Examples
        --------
        >>> evals, evects = Matrix([[2.0, 0.0], [0.0, 3.0]]).get_evect_evals()
        >>> sorted(round(v, 6) for v in evals)
        [2.0, 3.0]
        """
        evals_ = vector.init_dvector()
        evects_ = init_matrix()
        evect_eval(self.mtx, evals_, evects_)
        evals = vector.dvector_tolist(evals_)
        evects = matrix_to_list(evects_)
        del_matrix(evects_)
        vector.del_dvector(evals_)
        return evals, evects

    def svd(self):
        """Return the thin singular value decomposition ``M = U @ S @ VT``.

        Returns
        -------
        U : list[list[float]]
            Left singular vectors, ``n x k`` with ``k = min(n, m)``.
        S : list[list[float]]
            Diagonal ``k x k`` matrix of the singular values, largest first.
        VT : list[list[float]]
            Right singular vectors, transposed, ``k x m``.

        Examples
        --------
        >>> U, S, VT = Matrix([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]).svd()
        >>> len(U), len(U[0]), len(S), len(VT), len(VT[0])
        (3, 2, 2, 2, 2)
        >>> [round(S[i][i], 4) for i in range(2)]
        [9.5255, 0.5143]
        """
        mtx_u = init_matrix()
        mtx_s = init_matrix()
        mtx_vt = init_matrix()
        svd_lapack(self.mtx, mtx_u, mtx_s, mtx_vt)
        mtx_u_lst = matrix_to_list(mtx_u)
        mtx_s_lst = matrix_to_list(mtx_s)
        mtx_vt_lst = matrix_to_list(mtx_vt)
        del_matrix(mtx_u)
        del_matrix(mtx_s)
        del_matrix(mtx_vt)
        return mtx_u_lst, mtx_s_lst, mtx_vt_lst

    def debug(self):
        """Print the matrix to standard output."""
        print_matrix(self.mtx)
