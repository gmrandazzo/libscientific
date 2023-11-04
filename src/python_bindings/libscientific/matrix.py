"""matrix.py libscientific python binding

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
import csv
from libscientific.loadlibrary import load_libscientific_library
from libscientific import misc
from libscientific import vector

lsci = load_libscientific_library()

class MATRIX(ctypes.Structure):
    """
    matrix data structure
    """
    _fields_ = [
        ("data",    ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
        ("row",     ctypes.c_size_t),
        ("col",     ctypes.c_size_t)]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


lsci.initMatrix.argtypes = [ctypes.POINTER(ctypes.POINTER(MATRIX))]
lsci.initMatrix.restype = None


def init_matrix():
    """
    init_matrix: Allocate in memory an empty libscientific matrix
    """
    mtx = ctypes.POINTER(MATRIX)()
    # m = matrix()
    lsci.initMatrix(ctypes.pointer(mtx))
    return mtx


lsci.NewMatrix.argtypes = [ctypes.POINTER(ctypes.POINTER(MATRIX)),
                           ctypes.c_size_t,
                           ctypes.c_size_t]
lsci.NewMatrix.restype = None

def new_matrix_from_list(mx_lst):
    """
    Import a matrix from a list
    """
    nrows = None
    ncols = None
    try:
        nrows = len(mx_lst)
        try:
            ncols = len(mx_lst[0])
        except TypeError:
            ncols = 1
        except IndexError:
            ncols = 0
    except IndexError:
        nrows = 0
    mtx = ctypes.POINTER(MATRIX)()
    lsci.NewMatrix(ctypes.pointer(mtx),
                   nrows,
                   ncols)

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

def new_matrix_from_csv(fcsv):
    """
    Import a matrix directly from CSV file
    """
    mtx = init_matrix()
    with open(fcsv, "r", encoding='utf-8') as fcsv_:
        reader = csv.reader(fcsv_, delimiter=',', quotechar='"')
        for csvrow in reader:
            row_lst = vector.DVector(csvrow)
            matrix_append_row(mtx, row_lst)
            del row_lst
    return mtx

def new_matrix(mtx_input_):
    """
    new_matrix: Allocate in memory a libscientific matrix from a list of lists
    """
    if "str" in str(type(mtx_input_)):
        return new_matrix_from_csv(mtx_input_)
    return new_matrix_from_list(mtx_input_)


lsci.ResizeMatrix.argtypes = [ctypes.POINTER(MATRIX),
                              ctypes.c_size_t,
                              ctypes.c_size_t]
lsci.ResizeMatrix.restype = None


def resize_matrix(mtx, nrows, ncols):
    """
    resize_matrix: Resize an already allocated or reallocate a libscientific matrix
    """
    lsci.ResizeMatrix(mtx, nrows, ncols)


lsci.DelMatrix.argtypes = [ctypes.POINTER(ctypes.POINTER(MATRIX))]
lsci.DelMatrix.restype = None


def del_matrix(mtx):
    """
    del_matrix: Delete an allocated libscientific matrix
    """
    lsci.DelMatrix(ctypes.pointer(mtx))


lsci.MatrixCheck.argtypes = [ctypes.POINTER(MATRIX)]
lsci.MatrixCheck.restype = None


def matrix_check(mtx):
    """
    matrix_check: Find infinite and nan numbers and sobstitute with MISSING value.
    """
    lsci.MatrixCheck(mtx)


lsci.PrintMatrix.argtypes = [ctypes.POINTER(MATRIX)]
lsci.PrintMatrix.restype = None


def print_matrix(mtx):
    """
    print_matrix: Print to video a libscientific matrix
    """
    lsci.PrintMatrix(mtx)


lsci.ValInMatrix.argtypes = [ctypes.POINTER(MATRIX), ctypes.c_double]
lsci.ValInMatrix.restype = ctypes.c_int


def val_in_matrix(mtx, val):
    """
    val_in_matrix: Check if a libscientific matrix contains an exact value "val"
    and return 1 or 0 respectivelly for yes or no.
    """
    return lsci.ValInMatrix(mtx, val)


lsci.MatrixSet.argtypes = [ctypes.POINTER(MATRIX), ctypes.c_double]
lsci.MatrixSet.restype = None


def matrix_set(mtx, val):
    """
    matrix_set: Set all values of a libscientific matrix to "val"
    """
    lsci.MatrixSet(mtx, val)


lsci.MatrixCopy.argtypes = [ctypes.POINTER(MATRIX),
                            ctypes.POINTER(ctypes.POINTER(MATRIX))]
lsci.MatrixCopy.restype = None


def matrix_copy(msrc, mdst):
    """
    matrix_copy: Copy a libscientifi matrix to another allocated one
    """
    lsci.MatrixCopy(msrc, ctypes.pointer(mdst))


lsci.setMatrixValue.argtypes = [ctypes.POINTER(MATRIX),
                                ctypes.c_size_t,
                                ctypes.c_size_t,
                                ctypes.c_double]
lsci.setMatrixValue.restype = None


def set_matrix_value(mtx, irow, jcol, value):
    """
    set_matrix_value: Set/modify a value in the irow and jcol of a libscientific matrix
    """
    lsci.setMatrixValue(mtx, irow, jcol, value)


lsci.getMatrixValue.argtypes = [ctypes.POINTER(MATRIX),
                                ctypes.c_size_t,
                                ctypes.c_size_t]
lsci.getMatrixValue.restype = ctypes.c_double


def get_matrix_value(mtx, irow, jcol):
    """
    get_matrix_value: Get a value in the irow and jcol of a libscientific matrix
    """
    return lsci.getMatrixValue(mtx, irow, jcol)


def get_matrix_row(mtx, irow):
    """
    get_matrix_row: Python version of getting a row from a libscientific matrix
    """
    row_lst = []
    for j in range(mtx[0].col):
        row_lst.append(get_matrix_value(mtx, irow, j))
    return row_lst


def get_matrix_column(mtx, jcol):
    """
    getMatrixColumn: Python version of getting a column from a libscientific  matrix
    """
    col_lst = []
    for i in range(mtx[0].row):
        col_lst.append(get_matrix_value(mtx, i, jcol))
    return col_lst


def matrix_to_list(mtx):
    """
    matrix_to_list: Convert a libscientific matrix to list of list
    """
    mtx_list = []
    try:
        for i in range(mtx[0].row):
            row_lst = []
            for j in range(mtx[0].col):
                row_lst.append(mtx[0].data[i][j])
            mtx_list.append(row_lst)
    except TypeError:
        for i in range(mtx.row):
            row_lst = []
            for j in range(mtx.col):
                row_lst.append(mtx.data[i][j])
            mtx_list.append(row_lst)
    return mtx_list


def matrix_from_numpy(npm):
    """
    Convert a numpy matrix into a libscientific matrix
    """
    return new_matrix(npm.tolist())


def set_missing_matrix_value(mtx, row_id, col_id):
    """
    Set a value of a matrix into a missing value
    """
    mtx[0].data[row_id][col_id] = misc.missing_value()


lsci.MatrixAppendRow.argtypes = [ctypes.POINTER(MATRIX),
                                 ctypes.POINTER(vector.DVECTOR)]
lsci.MatrixAppendRow.restype = None

def matrix_append_row(mtx, row_vect):
    """
    Append to a libscientific matrix a row list
    """
    return lsci.MatrixAppendRow(mtx, row_vect.dvect)


lsci.MatrixAppendCol.argtypes = [ctypes.POINTER(MATRIX),
                                 ctypes.POINTER(vector.DVECTOR)]
lsci.MatrixAppendCol.restype = None

def matrix_append_col(mtx, col_vect):
    """
    Append to a libscientific matrix a column list
    """
    return lsci.MatrixAppendCol(mtx, col_vect.dvect)


lsci.MatrixDVectorDotProduct.argtypes = [ctypes.POINTER(MATRIX),
                                         ctypes.POINTER(vector.DVECTOR),
                                         ctypes.POINTER(vector.DVECTOR)]
lsci.MatrixDVectorDotProduct.restype = None

def matrix_dvector_dot_product(mtx, dvect, res):
    """
    Calculate a matrix-vector dot product
    """
    return lsci.MatrixDVectorDotProduct(mtx, dvect, res)


lsci.MT_MatrixDVectorDotProduct.argtypes = [ctypes.POINTER(MATRIX),
                                            ctypes.POINTER(vector.DVECTOR),
                                            ctypes.POINTER(vector.DVECTOR)]
lsci.MT_MatrixDVectorDotProduct.restype = None


def mt_matrix_dvector_dot_product(mtx, dvect, res):
    """
    Calculate a matrix-vector dot product in multi-threads
    """
    return lsci.MT_MatrixDVectorDotProduct(mtx, dvect, res)


lsci.DVectorMatrixDotProduct.argtypes = [ctypes.POINTER(MATRIX),
                                         ctypes.POINTER(vector.DVECTOR),
                                         ctypes.POINTER(vector.DVECTOR)]
lsci.DVectorMatrixDotProduct.restype = None

def dvector_matrix_dot_product(mtx, dvect, res):
    """
    Calculate a vector-matrix dot product
    """
    return lsci.DVectorMatrixDotProduct(mtx, dvect, res)



lsci.MT_DVectorMatrixDotProduct.argtypes = [ctypes.POINTER(MATRIX),
                                            ctypes.POINTER(vector.DVECTOR),
                                            ctypes.POINTER(vector.DVECTOR)]
lsci.MT_DVectorMatrixDotProduct.restype = None

def mt_dvector_matrix_dot_product(mtx, dvect, res):
    """
    Calculate a vector-matrix dot product in multi-threads
    """
    return lsci.MT_DVectorMatrixDotProduct(mtx, dvect, res)


lsci.DVectorTrasposedDVectorDotProduct.argtypes = [ctypes.POINTER(vector.DVECTOR),
                                                   ctypes.POINTER(vector.DVECTOR),
                                                   ctypes.POINTER(MATRIX)]
lsci.DVectorTrasposedDVectorDotProduct.restype = None

def dvector_transposed_dvector_dot_product(mtx, dvect, res):
    """
    Calculate a dvector transposed dvector dot product
    """
    return lsci.DVectorTrasposedDVectorDotProduct(mtx, dvect, res)


lsci.DVectorTransposedMatrixDivision.argtypes = [ctypes.POINTER(vector.DVECTOR),
                                                 ctypes.POINTER(MATRIX),
                                                 ctypes.POINTER(vector.DVECTOR)]
lsci.DVectorTransposedMatrixDivision.restype = None

def dvector_transposed_matrix_division(dvect, mtx, res):
    """
    Calculate a dvector transposed-matrix division
    """
    return lsci.DVectorTransposedMatrixDivision(dvect, mtx, res)


lsci.MatrixDotProduct.argtypes = [ctypes.POINTER(MATRIX),
                                  ctypes.POINTER(MATRIX),
                                  ctypes.POINTER(MATRIX)]
lsci.MatrixDotProduct.restype = None

def matrix_dot_product(mtx_t, mtx, res):
    """
    Calculate the matrix-matrix product
    """
    return lsci.MatrixDotProduct(mtx_t, mtx, res)


lsci.RowColOuterProduct.argtypes = [ctypes.POINTER(vector.DVECTOR),
                                    ctypes.POINTER(vector.DVECTOR),
                                    ctypes.POINTER(MATRIX)]
lsci.RowColOuterProduct.restype = None

def row_col_outer_product(dvect_a, dvect_b, res):
    """
    Calculate the row-column outer product of two dvectors
    """
    return lsci.RowColOuterProduct(dvect_a, dvect_b, res)


lsci.MatrixTranspose.argtypes = [ctypes.POINTER(MATRIX),
                                 ctypes.POINTER(MATRIX)]
lsci.MatrixTranspose.restype = None

def matrix_transpose(mtx, mtx_t):
    """
    Transpose a matrix
    """
    return lsci.MatrixTranspose(mtx, mtx_t)



lsci.MatrixInversion.argtypes = [ctypes.POINTER(MATRIX),
                                 ctypes.POINTER(MATRIX)]
lsci.MatrixInversion.restype = None

def matrix_inversion(mtx, mtx_inv):
    """
    Invert a matrix
    """
    return lsci.MatrixInversion(mtx, mtx_inv)


lsci.EVectEval.argtypes = [ctypes.POINTER(MATRIX),
                           ctypes.POINTER(vector.DVECTOR),
                           ctypes.POINTER(MATRIX)]
lsci.EVectEval.restype = None

def evect_eval(mtx, mtx_e_vect, mtx_e_val):
    """
    Eigenvectors and Eigenvalues using lapack
    """
    return lsci.EVectEval(mtx, mtx_e_vect, mtx_e_val)


lsci.SVD.argtypes = [ctypes.POINTER(MATRIX),
                     ctypes.POINTER(MATRIX),
                     ctypes.POINTER(MATRIX),
                     ctypes.POINTER(MATRIX)]
lsci.SVD.restype = None

def svd(mtx, u_mtx, s_mtx, vt_mtx):
    """
    Calculate the SVD of a matrix mtx using the eigenvector/values lapack method
    """
    lsci.SVD(mtx, u_mtx, s_mtx, vt_mtx)


lsci.SVDlapack.argtypes = [ctypes.POINTER(MATRIX),
                           ctypes.POINTER(MATRIX),
                           ctypes.POINTER(MATRIX),
                           ctypes.POINTER(MATRIX)]
lsci.SVDlapack.restype = None

def svd_lapack(mtx, u_mtx, s_mtx, vt_mtx):
    """
    Calculate the SVD of a matrix mtx using the SVD lapack method
    """
    lsci.SVDlapack(mtx, u_mtx, s_mtx, vt_mtx)


class Matrix():
    """
    A class for working with libscientific matrices.

    This class provides methods to perform various operations on matrices,
    including creating, accessing, modifying, and transforming matrices using
    the libscientific library.

    Attributes:
        mtx (CDataType): The underlying libscientific matrix data.

    Methods:
        __init__(self, mx_=None)
        __del__(self)
        __getitem__(self, keys)
        __setitem__(self, keys, value)
        nrow(self)
        ncol(self)
        data_ptr(self)
        tolist(self)
        fromlist(self, mtx_)
        fromnumpy(self, npmx)
        appendrow(self, row_lst_)
        appendcol(self, col_lst_)
        transpose(self)
        get_evect_evals(self)
        svd(self)
    """
    def __init__(self, mx_=None):
        """
        Initialize a Matrix instance.

        Args:
            mx_ (list[list]): A list of lists representing the initial matrix data. Default is None.
        """
        if mx_ is None:
            self.mtx = init_matrix()
        else:
            self.mtx = new_matrix(mx_)

    def __del__(self):
        """
        Clean up resources associated with the Matrix instance.
        """
        del_matrix(self.mtx)
        del self.mtx
        self.mtx = None

    def __getitem__(self, keys):
        """
        Get a value from the matrix using row and column indices.

        Args:
            keys (tuple): A tuple containing two integers representing the row and column indices.

        Returns:
            float: The value at the specified row and column in the matrix.
        """
        i, j = keys
        return self.data_ptr()[i][j]

    def __setitem__(self, keys, value):
        """
        Set a value in the matrix using row and column indices.

        Args:
            keys (tuple): A tuple containing two integers representing the row and column indices.
            value (float): The value to be set at the specified row and column in the matrix.
        """
        i, j = keys
        set_matrix_value(self.mtx, i, j, value)

    def nrow(self):
        """
        Get the number of rows in the matrix.

        Returns:
            int: The number of rows in the matrix.
        """
        return self.mtx[0].row

    def ncol(self):
        """
        Get the number of columns in the matrix.

        Returns:
            int: The number of columns in the matrix.
        """
        return self.mtx[0].col

    def data_ptr(self):
        """
        Get the data pointer of the matrix.

        Returns:
            CDataType: The data pointer of the matrix.
        """
        return self.mtx[0].data

    def tolist(self):
        """
        Convert the matrix to a list of lists.

        Returns:
            list[list]: A list of lists representing the matrix.
        """
        return matrix_to_list(self.mtx)

    def fromlist(self, mtx_):
        """
        Set the matrix data from a list of lists.

        Args:
            mtx_ (list[list]): A list of lists representing the new matrix data.
        """
        del_matrix(self.mtx)
        del self.mtx
        self.mtx = new_matrix(mtx_)

    def fromnumpy(self, npmx):
        """
        Set the matrix data from a NumPy array.

        Args:
            npmx (numpy.ndarray): A NumPy array representing the new matrix data.
        """
        del_matrix(self.mtx)
        del self.mtx
        self.mtx = matrix_from_numpy(npmx)

    def appendrow(self, row_lst_):
        """
        Append a row to the matrix.

        Args:
            row_lst_ (list): A list representing the row to be appended.
        """
        row_lst = vector.DVector(row_lst_)
        matrix_append_row(self.mtx, row_lst)
        del row_lst

    def appendcol(self, col_lst_):
        """
        Append a column vector to the matrix.

        Args:
            col_lst_ (list): A list representing the column vector to be appended.
        """
        col_lst = vector.DVector(col_lst_)
        matrix_append_col(self.mtx, col_lst)
        del col_lst

    def transpose(self):
        """
        Transpose the matrix in-place.
        """
        mtx_t = init_matrix()
        resize_matrix(mtx_t, self.mtx[0].col, self.mtx[0].row)
        matrix_transpose(self.mtx, mtx_t)
        matrix_copy(mtx_t, self.mtx)
        del_matrix(mtx_t)
        del mtx_t

    def get_evect_evals(self):
        """
        Get the eigenvectors and eigenvalues of the matrix.

        Returns:
            tuple: A tuple containing two lists - eigenvectors and eigenvalues.
        """
        mtx_e_vect_ = vector.init_dvector()
        mtx_e_vals_ = init_matrix()
        evect_eval(self.mtx,  mtx_e_vect_, mtx_e_vals_)
        mtx_e_vect = vector.dvector_tolist(mtx_e_vect_)
        mtx_e_vals = matrix_to_list(mtx_e_vals_)
        del_matrix(mtx_e_vals_)
        vector.del_dvector(mtx_e_vect_)
        return mtx_e_vect, mtx_e_vals

    def svd(self):
        """
        Perform Singular Value Decomposition (SVD) on the matrix.

        Returns:
            tuple: A tuple containing three lists - U, S, and VT matrices.
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
        """
        Debug the matrix content
        """
        print_matrix(self.mtx)
