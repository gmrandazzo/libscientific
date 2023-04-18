"""
matrix libscientific python binding

Copyright (C) <2019>  Giuseppe Marco Randazzo

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


def new_matrix(mtx_input_):
    """
    new_matrix: Allocate in memory a libscientific matrix from a list of lists
    """
    mtx_input = None
    if "numpy" in str(type(mtx_input_)):
        mtx_input = mtx_input_.tolist()
    else:
        mtx_input = mtx_input_
    nrows = None
    ncols = None
    try:
        nrows = len(mtx_input)
        try:
            ncols = len(mtx_input[0])
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
                val = float(mtx_input[i][j])
            except TypeError:
                val = float(mtx_input[i])
            except ValueError:
                val = None

            if val is None:
                set_missing_matrix_value(mtx, i, j)
            else:
                lsci.setMatrixValue(mtx, i, j, val)
    return mtx


lsci.ResizeMatrix.argtypes = [ctypes.POINTER(MATRIX),
                              ctypes.c_size_t,
                              ctypes.c_size_t]
lsci.ResizeMatrix.restype = None


def resize_matrix(mtx, nrows, ncols):
    """
    resize_matrix: Resize an already allocated or reallocate a libscientific
                  matrix
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
    matrix_check: Find infinite and nan numbers and sobstitute
                 with MISSING value.
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
    set_matrix_value: Set/modify a value in the irow and jcol of a libscientific
                    matrix
    """
    lsci.setMatrixValue(mtx, irow, jcol, value)


lsci.getMatrixValue.argtypes = [ctypes.POINTER(MATRIX),
                                ctypes.c_size_t,
                                ctypes.c_size_t]
lsci.getMatrixValue.restype = ctypes.c_double


def get_matrix_value(mtx, irow, jcol):
    """
    get_matrix_value: Get a value in the irow and jcol of a libscientific
                    matrix
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
    getMatrixColumn: Python version of getting a column from a libscientific
                     matrix
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
    /* Description:
    * matrix - row double vector product: the result is a row double vector
    * i.e.: X(10x5) * d(5x1) = r(10x1)
    */
    void matrix_dvector_dot_product(matrix *mtx, dvector *v, dvector *r);
    """
    return lsci.MatrixDVectorDotProduct(mtx, dvect, res)


lsci.MT_MatrixDVectorDotProduct.argtypes = [ctypes.POINTER(MATRIX),
                                            ctypes.POINTER(vector.DVECTOR),
                                            ctypes.POINTER(vector.DVECTOR)]
lsci.MT_MatrixDVectorDotProduct.restype = None


def mt_matrix_dvector_dot_product(mtx, dvect, res):
    """
    /* Multithread version of matrix_dvector_dot_product */
    void MT_MatrixDVectorDotProduct(matrix *mx, dvector *v, dvector *p);
    """
    return lsci.MT_MatrixDVectorDotProduct(mtx, dvect, res)


lsci.DVectorMatrixDotProduct.argtypes = [ctypes.POINTER(MATRIX),
                                         ctypes.POINTER(vector.DVECTOR),
                                         ctypes.POINTER(vector.DVECTOR)]
lsci.DVectorMatrixDotProduct.restype = None

def dvector_matrix_dot_product(mtx, dvect, res):
    """
    /* Description:
    * column double vector - matrix product: the result is a column double vector
    * i.e.: d(1x5) * X(5x10) = r(1x10)
    */
    void DVectorMatrixDotProduct(matrix *mx, dvector *v, dvector *p);
    """
    return lsci.DVectorMatrixDotProduct(mtx, dvect, res)



lsci.MT_DVectorMatrixDotProduct.argtypes = [ctypes.POINTER(MATRIX),
                                            ctypes.POINTER(vector.DVECTOR),
                                            ctypes.POINTER(vector.DVECTOR)]
lsci.MT_DVectorMatrixDotProduct.restype = None

def mt_dvector_matrix_dot_product(mtx, dvect, res):
    """
    /* Multithread version of DVectorMatrixDotProduct */
    void MT_DVectorMatrixDotProduct(matrix *mx, dvector *v, dvector *p);
    """
    return lsci.MT_DVectorMatrixDotProduct(mtx, dvect, res)


lsci.DVectorTrasposedDVectorDotProduct.argtypes = [ctypes.POINTER(vector.DVECTOR),
                                                   ctypes.POINTER(vector.DVECTOR),
                                                   ctypes.POINTER(MATRIX)]
lsci.DVectorTrasposedDVectorDotProduct.restype = None

def dvector_transposed_dvector_dot_product(mtx, dvect, res):
    """
    /* Description:
    * transposed double vector - double vecrtor product: the result is a matrix
    * i.e.: X  = d'd
    * i.e.: d'(5x1) * d(1x5) = X(1x10)
    */
    void DVectorTrasposedDVectorDotProduct(dvector *v1, dvector *v2, matrix *m);
    """
    return lsci.DVectorTrasposedDVectorDotProduct(mtx, dvect, res)


lsci.DVectorTransposedMatrixDivision.argtypes = [ctypes.POINTER(vector.DVECTOR),
                                                 ctypes.POINTER(MATRIX),
                                                 ctypes.POINTER(vector.DVECTOR)]
lsci.DVectorTransposedMatrixDivision.restype = None

def dvector_transposed_matrix_division(dvect, mtx, res):
    """
    /* r = v/mx  = (inv(mx^T)*v^T)^T*/
    void DVectorTransposedMatrixDivision(dvector *v, matrix *mx, dvector *r);
    """
    return lsci.DVectorTransposedMatrixDivision(dvect, mtx, res)


lsci.MatrixDotProduct.argtypes = [ctypes.POINTER(MATRIX),
                                  ctypes.POINTER(MATRIX),
                                  ctypes.POINTER(MATRIX)]
lsci.MatrixDotProduct.restype = None

def matrix_dot_product(mtx_t, mtx, res):
    """
    /*
    * Description:
    * Calculate the matrix matrix product
    */
    void MatrixDotProduct(matrix *m_t, matrix *mtx, matrix *r);
    """
    return lsci.MatrixDotProduct(mtx_t, mtx, res)


lsci.RowColOuterProduct.argtypes = [ctypes.POINTER(vector.DVECTOR),
                                    ctypes.POINTER(vector.DVECTOR),
                                    ctypes.POINTER(MATRIX)]
lsci.RowColOuterProduct.restype = None

def row_col_outer_product(dvect_a, dvect_b, res):
    """
    /*
    * Description:
    * Calculate the matrix matrix product
    */
    void RowColOuterProduct(dvector *a, dvector *b, matrix *m);
    """
    return lsci.RowColOuterProduct(dvect_a, dvect_b, res)


lsci.MatrixTranspose.argtypes = [ctypes.POINTER(MATRIX),
                                 ctypes.POINTER(MATRIX)]
lsci.MatrixTranspose.restype = None

def matrix_transpose(mtx, mtx_t):
    """
    /*
    * Description:
    * Generate a transpose matrix of m
    */
    void MatrixTranspose(matrix *mtx, matrix *r);
    """
    return lsci.MatrixTranspose(mtx, mtx_t)



lsci.MatrixInversion.argtypes = [ctypes.POINTER(MATRIX),
                                 ctypes.POINTER(MATRIX)]
lsci.MatrixInversion.restype = None

def matrix_inversion(mtx, mtx_inv):
    """
    /*
    * Description:
    * Matrix inversion using the Gauss-Jordan algorithm
    */
    void MatrixInversion(matrix *mtx, matrix **m_inv);
    """
    return lsci.MatrixInversion(mtx, mtx_inv)


lsci.EVectEval.argtypes = [ctypes.POINTER(MATRIX),
                           ctypes.POINTER(vector.DVECTOR),
                           ctypes.POINTER(MATRIX)]
lsci.EVectEval.restype = None

def evect_eval(mtx, mtx_e_vect, mtx_e_val):
    """
    Eigenvectors and Eigenvalues using a lapack method
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
    Translate a list of list into a libscientific matrix
    """
    def __init__(self, mx_):
        self.mtx = new_matrix(mx_)

    def __del__(self):
        del_matrix(self.mtx)
        del self.mtx
        self.mtx = None

    def __getitem__(self, keys):
        i, j = keys
        return self.data_ptr()[i][j]

    def __setitem__(self, keys, value):
        i, j = keys
        set_matrix_value(self.mtx, i, j, value)

    def nrow(self):
        """
        Return the number of rows of the matrix
        """
        return self.mtx[0].row

    def ncol(self):
        """
        Return the number of columns of the matrix
        """
        return self.mtx[0].col

    def data_ptr(self):
        """
        Return the matrix data pointer
        """
        return self.mtx[0].data

    def tolist(self):
        """
        Return the matrix into a list form
        """
        return matrix_to_list(self.mtx)

    def fromlist(self, mtx_):
        """
        Set the matrix from a list
        """
        del_matrix(self.mtx)
        del self.mtx
        self.mtx = new_matrix(mtx_)

    def fromnumpy(self, npmx):
        """
        Set the matrix from a numpy array
        """
        del_matrix(self.mtx)
        del self.mtx
        self.mtx = matrix_from_numpy(npmx)

    def appendrow(self, row_lst_):
        """
        Append a row to the matrix
        """
        row_lst = vector.DVector(row_lst_)
        matrix_append_row(self.mtx, row_lst)
        del row_lst

    def appendcol(self, col_lst_):
        """
        Append a columnt vector to the matrix
        """
        col_lst = vector.DVector(col_lst_)
        matrix_append_col(self.mtx, col_lst)
        del col_lst

    def transpose(self,):
        """
        Transpose the matrix
        """
        mtx_t = init_matrix()
        resize_matrix(mtx_t, self.mtx[0].col, self.mtx[0].row)
        matrix_transpose(self.mtx, mtx_t)
        matrix_copy(mtx_t, self.mtx)
        del_matrix(mtx_t)
        del mtx_t

    def get_evect_evals(self):
        """
        Return the eigenvectors and eigenvalues of the matrix
        """
        mtx_e_vect_ = vector.init_dvector()
        mtx_e_vals_ = init_matrix()
        evect_eval(self.mtx,  mtx_e_vect_, mtx_e_vals_)
        mtx_e_vect = vector.dvector_tolist( mtx_e_vect_)
        mtx_e_vals = matrix_to_list( mtx_e_vals_)
        del_matrix( mtx_e_vals_)
        vector.del_dvector( mtx_e_vect_)
        return  mtx_e_vect,  mtx_e_vals

    def svd(self):
        """
        Return the SVD U,S and VT of the matrix
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


if __name__ in "__main__":
    from random import random
    from random import seed
    seed(123)
    print("Create a matrix of 10x2 of random values")
    print("First prepare a list of list")
    a = [[random() for j in range(2)] for i in range(10)]
    print("Modify the row 4 column 1 with a string'aiuasd'")
    a[4][1] = "aiuasd"

    print("Convert the list of list matrix into matrix type")
    m = Matrix(a)
    m.debug()


    print("Get value at row 1 column 1")
    print(m[1, 1])
    print("Set value at row 1 column 1 with value -2 ")
    m[1, 1] = -2.
    m.debug()
    print("Convert the matrix into a list of list")
    mlst = m.tolist()
    for row in mlst:
        print(row)

    print("Append the row [99, 100]")
    row = [99, 100]
    m.appendrow(row)
    m.debug()

    print("Append the a col [80, 81, 82, 83, 84, 85, 86 ,87, 88, 89, 90]")
    col = [80, 81, 82, 83, 84, 85, 86 ,87, 88, 89, 90]
    m.appendcol(col)
    m.debug()

    print("Test matrix_dvector_dot_product")
    d = vector.DVector([1, 2, 3])
    r = vector.DVector([0 for i in range(11)])
    matrix_dvector_dot_product(m.mtx, d.dvect, r.dvect)
    r.debug()
    del r
    r = vector.DVector([0 for i in range(11)])
    print("Test MT_MatrixDVectorDotProduct")
    mt_matrix_dvector_dot_product(m.mtx, d.dvect, r.dvect)
    r.debug()
    del d
    del r

    print("Test DVectorMatrixDotProduct")
    d = vector.DVector(list(range(1, 12)))
    r = vector.DVector([0 for i in range(3)])
    dvector_matrix_dot_product(m.mtx, d.dvect, r.dvect)
    r.debug()
    del r
    r = vector.DVector([0 for i in range(3)])
    print("Test MT_DVectorMatrixDotProduct")
    mt_dvector_matrix_dot_product(m.mtx, d.dvect, r.dvect)
    r.debug()
    del d
    del r
    del m

    print("Test DVectorTrasposedDVectorDotProduct")
    d1 = vector.DVector(list(range(10)))
    d2 = vector.DVector(list(range(10)))
    m = Matrix([[0 for j in range(10)] for i in range(10)])
    dvector_transposed_dvector_dot_product(d1.dvect, d2.dvect, m.mtx)
    m.debug()
    del d1
    del d2
    del m

    print("Test DVectorTransposedMatrixDivision")
    d = vector.DVector([random() for i in range(10)])
    m = Matrix([[random() for j in range(10)] for i in range(10)])
    r = vector.DVector([0 for i in range(10)])
    dvector_transposed_matrix_division(d.dvect, m.mtx, r.dvect)
    r.debug()
    del d
    del r
    del m


    print("Test MatrixDotProduct")
    m_t = Matrix([[random() for j in range(10)] for i in range(3)])
    m = Matrix([[random() for j in range(3)] for i in range(10)])
    r = Matrix([[0 for j in range(3)] for i in range(3)])
    matrix_dot_product(m_t.mtx, m.mtx, r.mtx)
    r.debug()
    del m_t
    del m
    del r

    print("Test RowColOuterProduct == DVectorTrasposedDVectorDotProduct")
    a = vector.DVector(list(range(10)))
    b = vector.DVector(list(range(10)))
    r = Matrix([[0 for j in range(10)] for i in range(10)])
    row_col_outer_product(a.dvect, b.dvect, r.mtx)
    r.debug()
    del a
    del b
    del r

    print("Test MatrixTranspose")
    m = Matrix([[random() for j in range(3)] for i in range(10)])
    r = Matrix([[0 for j in range(10)] for i in range(3)])
    matrix_transpose(m.mtx, r.mtx)
    m.debug()
    r.debug()

    print("Test MatrixTranspose class inner method")
    m.transpose()
    m.debug()
    m.transpose()

    del m
    del r

    print("Test MatrixInversion")
    m = Matrix([[random() for j in range(4)] for i in range(4)])
    r = Matrix([])
    m.debug()
    matrix_inversion(m.mtx, r.mtx)
    r.debug()
    del r
    del m

    print("Test Eigen vectors/values")
    m = Matrix([[random() for j in range(10)] for i in range(10)])
    e_vect = vector.DVector([])
    e_vals = Matrix([])
    m.debug()
    evect_eval(m.mtx, e_vect.dvect, e_vals.mtx)
    e_vect.debug()
    e_vals.debug()
    del e_vect
    del e_vals
    print("Test Eigen vectors/values inner method")
    e_vect, e_vals = m.get_evect_evals()
    for row in e_vect:
        print(row)
    for row in e_vals:
        print(row)

    u, s, vt = m.svd()

    del m
