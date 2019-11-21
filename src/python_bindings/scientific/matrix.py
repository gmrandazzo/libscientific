# matrix libscientific python binding
#
# Copyright (C) <2019>  Giuseppe Marco Randazzo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
import ctypes
from scientific.loadlibrary import LoadLibrary
from scientific import misc

lsci = LoadLibrary()


class matrix(ctypes.Structure):
    _fields_ = [
        ("data",    ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
        ("row",     ctypes.c_size_t),
        ("col",     ctypes.c_size_t)]


lsci.initMatrix.argtypes = [ctypes.POINTER(ctypes.POINTER(matrix))]
lsci.initMatrix.restype = None


def initMatrix():
    """
    initMatrix: Allocate in memory an empty libscientific matrix
    """
    m = ctypes.POINTER(matrix)()
    # m = matrix()
    lsci.initMatrix(ctypes.pointer(m))
    return m


lsci.NewMatrix.argtypes = [ctypes.POINTER(ctypes.POINTER(matrix)),
                           ctypes.c_size_t,
                           ctypes.c_size_t]
lsci.NewMatrix.restype = None


def NewMatrix(a_):
    """
    NewMatrix: Allocate in memory a libscientific matrix from a list of lists
    """
    a = None
    if "numpy" in str(type(a_)):
        a = a_.tolist()
    else:
        a = a_
    nrows = None
    ncols = None
    try:
        nrows = len(a)
        try:
            ncols = len(a[0])
        except IndexError:
            ncols = 0
    except IndexError:
        nrows = 0
    m = ctypes.POINTER(matrix)()
    lsci.NewMatrix(ctypes.pointer(m),
                   nrows,
                   ncols)

    for i in range(nrows):
        for j in range(ncols):
            val = None
            try:
                val = float(a[i][j])
            except ValueError:
                val = None

            if val is None:
                setMissingMatrixValue(m, i, j)
            else:
                lsci.setMatrixValue(m, i, j, val)
    return m


lsci.ResizeMatrix.argtypes = [ctypes.POINTER(ctypes.POINTER(matrix)),
                              ctypes.c_size_t,
                              ctypes.c_size_t]
lsci.ResizeMatrix.restype = None


def ResizeMatrix(m, nrows, ncols):
    """
    ResizeMatrix: Resize an already allocated or reallocate a libscientific
                  matrix
    """
    lsci.ResizeMatrix(ctypes.pointer(m), nrows, ncols)


lsci.DelMatrix.argtypes = [ctypes.POINTER(ctypes.POINTER(matrix))]
lsci.DelMatrix.restype = None


def DelMatrix(m):
    """
    DelMatrix: Delete an allocated libscientific matrix
    """
    lsci.DelMatrix(ctypes.pointer(m))


lsci.MatrixCheck.argtypes = [ctypes.POINTER(matrix)]
lsci.MatrixCheck.restype = None


def MatrixCheck(m):
    """
    MatrixCheck: Find infinite and nan numbers and sobstitute
                 with MISSING value.
    """
    lsci.DelMatrix(ctypes.pointer(m))


lsci.PrintMatrix.argtypes = [ctypes.POINTER(matrix)]
lsci.PrintMatrix.restype = None


def PrintMatrix(m):
    """
    PrintMatrix: Print to video a libscientific matrix
    """
    lsci.PrintMatrix(m)


lsci.ValInMatrix.argtypes = [ctypes.POINTER(matrix), ctypes.c_double]
lsci.ValInMatrix.restype = ctypes.c_int


def ValInMatrix(m, val):
    """
    ValInMatrix: Check if a libscientific matrix contains an exact value "val"
                 and return 1 or 0 respectivelly for yes or no.
    """
    return lsci.ValInMatrix(m, val)


lsci.MatrixSet.argtypes = [ctypes.POINTER(matrix), ctypes.c_double]
lsci.MatrixSet.restype = None


def MatrixSet(m, val):
    """
    MatrixSet: Set all values of a libscientific matrix to "val"
    """
    lsci.MatrixSet(m, val)


lsci.MatrixCopy.argtypes = [ctypes.POINTER(matrix),
                            ctypes.POINTER(ctypes.POINTER(matrix))]
lsci.MatrixCopy.restype = None


def MatrixCopy(msrc, mdst):
    """
    MatrixCopy: Copy a libscientifi matrix to another allocated one
    """
    lsci.MatrixCopy(msrc, ctypes.pointer(mdst))


lsci.setMatrixValue.argtypes = [ctypes.POINTER(matrix),
                                ctypes.c_size_t,
                                ctypes.c_size_t,
                                ctypes.c_double]
lsci.setMatrixValue.restype = None


def setMatrixValue(m, irow, jcol, value):
    """
    setMatrixValue: Set/modify a value in the irow and jcol of a libscientific
                    matrix
    """
    lsci.setMatrixValue(m, irow, jcol, value)


lsci.getMatrixValue.argtypes = [ctypes.POINTER(matrix),
                                ctypes.c_size_t,
                                ctypes.c_size_t]
lsci.getMatrixValue.restype = ctypes.c_double


def getMatrixValue(m, irow, jcol):
    """
    getMatrixValue: Get a value in the irow and jcol of a libscientific
                    matrix
    """
    return lsci.getMatrixValue(m, irow, jcol)


def getMatrixRow(m, irow):
    """
    getMatrixRow: Python version of getting a row from a libscientific matrix
    """
    row = []
    for j in range(m[0].col):
        row.append(getMatrixValue(m, irow, j))
    return row


def getMatrixColumn(m, jcol):
    """
    getMatrixColumn: Python version of getting a column from a libscientific
                     matrix
    """
    col = []
    for i in range(m[0].row):
        col.append(getMatrixValue(m, i, jcol))
    return col


def MatrixToList(m):
    """
    MatrixToList: Convert a libscientific matrix to list of list
    """
    mlst = []
    try:
        for i in range(m[0].row):
            row = []
            for j in range(m[0].col):
                row.append(m[0].data[i][j])
            mlst.append(row)
    except TypeError:
        for i in range(m.row):
            row = []
            for j in range(m.col):
                row.append(m.data[i][j])
            mlst.append(row)
    return mlst


def MatrixFromNumpy(npm):
    return NewMatrix(npm.tolist())


def setMissingMatrixValue(m, row, col):
    m[0].data[row][col] = misc.missing_value()


"""
/* Description:
 * Append double vector uivector as row.
 */
void MatrixAppendRow(matrix **mx, dvector *row);

/* Description:
 * Append double vector uivector as column.
 */
void MatrixAppendCol(matrix **mx, dvector *col);
/* Description:
 * Append unsigned int vector uivector as row
 */
void MatrixAppendUIRow(matrix **mx, uivector *row);

/* Description:
 * Append unsigned int vector uivector as column
 */
void MatrixAppendUIRow(matrix **mx, uivector *row);

/*  Matrix Operations */

/* Description:
 * matrix - row double vector product: the result is a row double vector
 * i.e.: X(10x5) * d(5x1) = r(10x1)
 */
void MatrixDVectorDotProduct(matrix *m, dvector *v, dvector *r);

/* Multithread version of MatrixDVectorDotProduct */
void MT_MatrixDVectorDotProduct(matrix *mx, dvector *v, dvector *p);

/* Description:
 * column double vector - matrix product: the result is a column double vector
 * i.e.: d(1x5) * X(5x10) = r(1x10)
 */
void DVectorMatrixDotProduct(matrix *mx, dvector *v, dvector *p);

/* Multithread version of DVectorMatrixDotProduct */
void MT_DVectorMatrixDotProduct(matrix *mx, dvector *v, dvector *p);

/* Description:
 * transposed double vector - double vecrtor product: the result is a matrix
 * i.e.: X  = d'd
 * i.e.: d'(5x1) * d(1x5) = X(1x10)
 */
void DVectorTrasposedDVectorDotProduct(dvector *v1, dvector *v2, matrix *m);

/* r = v/mx  = (inv(mx^T)*v^T)^T*/
void DVectorTransposedMatrixDivision(dvector *v, matrix *mx, dvector *r);

/*
 * Description:
 * Calculate the matrix matrix product
 */
void MatrixDotProduct(matrix *m_t, matrix *m, matrix *r);
void RowColOuterProduct(dvector *a, dvector *b, matrix *m);

/*
 * Description:
 * Generate a transpose matrix of m
 */
void MatrixTranspose(matrix *m, matrix *r);

/*
 * Description:
 * Matrix inversion using the Gauss-Jordan algorithm
 */
void MatrixInversion(matrix *m, matrix **m_inv);

/*
 * Description:
 * Matrix inversion using the LU decompositio
 */
void MatrixLUInversion(matrix *m, matrix **m_inv);

/*
 * Description:
 * Matrix pseudo inversion using the SVD algorithm
 */
void MatrixPseudoinversion(matrix *m, matrix **m_inv);

/*
 * Description:
 * Matrix pseudo inversion using the Moore-Penrose pseudoinverse
 */
void MatrixMoorePenrosePseudoinverse(matrix *m, matrix **inv);

/*
 * Description:
 * Generate the identity matrix
 */
void GenIdentityMatrix(matrix **m);

/*
 * Description:
 * Calculate the mean centered matrix
 */
void MeanCenteredMatrix(matrix *mx, matrix *mxc);

/*
 * Description:
 * Calculate the pearson correlation matrix
 */
void PearsonCorrelMatrix(matrix *mxsrc, matrix *mxdst);

/*
 * Description:
 * Calculate the spearmann correlation matrix
 */
void SpearmanCorrelMatrix(matrix *mxsrc, matrix *mxdst);

/*Calculate the column Average for the matrix mx*/
void MatrixColAverage(matrix *mx, dvector **colaverage);

/*Calculate the row Average for the matrix mx*/
void MatrixRowAverage(matrix *mx, dvector **rowaverage);

/*Calculate the column Standard Deviation for the matrix mx*/
void MatrixColSDEV(matrix *mx, dvector **colsdev);

/*Calculate the column Root Mean Square for the matrix mx*/
void MatrixColRMS(matrix* mx, dvector** colrms);

/*Calculate the column variance for the matrix mx*/
void MatrixColVar(matrix *mx, dvector **colvar);

/* Calculate the matrix descriptive statistics:
 *  - Column Average
 *  - Column Median
 *  - Column Armonic Average
 *  - Column Variance Population
 *  - Column Variance Sample (Correcter Variance)
 *  - Column Standard Deviation
 *  - Column Standard Deviation Sample (Corrected Standard Deviation)
 *  - Column Max
 *  - Column Min
 */
void MatrixColDescStat(matrix *mx, matrix **ds);

/*Calculate the covariance matrix*/
void MatrixCovariance(matrix *mx, matrix **cm);

/* Transform a matrix into a logaritmic matrix */
void Matrix2LogMatrix(matrix *mx_in, matrix **mx_out);

/* Transform a matrix into a SQUARE matrix */
void Matrix2SquareMatrix(matrix *mx_in, matrix **mx_out);

/* Transform a matrix into a SQRT matrix */
void Matrix2SQRTMatrix(matrix *mx_in, matrix **mx_out);

/* Transform a matrix into ABS matrix */
void Matrix2ABSMatrix(matrix *mx_in, matrix **mx_out);

/*
 * Description:
 * Develop an interaction factors matrix
 * Es. Use in DOE
 */
void Matrix2IntFactorsMatrix(matrix *mx_in, size_t factors, matrix **mx_out);

/*
 * Description:
 * Transform a matrix into a row centered scaled matrix
 * Es. Use in Spectroscopy
 */
void MatrixRowCenterScaling(matrix *mx_in, matrix **mx_out);

/*
 * Description:
 * Transform a matrix into a SVN row scaled matrix
 * Es. Use in Spectroscopy
 */
void MatrixSVNScaling(matrix *mx_in, matrix **mx_out);

/*
 * Description:
 * calculate the square root of the sum of the squares
 * of all the elements in the matrix
 * ||X|| = double
 */
double Matrixnorm(matrix *mx);
double Matrix1norm(matrix *mx);

/*
 * Description:
 * calculate the determinant of a matrix
 */
double MatrixDeterminant(matrix *mx);

/*
 * Description:
 * Normalize the matrix for the Matrixnorm value.
 * Each value of mx is divided by double Matrixnorm(matrix *mx);
 */
void MatrixNorm(matrix *mx, matrix *nmx);

/*
 * Description:
 * Find the minimum and maximum of a column in matrix
 */
void MatrixColumnMinMax(matrix* mx, size_t col, double* min, double* max);

/*
 * Description:
 * Sort of a matrix ba a column number col_n
 */
void MatrixSort(matrix *mx, size_t col_n);

/*
 * Description:
 * Reverse sort of a matrix by a column number col_n
 */
void MatrixReverseSort(matrix* mx, size_t col_n);

/*
 * Description:
 * find the maximum value in matrix and return the row and col indexes
 */
void MatrixGetMaxValueIndex(matrix *mx, size_t *row, size_t *col);

/*
 * Description:
 * find the minimum value in matrix and return the row and col indexes
 */
void MatrixGetMinValueIndex(matrix *mx, size_t *row, size_t *col);


/*
 * Description:
 * Singular Value Decomposition local implementation
 */
void SVD(matrix* mx, matrix **U, matrix **S, matrix **VT);

/*
 * Description:
 * Singular Value Decomposition lapack implementation
 */
void SVDlapack(matrix *mx, matrix **u, matrix **s, matrix **vt);

/*
 * Description:
 * Eigenvectors and Eigenvalues with the QR Method
 */
void EVectEval(matrix *mx, dvector **eval, matrix **evect);

/*
 * Description:
 * QR Decomposition with the householder method
 */
void QRDecomposition(matrix *mx, matrix **Q, matrix **R);

/*
 * Description:
 * LU Decomposition with the householder method
 */
void LUDecomposition(matrix *mx, matrix **L, matrix **U);

/*
 * Description:
 * Householder Reflector vector
 */
void HouseReflectorVect(dvector *x, dvector *u);

/*
 * Description:
 * Householder Vector Matrix Multiplication
 */
void HouseholderVectMatrixProduct(dvector *h, matrix *A, matrix *P);

/*
 * Description:
 * Householder matrix
 */
void HouseholderMatrix(dvector *v, matrix *h);

/*
 * Description:
 *  Reduce a matrix to an Hessenberg form
 * - Householder implementation
 * - Cholesky implementation
 */
void HouseholderReduction(matrix *mx);
void CholeskyReduction(matrix *m);
"""


class Matrix(object):
    """
    Translate a list of list into a libscientific matrix
    """
    def __init__(self, mx_):
        self.mx = NewMatrix(mx_)

    def __del__(self):
        DelMatrix(self.mx)
        del self.mx
        self.mx = None

    def __getitem__(self, keys):
        i, j = keys
        return self.data_ptr()[i][j]

    def __setitem__(self, keys, value):
        i, j = keys
        setMatrixValue(self.mx, i, j, value)

    def nrow(self):
        return self.mx[0].row

    def ncol(self):
        return self.mx[0].col

    def data_ptr(self):
        return self.mx[0].data

    def tolist(self):
        return MatrixToList(self.mx)

    def fromlist(self, mx_):
        DelMatrix(self.mx)
        del self.mx
        self.mx = NewMatrix(mx_)

    def fromnumpy(self, npmx):
        DelMatrix(self.mx)
        del self.mx
        self.mx = MatrixFromNumpy(npmx)

    def debug(self):
        PrintMatrix(self.mx)


if __name__ in "__main__":
    from random import random
    a = [[random() for j in range(2)] for i in range(10)]
    a[4][1] = "aiuasd"
    m = Matrix(a)
    print("Get value example")
    print(m[1, 1])
    print("Set value example")
    m[1, 1] = -2.
    m.debug()
    mlst = m.tolist()
    for row in mlst:
        print(row)
