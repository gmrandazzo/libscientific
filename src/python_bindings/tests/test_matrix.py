from libscientific.matrix import *
from libscientific.vector import *
from libscientific.misc import *
from utils import(matrix_sum,
                  raw_vector_sum)
import random

def test_matrix():
    random.seed(123)
    #print("Create a matrix of 10x2 of random values")
    #print("First prepare a list of list")
    a = [[random.random() for j in range(2)] for i in range(10)]
    #print("Modify the row 4 column 1 with a string'aiuasd'")
    a[4][1] = "aiuasd"

    #print("Convert the list of list matrix into matrix type")
    m = Matrix(a)
    #m.debug()
    assert abs(m[4,1]-missing_value()) <= 1e-14

    #print("Get value at row 1 column 1")
    print(m[1, 1])
    #print("Set value at row 1 column 1 with value -2 ")
    m[1, 1] = -2.
    #m.debug()
    #print("Convert the matrix into a list of list")
    mlst = m.tolist()
    for row in mlst:
        print(row)

    #print("Append the row [99, 100]")
    row = [99, 100]
    m.appendrow(row)
    #m.debug()

    #print("Append the a col [80, 81, 82, 83, 84, 85, 86 ,87, 88, 89, 90]")
    col = [80, 81, 82, 83, 84, 85, 86 ,87, 88, 89, 90]
    m.appendcol(col)
    #m.debug()

    #print("Test matrix_dvector_dot_product")
    d = vector.DVector([1, 2, 3])
    r = vector.DVector([0 for i in range(11)])
    matrix_dvector_dot_product(m.mtx, d.dvect, r.dvect)
    #r.debug()
    del r
    r = vector.DVector([0 for i in range(11)])
    #print("Test MT_MatrixDVectorDotProduct")
    mt_matrix_dvector_dot_product(m.mtx, d.dvect, r.dvect)
    #r.debug()
    del d
    del r

    #print("Test DVectorMatrixDotProduct")
    d = vector.DVector(list(range(1, 12)))
    r = vector.DVector([0 for i in range(3)])
    dvector_matrix_dot_product(m.mtx, d.dvect, r.dvect)
    #r.debug()
    del r
    r = vector.DVector([0 for i in range(3)])
    #print("Test MT_DVectorMatrixDotProduct")
    mt_dvector_matrix_dot_product(m.mtx, d.dvect, r.dvect)
    #r.debug()
    del d
    del r
    del m

    #print("Test DVectorTrasposedDVectorDotProduct")
    d1 = vector.DVector(list(range(10)))
    d2 = vector.DVector(list(range(10)))
    m = Matrix([[0 for j in range(10)] for i in range(10)])
    dvector_transposed_dvector_dot_product(d1.dvect, d2.dvect, m.mtx)
    #m.debug()
    del d1
    del d2
    del m

    #print("Test DVectorTransposedMatrixDivision")
    d = vector.DVector([random.random() for i in range(10)])
    m = Matrix([[random.random() for j in range(10)] for i in range(10)])
    r = vector.DVector([0 for i in range(10)])
    dvector_transposed_matrix_division(d.dvect, m.mtx, r.dvect)
    #r.debug()
    del d
    del r
    del m


    #print("Test MatrixDotProduct")
    m_t = Matrix([[random.random() for j in range(10)] for i in range(3)])
    m = Matrix([[random.random() for j in range(3)] for i in range(10)])
    r = Matrix([[0 for j in range(3)] for i in range(3)])
    matrix_dot_product(m_t.mtx, m.mtx, r.mtx)
    #r.debug()
    del m_t
    del m
    del r

    #print("Test RowColOuterProduct == DVectorTrasposedDVectorDotProduct")
    a = vector.DVector(list(range(10)))
    b = vector.DVector(list(range(10)))
    r = Matrix([[0 for j in range(10)] for i in range(10)])
    row_col_outer_product(a.dvect, b.dvect, r.mtx)
    #r.debug()
    del a
    del b
    del r

    #print("Test MatrixTranspose")
    m = Matrix([[random.random() for j in range(3)] for i in range(10)])
    r = Matrix([[0 for j in range(10)] for i in range(3)])
    matrix_transpose(m.mtx, r.mtx)
    #m.debug()
    #r.debug()

    #print("Test MatrixTranspose class inner method")
    m.transpose()
    #m.debug()
    m.transpose()

    del m
    del r

    #print("Test MatrixInversion")
    m = Matrix([[random.random() for j in range(4)] for i in range(4)])
    r = Matrix([])
    #m.debug()
    matrix_inversion(m.mtx, r.mtx)
    #r.debug()
    del r
    del m

    #print("Test Eigen vectors/values")
    m = Matrix([[random.random() for j in range(10)] for i in range(10)])
    e_vect = vector.DVector([])
    e_vals = Matrix([])
    #m.debug()
    evect_eval(m.mtx, e_vect.dvect, e_vals.mtx)
    #e_vect.debug()
    #e_vals.debug()
    del e_vect
    del e_vals
    #print("Test Eigen vectors/values inner method")
    e_vect, e_vals = m.get_evect_evals()
    #for row in e_vect:
    #    print(row)
    #for row in e_vals:
    #    print(row)
    u, s, vt = m.svd()
    del m
    m = Matrix()
    m.appendrow([1,2,3,4])
    m.appendrow([5,6,7,8])
    #m.debug()
