from libscientific.interpolate import * 
from libscientific.matrix import *
from libscientific.vector import *
from utils import(matrix_sum,
                  raw_vector_sum)


def test_interpolate():
    xy = [[-1.82,0.63],
           [-0.73,0.19],
           [-0.17,0.01],
           [-0.09,0.00],
           [0.15,0.01],
           [0.39,0.06],
           [0.86,0.24],
           [1.44,0.49]]
    x = [-1.82, -0.73, -0.17, -0.09, 0.15, 0.39, 0.86, 1.44]
    xy = mx.new_matrix(xy)
    S = mx.init_matrix()
    cubic_spline_interpolation(xy, S)
    check_value_1 = 1.6251152850133066
    assert abs(matrix_sum(matrix_to_list(S))-check_value_1) <= 1e-14
    x = vect.DVector(x)
    yp = vect.init_dvector()
    cubic_spline_predict(x, S, yp)
    check_value_2 = 1.6300000000000001
    assert abs(raw_vector_sum(yp)-check_value_2) <= 1e-14
    mx.del_matrix(S)
    vect.del_dvector(yp)
    mx.del_matrix(xy)
