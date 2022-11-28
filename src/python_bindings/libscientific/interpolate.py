"""
interpolate - libscientific python binding

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
import libscientific.matrix as mx
import libscientific.vector as vect

lsci = load_libscientific_library()


lsci.cubic_spline_interpolation.argtypes = [ctypes.POINTER(mx.MATRIX),
                                            ctypes.POINTER(mx.MATRIX)]
lsci.cubic_spline_interpolation.restype = None


def cubic_spline_interpolation(x_y, s_interp):
    """
    CubicSplineInterpolation: Given a matrix x_y, interpolate the npoints
                              using a natural cubic spline approach
    """
    x_y_ = None
    if isinstance(x_y, mx.Matrix):
        x_y_ = x_y.mx
    else:
        x_y_ = x_y

    s_interp_ = None
    if isinstance(S, mx.Matrix):
        s_interp_ = s_interp.mx
    else:
        s_interp_ = s_interp

    lsci.cubic_spline_interpolation(x_y_, s_interp_)


lsci.cubic_spline_predict.argtypes = [ctypes.POINTER(vect.DVECTOR),
                                      ctypes.POINTER(mx.MATRIX),
                                      ctypes.POINTER(vect.DVECTOR)]
lsci.cubic_spline_predict.restype = None


def cubic_spline_predict(x_cc, s_interp, y_pred):
    """
    Natural cubic spline prediction
    """
    x_cc_ = None
    if isinstance(x, vect.DVector):
        x_cc_ = x_cc.dvect
    else:
        x_cc_ = x_cc

    s_interp_ = None
    if isinstance(s_interp, mx.Matrix):
        s_interp_ = s_interp.mx
    else:
        s_interp_ = s_interp

    y_pred_ = None
    if isinstance(y_pred, vect.DVector):
        y_pred_ = y_pred.d
    else:
        y_pred_ = y_pred

    lsci.cubic_spline_predict(x_cc_, s_interp_, y_pred_)


if __name__ in "__main__":
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
    print("Origin matrix x->y to interpolate")
    mx.print_matrix(xy)
    S = mx.init_matrix()
    cubic_spline_interpolation(xy, S)
    print("Spline interpolation coefficient matrix")
    mx.print_matrix(S)
    x = vect.DVector(x)
    print("Vector to predict")
    x.debug()
    yp = vect.init_dvector()
    cubic_spline_predict(x, S, yp)

    print("Prediction")
    vect.print_dvector(yp)
    mx.del_matrix(S)
    vect.del_dvector(yp)
    mx.del_matrix(xy)
