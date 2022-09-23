# interpolate - libscientific python binding
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
from libscientific.loadlibrary import LoadLibrary
import libscientific.matrix as mx
import libscientific.vector as vect

lsci = LoadLibrary()


lsci.cubic_spline_interpolation.argtypes = [ctypes.POINTER(mx.matrix),
                                            ctypes.POINTER(mx.matrix)]
lsci.cubic_spline_interpolation.restype = None


def cubic_spline_interpolation(xy, S):
    """
    CubicSplineInterpolation: Given a matrix xy, interpolate the npoints
                              using a natural cubic spline approach
    """
    xy_ = None
    if "LP_matrix" in str(type(xy)):
        xy_ = xy
    elif "Matrix" in str(type(xy)):
        xy_ = xy.mx

    S_ = None
    if "LP_matrix" in str(type(S)):
        S_ = S
    elif "Matrix" in str(type(S)):
        S_ = S.mx

    lsci.cubic_spline_interpolation(xy_, S_)


lsci.cubic_spline_predict.argtypes = [ctypes.POINTER(vect.dvector),
                                      ctypes.POINTER(mx.matrix),
                                      ctypes.POINTER(vect.dvector)]
lsci.cubic_spline_predict.restype = None


def cubic_spline_predict(x, S, yp):
    """
    Natural cubic spline prediction
    """
    x_ = None
    if "DVector" in str(type(x)):
        x_ = x.d
    else:
        x_ = x

    S_ = None
    if "LP_matrix" in str(type(S)):
        S_ = S
    elif "Matrix" in str(type(S)):
        S_ = S.mx

    yp_ = None
    if "DVector" in str(type(yp)):
        yp_ = yp.d
    else:
        yp_ = yp

    lsci.cubic_spline_predict(x_, S_, yp_)


if __name__ in "__main__":
    xy_ = [[-1.82,0.63], [-0.73,0.19], [-0.17,0.01], [-0.09,0.00], [0.15,0.01], [0.39,0.06], [0.86,0.24], [1.44,0.49]]
    x_ = [-1.82, -0.73, -0.17, -0.09, 0.15, 0.39, 0.86, 1.44]
    xy = mx.NewMatrix(xy_)
    print("Origin matrix x->y to interpolate")
    mx.PrintMatrix(xy)
    S = mx.initMatrix()
    cubic_spline_interpolation(xy, S)
    print("Spline interpolation coefficient matrix")
    mx.PrintMatrix(S)
    x = vect.DVector(x_)
    print("Vector to predict")
    x.debug()
    yp = vect.initDVector()
    cubic_spline_predict(x, S, yp)

    print("Prediction")
    vect.PrintDVector(yp)
    mx.DelMatrix(S)
    vect.DelDVector(yp)
