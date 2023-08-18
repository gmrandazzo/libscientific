"""interpolate.py libscientific python binding

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
    if isinstance(s_interp, mx.Matrix):
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
    if isinstance(x_cc, vect.DVector):
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

class CubicSplineInterpolation:
    """
    Cubic Spline Interpolation

    This class provides a cubic spline interpolation model for fitting and predicting values.

    Examples
    --------
    >>> from libscientific import interpolate
    >>> xy = [[-1.82, 0.63],
    ...       [-0.73, 0.19],
    ...       [-0.17, 0.01],
    ...       [-0.09, 0.00],
    ...       [0.15, 0.01],
    ...       [0.39, 0.06],
    ...       [0.86, 0.24],
    ...       [1.44, 0.49]]
    >>> x = [-1.82, -0.73, -0.17, -0.09, 0.15, 0.39, 0.86, 1.44]
    >>> csi = interpolate.CubicSplineInterpolation()
    >>> csi.fit(xy)
    >>> yp = csi.predict(x)
    >>> yp
    [0.63, 0.19, 0.010000000000000002, 1.463672932855431e-18, 0.01, 0.06, 0.23999999999999996, 0.49]

    Methods
    -------
    fit(XY):
        Fit the cubic spline interpolation model with a 2D matrix x and y.
    predict(X):
        Predict targets of given samples.

    Attributes
    ----------
    s_interp : CDataType
        The interpolation matrix.

    """

    def __init__(self) -> None:
        """
        Initialize a CubicSplineInterpolation instance.
        """
        self.s_interp = None

    def __del__(self):
        """
        Clean up resources associated with the CubicSplineInterpolation instance.
        """
        if self.s_interp is not None:
            mx.del_matrix(self.s_interp)
            del self.s_interp
        self.s_interp = None

    def fit(self, x_y_) -> None:
        """
        Fit the cubic spline interpolation model.

        Parameters
        ----------
        x_y_ : Matrix or List[List]
            2D matrix with x and y values for interpolation.

        Returns
        -------
        None
        """
        self.s_interp = mx.init_matrix()
        if isinstance(x_y_, mx.Matrix):
            x_y = x_y_.mx
        else:
            x_y = mx.new_matrix(x_y_)
        cubic_spline_interpolation(x_y, self.s_interp)
        mx.del_matrix(x_y)

    def predict(self, x_pred_):
        """
        Predict values using the cubic spline interpolation model.

        Parameters
        ----------
        x_pred_ : Matrix or List[float]
            Matrix or list of x values for prediction.

        Returns
        -------
        List[float]
            Predicted y values for the given x values.
        """
        if isinstance(x_pred_, mx.Matrix):
            x_pred = x_pred_.d
        else:
            x_pred = vect.new_dvector(x_pred_)

        ypred_ = vect.init_dvector()
        cubic_spline_predict(x_pred, self.s_interp, ypred_)
        ypred = vect.dvector_tolist(ypred_)
        vect.del_dvector(x_pred)
        return ypred
