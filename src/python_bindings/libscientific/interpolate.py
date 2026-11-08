# interpolate.py libscientific python binding
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
"""Natural cubic spline interpolation.

>>> from libscientific.interpolate import CubicSplineInterpolation
>>> spline = CubicSplineInterpolation()
>>> spline.fit([[0.0, 0.0], [1.0, 1.0], [2.0, 4.0], [3.0, 9.0]])
>>> round(spline.predict([2.0])[0], 6)
4.0
"""

import ctypes

import libscientific.matrix as mx
import libscientific.vector as vect
from libscientific.loadlibrary import load_libscientific_library

lsci = load_libscientific_library()


lsci.cubic_spline_interpolation.argtypes = [ctypes.POINTER(mx.MATRIX), ctypes.POINTER(mx.MATRIX)]
lsci.cubic_spline_interpolation.restype = None


def cubic_spline_interpolation(x_y, s_interp):
    """Compute the coefficients of a natural cubic spline (low-level API).

    Parameters
    ----------
    x_y : Matrix or ctypes.POINTER(MATRIX)
        The points, shape ``(n_points, 2)``: column 0 is x (**sorted in
        increasing order**, no duplicates), column 1 is y.
    s_interp : Matrix or ctypes.POINTER(MATRIX)
        Output matrix; receives one row per interval, ``(n_points - 1, 5)``:
        ``[x_i, a, b, c, d]`` with
        ``y = a + b (x - x_i) + c (x - x_i)^2 + d (x - x_i)^3`` on
        ``[x_i, x_{i+1}]``.
    """
    x_y_ = None
    if isinstance(x_y, mx.Matrix):
        x_y_ = x_y.mtx
    else:
        x_y_ = x_y

    s_interp_ = None
    if isinstance(s_interp, mx.Matrix):
        s_interp_ = s_interp.mtx
    else:
        s_interp_ = s_interp

    lsci.cubic_spline_interpolation(x_y_, s_interp_)


lsci.cubic_spline_predict.argtypes = [
    ctypes.POINTER(vect.DVECTOR),
    ctypes.POINTER(mx.MATRIX),
    ctypes.POINTER(vect.DVECTOR),
]
lsci.cubic_spline_predict.restype = None


def cubic_spline_predict(x_cc, s_interp, y_pred):
    """Evaluate a natural cubic spline (low-level API).

    Parameters
    ----------
    x_cc : DVector or ctypes.POINTER(DVECTOR)
        The x values to evaluate.
    s_interp : Matrix or ctypes.POINTER(MATRIX)
        Spline coefficients from :func:`cubic_spline_interpolation`.
    y_pred : DVector or ctypes.POINTER(DVECTOR)
        Output vector; resized and filled with one y value per x.
    """
    x_cc_ = None
    if isinstance(x_cc, vect.DVector):
        x_cc_ = x_cc.dvect
    else:
        x_cc_ = x_cc

    s_interp_ = None
    if isinstance(s_interp, mx.Matrix):
        s_interp_ = s_interp.mtx
    else:
        s_interp_ = s_interp

    y_pred_ = None
    if isinstance(y_pred, vect.DVector):
        y_pred_ = y_pred.dvect
    else:
        y_pred_ = y_pred

    lsci.cubic_spline_predict(x_cc_, s_interp_, y_pred_)


class CubicSplineInterpolation:
    """Natural cubic spline interpolation of 2D points.

    The spline goes exactly through every given point, with continuous first
    and second derivatives, and zero second derivative at the two ends.

    Attributes
    ----------
    s_interp : ctypes.POINTER(MATRIX) or None
        The spline coefficients (see :func:`cubic_spline_interpolation`),
        ``None`` before :meth:`fit`.

    Examples
    --------
    >>> import math
    >>> from libscientific.interpolate import CubicSplineInterpolation
    >>> points = [[x / 2, math.sin(x / 2)] for x in range(13)]   # x from 0 to 6
    >>> spline = CubicSplineInterpolation()
    >>> spline.fit(points)
    >>> [round(y, 3) for y in spline.predict([0.25, 1.1, 3.3])]
    [0.247, 0.891, -0.158]
    """

    def __init__(self) -> None:
        """Create an empty (not fitted) spline."""
        self.s_interp = None

    def __repr__(self):
        return f"CubicSplineInterpolation(fitted={self.s_interp is not None})"

    def __str__(self):
        return f"CubicSplineInterpolation(fitted={self.s_interp is not None})"

    def __del__(self):
        """Release the C spline coefficients."""
        if getattr(self, "s_interp", None) is not None:
            mx.del_matrix(self.s_interp)
        self.s_interp = None

    def fit(self, x_y_) -> None:
        """Compute the spline through a set of points.

        Calling ``fit`` again replaces the previous spline.

        Parameters
        ----------
        x_y_ : MatrixLike
            The points, shape ``(n_points, 2)``: ``[[x0, y0], [x1, y1], ...]``.
            Lists and numpy arrays are sorted by x automatically; a
            :class:`~libscientific.matrix.Matrix` must already be sorted. The x values
            must be distinct.

        Returns
        -------
        None
        """
        if self.s_interp is not None:
            mx.del_matrix(self.s_interp)
        self.s_interp = mx.init_matrix()
        if isinstance(x_y_, (mx.Matrix, ctypes.POINTER(mx.MATRIX))):
            x_y, xalloc = mx.as_c_matrix(x_y_)
        else:
            # the spline needs the points sorted by x
            points = x_y_.tolist() if "numpy" in str(type(x_y_)) else list(x_y_)
            x_y, xalloc = mx.as_c_matrix(sorted(points, key=lambda p: p[0]))
        try:
            cubic_spline_interpolation(x_y, self.s_interp)
        finally:
            if xalloc:
                mx.del_matrix(x_y)

    def predict(self, x_pred_):
        """Evaluate the spline at new x values.

        Parameters
        ----------
        x_pred_ : list[float], DVector or ctypes.POINTER(DVECTOR)
            The x values. Values outside the range of the fitted points are
            extrapolated with the first/last polynomial and are not reliable.

        Returns
        -------
        list[float]
            One y value per x.

        Raises
        ------
        RuntimeError
            If the spline has not been fitted.
        """
        if self.s_interp is None:
            raise RuntimeError("the model is not fitted: call fit() first")
        xalloc = False
        if isinstance(x_pred_, vect.DVector):
            x_pred = x_pred_.dvect
        elif isinstance(x_pred_, ctypes.POINTER(vect.DVECTOR)):
            x_pred = x_pred_
        else:
            x_pred = vect.new_dvector(list(x_pred_))
            xalloc = True

        ypred_ = vect.init_dvector()
        cubic_spline_predict(x_pred, self.s_interp, ypred_)
        ypred = vect.dvector_tolist(ypred_)
        vect.del_dvector(ypred_)
        if xalloc:
            vect.del_dvector(x_pred)
        return ypred
