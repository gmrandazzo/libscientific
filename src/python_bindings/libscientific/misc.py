# misc.py libscientific python binding
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
"""Miscellaneous utilities: the missing value marker."""

import ctypes

from libscientific.loadlibrary import load_libscientific_library

lsci = load_libscientific_library()


lsci.missing_value.argtypes = []
lsci.missing_value.restype = ctypes.c_double


def missing_value():
    """Return the value that marks a missing data point in libscientific.

    Data points equal to this value (within a tolerance) are skipped by the
    algorithms (PCA, PLS, statistics, ...). In Python you can also pass
    ``float("nan")``: NaN and infinite values are converted to this value by the
    C library. Non-numeric entries of lists (e.g. ``"NA"``) are converted to it
    when a matrix or vector is created.

    Returns
    -------
    float
        The missing value marker (``99999999.0``).

    Examples
    --------
    >>> missing_value()
    99999999.0
    """
    return lsci.missing_value()
