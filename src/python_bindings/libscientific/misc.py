# misc libscientific python binding
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

import ctypes
from libscientific.loadlibrary import load_libscientific_library

lsci = load_libscientific_library()


lsci.missing_value.argtypes = None
lsci.missing_value.restype = ctypes.c_double


def missing_value():
    """
    missing_value: return the standard missing value in
    libscientific. T^his variable is useful to
    compute PCA or PLS with missing values.
    """
    return lsci.missing_value()
