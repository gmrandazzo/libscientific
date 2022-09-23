# tensor libscientific python binding
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
from libscientific import matrix as mx
from libscientific import misc

lsci = LoadLibrary()


class tensor(ctypes.Structure):
    _fields_ = [
        ("m",    ctypes.POINTER(ctypes.POINTER(mx.matrix))),
        ("order",     ctypes.c_size_t)]


lsci.initTensor.argtypes = [ctypes.POINTER(ctypes.POINTER(tensor))]
lsci.initTensor.restype = None


def initTensor():
    """
    initTensor: Allocate in memory an empty libscientific tensor
    """
    t = ctypes.POINTER(tensor)()
    # m = matrix()
    lsci.initTensor(ctypes.pointer(t))
    return t

# Craking the code interview


lsci.NewTensorMatrix.argtypes = [ctypes.POINTER(ctypes.POINTER(tensor)),
                                 ctypes.c_size_t,
                                 ctypes.c_size_t,
                                 ctypes.c_size_t]
lsci.NewTensorMatrix.restype = None


def NewTensorMatrix(t, k, row, col):
    lsci.NewTensorMatrix(ctype.pointer(t), k, row, col);

lsci.AddTensorMatrix.argtypes = [ctypes.POINTER(ctypes.POINTER(tensor)),
                                 ctypes.c_size_t,
                                 ctypes.c_size_t]
lsci.AddTensorMatrix.restype = None

lsci.NewTensor.argtypes = [ctypes.POINTER(ctypes.POINTER(tensor)),
                           ctypes.c_size_t]
lsci.NewTensor.restype = None


def NewTensor(a_):
    """
    NewTensorMatrix: Allocate in memory a libscientific tensor
    from a list of lists of lists
    """
    a = None
    if "numpy" in str(type(a_)):
        a = a_.tolist()
    else:
        a = a_
    order = None
    nrows = None
    ncols = None
    try:
        order = len(a)
        try:
            nrows = len(a[0])
            ncols = len(a[0][0])
        except IndexError:
            nrows = 0
            ncols = 0
    except IndexError:
        order = 0
    t = ctypes.POINTER(tensor)()
    lsci.NewTensor(ctypes.pointer(t),
                   order)
    for k in range(order):
        lsci.NewTensorMatrix()
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
