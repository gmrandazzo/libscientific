"""
tensor libscientific python binding

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
from libscientific.loadlibrary import LoadLibrary
from libscientific import misc
from libscientific import matrix as mx
from libscientific import vector as vect
from libscientific import misc

lsci = LoadLibrary()

class tensor(ctypes.Structure):
    """
    tensor data structure
    """
    _fields_ = [
        ("m", ctypes.POINTER(ctypes.POINTER(mx.matrix))),
        ("order", ctypes.c_size_t)]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__

lsci.initTensor.argtypes = [ctypes.POINTER(ctypes.POINTER(tensor))]
lsci.initTensor.restype = None

def initTensor():
    """
    initTensor: Allocate in memory an empty libscientific tensor
    """
    t = ctypes.POINTER(tensor)()
    lsci.initTensor(ctypes.pointer(t))
    return t


lsci.NewTensorMatrix.argtypes = [ctypes.POINTER(tensor),
                                 ctypes.c_size_t,
                                 ctypes.c_size_t,
                                 ctypes.c_size_t]
lsci.NewTensorMatrix.restype = None

def NewTensorMatrix(t, k, row, col):
    lsci.NewTensorMatrix(t, k, row, col)


lsci.AddTensorMatrix.argtypes = [ctypes.POINTER(tensor),
                                 ctypes.c_size_t,
                                 ctypes.c_size_t]
lsci.AddTensorMatrix.restype = None

def AddTensorMatrix(t, row, col):
    lsci.AddTensorMatrix(t, row, col)


lsci.NewTensor.argtypes = [ctypes.POINTER(ctypes.POINTER(tensor)),
                           ctypes.c_size_t]
lsci.NewTensor.restype = None

def NewTensor(a_):
    """
    NewTensor: Allocate in memory a libscientific tensor
    from a list of lists of lists
    """
    a = None
    if "numpy" in str(type(a_)):
        a = a_.tolist()
    else:
        a = a_
    order = len(a)
    t = ctypes.POINTER(tensor)()
    lsci.NewTensor(ctypes.pointer(t), order)
    for k in range(order):
        nrows = len(a[k])
        ncols = len(a[k][0])
        lsci.NewTensorMatrix(t, k, nrows, ncols)
        for i in range(nrows):
            for j in range(ncols):
                val = None
                try:
                    val = float(a[k][i][j])
                except ValueError:
                    val = None

                if val is None:
                    t.contents.m[k].contents.data[i][j] = misc.missing_value()
                else:
                    t.contents.m[k].contents.data[i][j] = val
    return t


lsci.DelTensor.argtypes = [ctypes.POINTER(ctypes.POINTER(tensor))]
lsci.DelTensor.restype = None

def DelTensor(t):
    lsci.DelTensor(ctypes.pointer(t))


lsci.PrintTensor.argtypes = [ctypes.POINTER(tensor)]
lsci.PrintTensor.restype = None

def PrintTensor(t):
    lsci.PrintTensor(t)


lsci.setTensorValue.argtypes = [ctypes.POINTER(tensor)]
lsci.setTensorValue.restype = None

def setTensorValue(t, k, i, j):
    lsci.setTensorValue(t, k, i, j, val)


lsci.TensorAppendColumn.argtypes = [ctypes.POINTER(tensor),
                                    ctypes.c_size_t,
                                    ctypes.POINTER(vect.dvector)]
lsci.TensorAppendColumn.restype = None

def TensorAppendColumn(t, k, cvect):
    lsci.TensorAppendColumn(t, k, dvect)


lsci.TensorAppendRow.argtypes = [ctypes.POINTER(tensor),
                                 ctypes.c_size_t,
                                 ctypes.POINTER(vect.dvector)]
lsci.TensorAppendRow.restype = None

def TensorAppendRow(t, k, rvect):
    lsci.TensorAppendRow(t, k, rvect)


lsci.TensorSet.argtypes = [ctypes.POINTER(tensor),
                           ctypes.c_double]
lsci.TensorSet.restype = None

def TensorSet(t, val):
    lsci.TensorSet(t, val)


lsci.TensorCopy.argtypes = [ctypes.POINTER(tensor),
                           ctypes.POINTER(ctypes.POINTER(tensor))]
lsci.TensorCopy.restype = None

def TensorCopy(src, dst):
    lsci.TensorCopy(src, ctypes.pointer(dst))


def TensorToList(t):
    lst = []
    for k in range(t.contents.order):
        lst.append(mx.MatrixToList(t.contents.m[k]))
    return lst


class Tensor(object):
    """
    Translate a list of list of list into a libscientific tensor
    """
    def __init__(self, t_):
        self.t = NewTensor(t_)

    def __del__(self):
        DelTensor(self.t)
        del self.t
        self.t = None

    def __getitem__(self, keys):
        k, i, j = keys
        return self.data_ptr().m[k].contents.data[i][j]

    def __setitem__(self, keys, value):
        k, i, j = keys
        self.data_ptr().m[k].contents.data[i][j] = value
    
    def order(self):
        return self.data_ptr().order
    
    def nrow(self, k):
        return self.data_ptr().m[k].contents.row

    def ncol(self, k):
        return self.data_ptr().m[k].contents.col

    def data_ptr(self):
        return self.t.contents

    def tolist(self):
        return TensorToList(self.t)

    def fromlist(self, t_):
        DelTensor(self.t)
        del self.t
        self.t = NewTensor(t_)

    def fromnumpy(self, npt):
        t_ = npt.tolist()
        DelTensor(self.t)
        del self.t
        self.t = NewTensor(t_)

    def debug(self):
        PrintTensor(self.t)



if __name__ in "__main__":
    print("Tensort Test")
    a = [[[1,2,3],[1,2,3]],[[4,5,6],[4,5,6]]]
    t = NewTensor(a)
    PrintTensor(t)
    lst = TensorToList(t)
    print(lst)
    AddTensorMatrix(t, 5, 6)
    PrintTensor(t)
    b = initTensor()
    TensorCopy(t, b)
    PrintTensor(b)
    DelTensor(t)
    DelTensor(b)
    
    t = initTensor()
    AddTensorMatrix(t, 3, 2)
    PrintTensor(t)
    DelTensor(t)
    
    t = Tensor(a)
    print(type(t))
    print(t.order())
    print(t.nrow(0))
    print(t.ncol(0))
    print(t.tolist())
    a[0][0][0] = -9876
    t.fromlist(a)
    t.debug()
    
    import numpy as np
    a = np.array(a)
    a[0][1][1] = 765
    t.fromnumpy(a)
    t.debug()
