"""
vector libscientific python binding

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

lsci = LoadLibrary()


class strvector(ctypes.Structure):
    """
    string vector class
    """
    _fields_ = [
        ("data",    ctypes.POINTER(ctypes.POINTER(ctypes.c_char))),
        ("size",     ctypes.c_size_t)]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


class dvector(ctypes.Structure):
    """
    double vector class
    """
    _fields_ = [
        ("data",    ctypes.POINTER(ctypes.c_double)),
        ("size",     ctypes.c_size_t)]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


lsci.initDVector.argtypes = [ctypes.POINTER(ctypes.POINTER(dvector))]
lsci.initDVector.restype = None


def initDVector():
    """
    initDVector: Allocate in memory an empty libscientific double vector
    """
    d = ctypes.POINTER(dvector)()
    lsci.initDVector(ctypes.pointer(d))
    return d


lsci.NewDVector.argtypes = [ctypes.POINTER(ctypes.POINTER(dvector)),
                            ctypes.c_size_t]
lsci.NewDVector.restype = None


def NewDVector(vlst):
    """
    NewDVector: Allocate in memory a libscientific double vector from a list
    """
    size = len(vlst)
    d = ctypes.POINTER(dvector)()
    lsci.NewDVector(ctypes.pointer(d), size)

    for i in range(size):
        val = None
        try:
            val = float(vlst[i])
        except ValueError:
            val = None
        if val is None:
            lsci.setDVectorValue(d, i, misc.missing_value())
        else:
            lsci.setDVectorValue(d, i, val)
    return d


lsci.DVectorResize.argtypes = [ctypes.POINTER(dvector),
                               ctypes.c_size_t]
lsci.DVectorResize.restype = None


def DVectorResize(d, size):
    """
    DVectorResize: Resize an already allocated or reallocate a libscientific
                   double vector
    """
    lsci.DVectorResize(d, size)


lsci.DelDVector.argtypes = [ctypes.POINTER(ctypes.POINTER(dvector))]
lsci.DelDVector.restype = None


def DelDVector(d):
    """
    DelDVector: Delete an allocated libscientific double vector
    """
    lsci.DelDVector(ctypes.pointer(d))


lsci.PrintDVector.argtypes = [ctypes.POINTER(dvector)]
lsci.PrintDVector.restype = None


def PrintDVector(d):
    """
    PrintDVector: Print to video a libscientific double vector
    """
    lsci.PrintDVector(d)


lsci.DVectorHasValue.argtypes = [ctypes.POINTER(dvector), ctypes.c_double]
lsci.DVectorHasValue.restype = ctypes.c_int


def DVectorHasValue(v, val):
    """
    DVectorHasValue: Check if a libscientific double vector contains an exact
                     value "val" and return 0 or 1 respectivelly for have or not have.
    """
    return lsci.DVectorHasValue(v, val)


lsci.DVectorSet.argtypes = [ctypes.POINTER(dvector), ctypes.c_double]
lsci.DVectorSet.restype = None


def DVectorSet(d, val):
    """
    DVectorSet: Set all values of a libscientific double vector to "val"
    """
    lsci.DVectorSet(d, val)


lsci.setDVectorValue.argtypes = [ctypes.POINTER(dvector),
                                 ctypes.c_size_t,
                                 ctypes.c_double]
lsci.setDVectorValue.restype = None


def setDVectorValue(d, row, value):
    """
    setDVectorValue: Set/modify a value in the row of a libscientific
                     double vector
    """
    lsci.setDVectorValue(d, row, value)


lsci.getDVectorValue.argtypes = [ctypes.POINTER(dvector),
                                 ctypes.c_size_t]
lsci.getDVectorValue.restype = ctypes.c_double


def getDVectorValue(d, row):
    """
    getDVectorValue: Get a value in the row of a libscientific
                     double vector
    """
    return lsci.getDVectorValue(d, row)


def DVectorToList(d):
    dlst = []
    try:
        for i in range(d[0].size):
            dlst.append(d[0].data[i])
    except TypeError:
        for i in range(d.size):
            dlst.append(d.data[i])
    return dlst


lsci.DVectorAppend.argtypes = [ctypes.POINTER(dvector),
                               ctypes.c_double]
lsci.DVectorAppend.restype = None

def DVectorAppend(d, val):
    """
    Append a value to a double vector d
    """
    return lsci.DVectorAppend(d, val)


lsci.DVectorRemoveAt.argtypes = [ctypes.POINTER(dvector),
                                 ctypes.c_size_t]
lsci.DVectorRemoveAt.restype = None

def DVectorRemoveAt(d, indx):
    """
    Remove a value from a double vector d at index indx
    """
    return lsci.DVectorRemoveAt(d, indx)

lsci.DVectorCopy.argtypes = [ctypes.POINTER(dvector),
                             ctypes.POINTER(dvector)]
lsci.DVectorCopy.restype = None

def DVectorCopy(src):
    """
    Create a copy of dvector d to a
    """
    dst = initDVector()
    lsci.DVectorCopy(src, dst)
    return dst


"""
TODO: Implement this functions

/* Append to a dvector an other dvector */
dvector *DVectorExtend(dvector *d1, dvector *d2);

/*Vector operations*/
double DVectorDVectorDotProd(dvector *v1, dvector *v2); /* product between two vector */

double DvectorModule(dvector *v); /* get the Dvector Module */
void DVectNorm(dvector *v, dvector *nv); /* vector normalizing */
void DVectorDVectorDiff(dvector *v1, dvector *v2, dvector **v3);
void DVectorDVectorSum(dvector *v1, dvector *v2, dvector **v3);
void DVectorMinMax(dvector *v, double *min, double *max);
void DVectorMean(dvector *d, double *mean);
void DVectorMedian(dvector *d, double *mean);
void DVectorSDEV(dvector *d, double *sdev);
void DVectorSort(dvector *v);


"""


class DVector(object):
    """
    Translate a list  into a libscientific double vector
    """
    def __init__(self, d_):
        self.d = NewDVector(d_)

    def __del__(self):
        DelDVector(self.d)
        del self.d

    def __getitem__(self, key):
        return self.data_ptr()[key]

    def __setitem__(self, key, value):
        setDVectorValue(self.d, key, value)

    def size(self):
        """
        return the size of the dvector
        """
        return self.d[0].size

    def data_ptr(self):
        """
        return the pointer to data
        """
        return self.d[0].data

    def append(self, value):
        """
        Append a value to the dvector
        """
        return DVectorAppend(self.d, value)

    def extend(self, lst):
        """
        Extend the dvector by adding a list
        """
        for item in lst:
            DVectorAppend(self.d, item)

    def tolist(self):
        """
        Convert the dvector to a list
        """
        return DVectorToList(self.d)

    def fromlist(self, vlst_):
        """
        Convert a list to a dvector
        """
        DelDVector(self.d)
        del self.d
        self.d = NewDVector(vlst_)

    def debug(self):
        """
        Debug the dvector
        """
        PrintDVector(self.d)


if __name__ in "__main__":
    from random import random
    a = [random() for j in range(10)]
    d = DVector(a)
    d.debug()
    print("get value")
    print(d[1])
    print("set value")
    d[1] = -2
    dlst = d.tolist()
    print("print the list converted from the double vector")
    for item in dlst:
        print(item)

    print("Add at the end the value -123")
    DVectorAppend(d.d, -123)
    d.debug()

    print("Append in a different way -123 at the end")
    d.append(-123)
    d.debug()

    print("remove at index 1, then value -2")
    DVectorRemoveAt(d.d, 1)
    d.debug()

    print("Extend d with b")
    print("b:")
    b = [random() for j in range(4)]
    print(b)
    d.extend(b)
    print("d extended:")
    d.debug()

    print("Create a copy of d.d in q")
    q = DVectorCopy(d.d)
    PrintDVector(q)
    DelDVector(q)

    print("Check if the double vector d have the value -123")
    print(DVectorHasValue(d.d, -13.0000))
