# vector libscientific python binding
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
from scientific.loadlibrary import LoadLibrary

lsci = LoadLibrary()


class strvector(ctypes.Structure):
    _fields_ = [
        ("data",    ctypes.POINTER(ctypes.POINTER(ctypes.c_char))),
        ("size",     ctypes.c_size_t)]


class dvector(ctypes.Structure):
    _fields_ = [
        ("data",    ctypes.POINTER(ctypes.c_double)),
        ("size",     ctypes.c_size_t)]


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
        lsci.setDVectorValue(d, i, vlst[i])
    return d


lsci.DVectorResize.argtypes = [ctypes.POINTER(ctypes.POINTER(dvector)),
                               ctypes.c_size_t]
lsci.DVectorResize.restype = None


def DVectorResize(d, size):
    """
    DVectorResize: Resize an already allocated or reallocate a libscientific
                   double vector
    """
    lsci.DVectorResize(ctypes.pointer(d), size)


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
                     value "val" and return 1 or 0 respectivelly for yes or no.
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


"""
TODO: Implement this functions
/* Append a value to a dvector */
void DVectorAppend(dvector **d, double val);

/* Remove a value to a dvector */
void DVectorRemoveAt(dvector **d, size_t indx);

/* Copy a Dvector from dsrc: source to ddst: destination */
void DVectorCopy(dvector *dsrc, dvector **ddst);

/* Append to a dvector an other dvector */
dvector *DVectorExtend(dvector *d1, dvector *d2);

void setDVectorValue(dvector *d, size_t id, double val);
double getDVectorValue(dvector *d, size_t id);

/*check if dvector has a value val. Return 0 if is present, 1 if is not present*/
int DVectorHasValue(dvector *d, double val);

/*Vector operations*/
void DVectorSet(dvector *v, double val);
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
        return self.d[0].size

    def data_ptr(self):
        return self.d[0].data

    def tolist(self):
        return DVectorToList(self.d)

    def fromlist(self, vlst_):
        DelDVector(self.d)
        del self.d
        self.d = NewDVector(vlst_)

    def debug(self):
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
    for item in dlst:
        print(item)
