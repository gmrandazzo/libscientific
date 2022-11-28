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
from libscientific.loadlibrary import load_libscientific_library
from libscientific import misc

lsci = load_libscientific_library()


class STRVECTOR(ctypes.Structure):
    """
    string vector data structure
    """
    _fields_ = [
        ("data", ctypes.POINTER(ctypes.POINTER(ctypes.c_char))),
        ("size", ctypes.c_size_t)]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


class DVECTOR(ctypes.Structure):
    """
    double vector data structure
    """
    _fields_ = [
        ("data", ctypes.POINTER(ctypes.c_double)),
        ("size", ctypes.c_size_t)]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


lsci.initDVector.argtypes = [ctypes.POINTER(ctypes.POINTER(DVECTOR))]
lsci.initDVector.restype = None


def init_dvector():
    """
    initDVector: Allocate in memory an empty libscientific double vector
    """
    d_vect = ctypes.POINTER(DVECTOR)()
    lsci.initDVector(ctypes.pointer(d_vect))
    return d_vect


lsci.NewDVector.argtypes = [ctypes.POINTER(ctypes.POINTER(DVECTOR)),
                            ctypes.c_size_t]
lsci.NewDVector.restype = None


def new_dvector(v_lst):
    """
    NewDVector: Allocate in memory a libscientific double vector from a list
    """
    size = len(v_lst)
    d_vect = ctypes.POINTER(DVECTOR)()
    lsci.NewDVector(ctypes.pointer(d_vect), size)

    for i in range(size):
        val = None
        try:
            val = float(v_lst[i])
        except ValueError:
            val = None
        if val is None:
            lsci.setDVectorValue(d_vect, i, misc.missing_value())
        else:
            lsci.setDVectorValue(d_vect, i, val)
    return d_vect

lsci.DelDVector.argtypes = [ctypes.POINTER(ctypes.POINTER(DVECTOR))]
lsci.DelDVector.restype = None

def del_dvector(dvect):
    """
    DelDVector: Delete an allocated libscientific double vector
    """
    lsci.initDVector(ctypes.pointer(dvect))

lsci.DVectorResize.argtypes = [ctypes.POINTER(DVECTOR),
                               ctypes.c_size_t]
lsci.DVectorResize.restype = None


def dvector_resize(dvect, size):
    """
    DVectorResize: Resize an already allocated or reallocate a libscientific
                   double vector
    """
    lsci.DVectorResize(dvect, size)


lsci.PrintDVector.argtypes = [ctypes.POINTER(DVECTOR)]
lsci.PrintDVector.restype = None


def print_dvector(dvect):
    """
    PrintDVector: Print to video a libscientific double vector
    """
    lsci.PrintDVector(dvect)


lsci.DVectorHasValue.argtypes = [ctypes.POINTER(DVECTOR), ctypes.c_double]
lsci.DVectorHasValue.restype = ctypes.c_int


def dvector_has_value(dvect, val):
    """
    DVectorHasValue: Check if a libscientific double vector contains an exact
                     value "val" and return 0 or 1 respectivelly for have or not have.
    """
    return lsci.DVectorHasValue(dvect, val)


lsci.DVectorSet.argtypes = [ctypes.POINTER(DVECTOR), ctypes.c_double]
lsci.DVectorSet.restype = None


def dvector_set(dvect, val):
    """
    DVectorSet: Set all values of a libscientific double vector to "val"
    """
    lsci.DVectorSet(dvect, val)


lsci.setDVectorValue.argtypes = [ctypes.POINTER(DVECTOR),
                                 ctypes.c_size_t,
                                 ctypes.c_double]
lsci.setDVectorValue.restype = None


def set_dvector_value(dvect, indx, val):
    """
    setDVectorValue: Set/modify a value in the row of a libscientific
                     double vector
    """
    lsci.setDVectorValue(dvect, indx, val)


lsci.getDVectorValue.argtypes = [ctypes.POINTER(DVECTOR),
                                 ctypes.c_size_t]
lsci.getDVectorValue.restype = ctypes.c_double


def get_dvector_value(dvect, indx):
    """
    getDVectorValue: Get a value in the row of a libscientific
                     double vector
    """
    return lsci.getDVectorValue(dvect, indx)


def dvector_tolist(dvect):
    """
    Convert a libscientific double vector into a python list
    """
    d_lst = []
    try:
        for i in range(dvect[0].size):
            d_lst.append(dvect[0].data[i])
    except TypeError:
        for i in range(dvect.size):
            d_lst.append(dvect.data[i])
    return d_lst


lsci.DVectorAppend.argtypes = [ctypes.POINTER(DVECTOR),
                               ctypes.c_double]
lsci.DVectorAppend.restype = None

def dvector_append(dvect, val):
    """
    Append a value to a double vector d
    """
    return lsci.DVectorAppend(dvect, val)


lsci.DVectorRemoveAt.argtypes = [ctypes.POINTER(DVECTOR),
                                 ctypes.c_size_t]
lsci.DVectorRemoveAt.restype = None

def dvector_remove_at(dvect, indx):
    """
    Remove a value from a double vector d at index indx
    """
    return lsci.DVectorRemoveAt(dvect, indx)

lsci.DVectorCopy.argtypes = [ctypes.POINTER(DVECTOR),
                             ctypes.POINTER(DVECTOR)]
lsci.DVectorCopy.restype = None

def dvector_copy(dvect_src):
    """
    Create a copy of dvector d to a
    """
    dvect_dst = init_dvector()
    lsci.DVectorCopy(dvect_src, dvect_dst)
    return dvect_dst


class DVector():
    """
    Translate a list  into a libscientific double vector
    """
    def __init__(self, dvect_):
        self.dvect = new_dvector(dvect_)

    def __del__(self):
        del_dvector(self.dvect)
        del self.dvect

    def __getitem__(self, key):
        return self.data_ptr()[key]

    def __setitem__(self, key, value):
        set_dvector_value(self.dvect, key, value)

    def size(self):
        """
        return the size of the dvector
        """
        return self.dvect[0].size

    def data_ptr(self):
        """
        return the pointer to data
        """
        return self.dvect[0].data

    def append(self, value):
        """
        Append a value to the dvector
        """
        return dvector_append(self.dvect, value)

    def extend(self, v_lst):
        """
        Extend the dvector by adding a list
        """
        for val in v_lst:
            dvector_append(self.dvect, val)

    def tolist(self):
        """
        Convert the dvector to a list
        """
        return dvector_tolist(self.dvect)

    def fromlist(self, vlst_):
        """
        Convert a list to a dvector
        """
        del_dvector(self.dvect)
        del self.dvect
        self.dvect = new_dvector(vlst_)

    def debug(self):
        """
        Debug the dvector
        """
        print_dvector(self.dvect)


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
    dvector_append(d.dvect, -123)
    d.debug()

    print("Append in a different way -123 at the end")
    d.append(-123)
    d.debug()

    print("remove at index 1, then value -2")
    dvector_remove_at(d.dvect, 1)
    d.debug()

    print("Extend d with b")
    print("b:")
    b = [random() for j in range(4)]
    print(b)
    d.extend(b)
    print("d extended:")
    d.debug()

    print("Create a copy of d.d in q")
    q = dvector_copy(d.dvect)
    print_dvector(q)
    del_dvector(q)

    print("Check if the double vector d have the value -123")
    print(dvector_has_value(d.dvect, -13.0000))
