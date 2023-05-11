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


class UIVECTOR(ctypes.Structure):
    """
    unsigned int vector data structure
    """
    _fields_ = [
        ("data", ctypes.POINTER(ctypes.c_size_t)),
        ("size", ctypes.c_size_t)]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


class IVECTOR(ctypes.Structure):
    """
    unsigned int vector data structure
    """
    _fields_ = [
        ("data", ctypes.POINTER(ctypes.c_int)),
        ("size", ctypes.c_size_t)]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


def xvector_tolist(dvect):
    """
    Convert a libscientific double/unsigned int/int vector into a python list
    """
    d_lst = []
    try:
        for i in range(dvect[0].size):
            d_lst.append(dvect[0].data[i])
    except TypeError:
        for i in range(dvect.size):
            d_lst.append(dvect.data[i])
    return d_lst

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
    lsci.DelDVector(ctypes.pointer(dvect))

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


lsci.DVectorAppend.argtypes = [ctypes.POINTER(DVECTOR),
                               ctypes.c_double]
lsci.DVectorAppend.restype = None

def dvector_append(dvect, val):
    """
    Append a value to a double vector d
    """
    return lsci.DVectorAppend(dvect, val)


def dvector_tolist(dvect):
    """
    Double vector to python list
    """
    return xvector_tolist(dvect)


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


lsci.initUIVector.argtypes = [ctypes.POINTER(ctypes.POINTER(UIVECTOR))]
lsci.initUIVector.restype = None


def init_uivector():
    """
    initUIVector: Allocate in memory an empty libscientific double vector
    """
    uivect = ctypes.POINTER(UIVECTOR)()
    lsci.initUIVector(ctypes.pointer(uivect))
    return uivect


lsci.NewUIVector.argtypes = [ctypes.POINTER(ctypes.POINTER(UIVECTOR)),
                              ctypes.c_size_t]
lsci.NewUIVector.restype = None


def new_uivector(v_lst):
    """
    NewUIVector: Allocate in memory a libscientific unsigned int vector from a list
    """
    size = len(v_lst)
    uivect = ctypes.POINTER(UIVECTOR)()
    lsci.NewUIVector(ctypes.pointer(uivect), size)

    for i in range(size):
        val = None
        try:
            val = int(v_lst[i])
        except ValueError:
            val = None
        if val is None:
            lsci.setUIVectorValue(uivect, i, misc.missing_value())
        else:
            lsci.setUIVectorValue(uivect, i, val)
    return uivect

lsci.DelUIVector.argtypes = [ctypes.POINTER(ctypes.POINTER(UIVECTOR))]
lsci.DelUIVector.restype = None

def del_uivector(uivect):
    """
    DelUIVector: Delete an allocated libscientific unsigned int vector
    """
    lsci.DelUIVector(ctypes.pointer(uivect))

lsci.UIVectorResize.argtypes = [ctypes.POINTER(UIVECTOR),
                                ctypes.c_size_t]
lsci.UIVectorResize.restype = None


def uivector_resize(uivect, size):
    """
    DVectorResize: Resize an already allocated or reallocate a libscientific
                   unsigned int vector
    """
    lsci.UIVectorResize(uivect, size)


lsci.PrintUIVector.argtypes = [ctypes.POINTER(UIVECTOR)]
lsci.PrintUIVector.restype = None


def print_uivector(uivect):
    """
    PrintUIVector: Print to video a libscientific unsigned int vector
    """
    lsci.PrintUIVector(uivect)


lsci.setUIVectorValue.argtypes = [ctypes.POINTER(UIVECTOR),
                                  ctypes.c_size_t,
                                  ctypes.c_size_t]
lsci.setUIVectorValue.restype = None


def set_uivector_value(uivect, indx, val):
    """
    setUIVectorValue: Set/modify a value in the row of a libscientific
                     unsigned int vector
    """
    lsci.setUIVectorValue(uivect, indx, val)


lsci.getUIVectorValue.argtypes = [ctypes.POINTER(UIVECTOR),
                                  ctypes.c_size_t]
lsci.getUIVectorValue.restype = ctypes.c_size_t


def get_uivector_value(uivect, indx):
    """
    getUIVectorValue: Get a value in the row of a libscientific
                     double vector
    """
    return lsci.getUIVectorValue(uivect, indx)

lsci.UIVectorAppend.argtypes = [ctypes.POINTER(UIVECTOR),
                                ctypes.c_size_t]
lsci.UIVectorAppend.restype = None

def uivector_append(uivect, val):
    """
    Append a value to an unsigned int vector
    """
    return lsci.UIVectorAppend(uivect, val)

def uivector_tolist(uivect):
    """
    Unsigned int vector to python list
    """
    return xvector_tolist(uivect)

lsci.UIVectorRemoveAt.argtypes = [ctypes.POINTER(UIVECTOR),
                                  ctypes.c_size_t]
lsci.UIVectorRemoveAt.restype = None

def uivector_remove_at(uivect, indx):
    """
    Remove a value from an unsigned int vector at index indx
    """
    return lsci.UIVectorRemoveAt(uivect, indx)

lsci.initIVector.argtypes = [ctypes.POINTER(ctypes.POINTER(IVECTOR))]
lsci.initIVector.restype = None


def init_ivector():
    """
    initIVector: Allocate in memory an empty libscientific double vector
    """
    ivect = ctypes.POINTER(IVECTOR)()
    lsci.initIVector(ctypes.pointer(ivect))
    return ivect


lsci.NewIVector.argtypes = [ctypes.POINTER(ctypes.POINTER(IVECTOR)),
                              ctypes.c_size_t]
lsci.NewIVector.restype = None


def new_ivector(v_lst):
    """
    NewIVector: Allocate in memory a libscientific int vector from a list
    """
    size = len(v_lst)
    ivect = ctypes.POINTER(DVECTOR)()
    lsci.NewUIDVector(ctypes.pointer(ivect), size)

    for i in range(size):
        val = None
        try:
            val = int(v_lst[i])
        except ValueError:
            val = None
        if val is None:
            lsci.setIVectorValue(ivect, i, misc.missing_value())
        else:
            lsci.setIVectorValue(ivect, i, val)
    return ivect


lsci.DelIVector.argtypes = [ctypes.POINTER(ctypes.POINTER(IVECTOR))]
lsci.DelIVector.restype = None


def del_ivector(ivect):
    """
    DelIVector: Delete an allocated libscientific int vector
    """
    lsci.DelIVector(ctypes.pointer(ivect))


#lsci.PrintIVector.argtypes = [ctypes.POINTER(IVECTOR)]
#lsci.PrintIVector.restype = None


#def print_ivector(dvect):
#    """
#    PrintIVector: Print to video a libscientific int vector
#    """
#    lsci.PrintIVector(dvect)


lsci.setIVectorValue.argtypes = [ctypes.POINTER(IVECTOR),
                                  ctypes.c_size_t,
                                  ctypes.c_size_t]
lsci.setIVectorValue.restype = None


def set_ivector_value(ivect, indx, val):
    """
    setIVectorValue: Set/modify a value in the row of a libscientific
                     int vector
    """
    lsci.setIVectorValue(ivect, indx, val)


lsci.getIVectorValue.argtypes = [ctypes.POINTER(IVECTOR),
                                  ctypes.c_size_t]
lsci.getIVectorValue.restype = ctypes.c_size_t


def get_ivector_value(ivect, indx):
    """
    getIVectorValue: Get a value in the row of a libscientific
                     double vector
    """
    return lsci.getIVectorValue(ivect, indx)


def ivector_tolist(ivect):
    """
    Int vector to python list
    """
    return xvector_tolist(ivect)

lsci.IVectorAppend.argtypes = [ctypes.POINTER(IVECTOR),
                                ctypes.c_size_t]
lsci.IVectorAppend.restype = None

def ivector_append(ivect, val):
    """
    Append a value to an int vector
    """
    return lsci.IVectorAppend(ivect, val)


lsci.IVectorRemoveAt.argtypes = [ctypes.POINTER(IVECTOR),
                                  ctypes.c_size_t]
lsci.IVectorRemoveAt.restype = None

def ivector_remove_at(ivect, indx):
    """
    Remove a value from an int vector at index indx
    """
    return lsci.IVectorRemoveAt(ivect, indx)


class DVector():
    """
    Translate a list  into a libscientific double vector
    """
    def __init__(self, dvect_ = None):
        if dvect_ is None:
            self.dvect = init_dvector()
        else:
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
        return xvector_tolist(self.dvect)

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


class UIVector():
    """
    Translate a list  into a libscientific unsigned int vector
    """
    def __init__(self, uivect_):
        self.uivect = new_uivector(uivect_)

    def __del__(self):
        del_uivector(self.uivect)
        del self.uivect

    def __getitem__(self, key):
        return self.data_ptr()[key]

    def __setitem__(self, key, value):
        set_uivector_value(self.uivect, key, value)

    def size(self):
        """
        return the size of the uivector
        """
        return self.uivect[0].size

    def data_ptr(self):
        """
        return the pointer to data
        """
        return self.uivect[0].data

    def append(self, value):
        """
        Append a value to the uivector
        """
        return uivector_append(self.uivect, value)

    def extend(self, v_lst):
        """
        Extend the dvector by adding a list
        """
        for val in v_lst:
            uivector_append(self.uivect, val)

    def tolist(self):
        """
        Convert the dvector to a list
        """
        return xvector_tolist(self.uivect)

    def fromlist(self, vlst_):
        """
        Convert a list to a dvector
        """
        del_uivector(self.uivect)
        del self.uivect
        self.uivect = new_dvector(vlst_)

    def debug(self):
        """
        Debug the iuvector
        """
        print_uivector(self.uivect)


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
    del d

    print("Start test on uivector")
    a = [0,1,2,3,4,5,6]
    u = UIVector(a)
    u.debug()
    print("get value")
    print(u[1])
    print("set value 2 in poistion 1")
    u[1] = 2
    lst = u.tolist()
    print("print the list converted from the unsigned int vector")
    for item in lst:
        print(item)

    print("Add at the end the value 123")
    uivector_append(u.uivect, 123)
    u.debug()

    print("Append 123 at the end in a different way ")
    u.append(123)
    u.debug()

    print("remove at index 1")
    uivector_remove_at(u.uivect, 1)
    u.debug()

    print("Extend d with b")
    print("b:")
    b = [7, 8, 9, 10]
    print(b)
    u.extend(b)
    print("u extended:")
    u.debug()
