"""
vector libscientific python binding

Copyright (C) <2022>  Giuseppe Marco Randazzo

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
from libscientific import vector as vect
from libscientific.loadlibrary import load_libscientific_library

lsci = load_libscientific_library()


class DVECTLIST(ctypes.Structure):
    """
    DVECTLIST data structure
    """
    _fields_ = [
        ("dvector", ctypes.POINTER(ctypes.POINTER(vect.DVECTOR))),
        ("size", ctypes.c_size_t)]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


lsci.initDVectorList.argtypes = [ctypes.POINTER(ctypes.POINTER(DVECTLIST))]
lsci.initDVectorList.restype = None


def init_dvector_list():
    """
    init_dvector_list: Allocate in memory an empty libscientific double vector list
    """
    dvl = ctypes.POINTER((DVECTLIST))()
    lsci.initDVectorList(ctypes.pointer(dvl))
    return dvl


def new_dvector_list(dvl_):
    """
    new_dvector_list: Allocate in memory a libscientific double vector list
    from a list of lists
    """
    dvl = init_dvector_list()
    for v_lst in dvl_:
        dvector_list_append(dvl, v_lst)
    return dvl


lsci.DelDVectorList.argtypes = [ctypes.POINTER(ctypes.POINTER(DVECTLIST))]
lsci.DelDVectorList.restype = None


def del_dvector_list(dvl):
    """
    del_dvector_list: Delete an allocated libscientific double vector
    """
    lsci.DelDVectorList(ctypes.pointer(dvl))


lsci.DVectorListAppend.argtypes = [ctypes.POINTER(DVECTLIST),
                                   ctypes.POINTER(vect.DVECTOR)]
lsci.DVectorListAppend.restype = None

def dvector_list_append(dvl, dvect):
    """
    dvector_list_append: Append a double vector to a DvectorList
    """
    if isinstance(dvect, list):
        dvect_ = vect.new_dvector(dvect)
        lsci.DVectorListAppend(dvl, dvect_)
        vect.del_dvector(dvect_)
        del dvect_
    else:
        lsci.DVectorListAppend(dvl, dvect)
    return 0


def dvector_list_tolist(dvl):
    """
    convert a divector list into a list of lists
    """
    lsts = []
    for i in range(dvl[0].size):
        lsts.append(vect.dvector_tolist(dvl[0].dvector[i].contents))
    return lsts


class DVectorList():
    """
    Translate a list  into a libscientific double vector
    """
    def __init__(self, dvl_):
        self.dvl = new_dvector_list(dvl_)

    def __del__(self):
        del_dvector_list(self.dvl)
        del self.dvl

    def __getitem__(self, key):
        return self.data_ptr()[key].contents

    def __setitem__(self, key, v_lst):
        for i, val in enumerate(v_lst):
            if i < self.data_ptr()[key].contents.size:
                self.data_ptr()[key].contents.data[i] = val
            else:
                break
        return 0

    def size(self):
        """
        return the size of the divector
        """
        return self.dvl[0].size

    def data_ptr(self):
        """
        return the pointer to data
        """
        return self.dvl[0].dvector

    def append(self, v_lst):
        """
        Append a value to the dvector
        """
        if isinstance(lst, vect.DVECTOR):
            dvector_list_append(self.dvl, v_lst)
        else:
            dv_lst = vect.new_dvector(lst)
            dvector_list_append(self.dvl, dv_lst)
            vect.del_dvector(dv_lst)
            del dv_lst
        return 0

    def tolist(self):
        """
        Convert the dvector list to a list of list
        """
        return dvector_list_tolist(self.dvl)

    def fromlist(self, v_lists):
        """
        Convert a list of list to a dvector list
        """
        for v_lst in v_lists:
            dvector_list_append(self.dvl, v_lst)
        return 0

    def debug(self):
        """
        Debug the double list vector
        """
        for i in range(self.dvl[0].size):
            vect.print_dvector(self.dvl[0].dvector[i].contents)


if __name__ in "__main__":
    from random import random
    a = [[random() for j in range(10)] for i in range(3)]
    d = DVectorList(a)
    d.debug()
    print("get value")
    vect.print_dvector(d[1])
    print("set list")
    d[1] = [8, 8, 8, 8, 8, 8, 8, 8, 8, 8]
    dlst = d.tolist()
    print("print the list converted")
    for item in dlst:
        print(item)

    print("Append a list - high level")
    lst = [9, 8, 7, 6]
    d.append(lst)
    print("Reappend ppend a list - low level")
    dvector_list_append(d.dvl, lst)
    d.debug()
