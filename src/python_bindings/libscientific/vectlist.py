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
from libscientific.loadlibrary import LoadLibrary
from libscientific import misc

lsci = LoadLibrary()


class dvectlist(ctypes.Structure):
    """
    dvector list class
    """
    _fields_ = [
        ("dvector", ctypes.POINTER(ctypes.POINTER(vect.dvector))),
        ("size", ctypes.c_size_t)]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


lsci.initDVectorList.argtypes = [ctypes.POINTER(ctypes.POINTER(dvectlist))]
lsci.initDVectorList.restype = None


def initDVectorList():
    """
    initDVectorList: Allocate in memory an empty libscientific double vector list
    """
    dl = ctypes.POINTER(dvectlist)()
    lsci.initDVector(ctypes.pointer(dl))
    return dl


def NewDVectorList(dl):
    """
    NewDVectorList: Allocate in memory a libscientific double vector list
    from a list of lists
    """
    size = len(dl)
    dl_ = initDVectorList()
    for lst in dl:
        DVectorListAppend(dl_, lst)
    return dl_


lsci.DelDVectorList.argtypes = [ctypes.POINTER(ctypes.POINTER(dvectlist))]
lsci.DelDVectorList.restype = None


def DelDVectorList(dl):
    """
    DelDVectorList: Delete an allocated libscientific double vector
    """
    lsci.DelDVectorList(ctypes.pointer(dl))


lsci.DVectorListAppend.argtypes = [ctypes.POINTER(dvectlist),
                                   ctypes.POINTER(vect.dvector)]
lsci.DVectorListAppend.restype = None

def DVectorListAppend(dl, dvect):
    """
    DVectorListAppend: Append a double vector to a DvectorList
    """
    if type(dvect) == list:
        lst_ = vect.NewDVector(dvect)
        lsci.DVectorListAppend(dl, lst_)
        vect.DelDVector(lst_)        
    else:
        lsci.DVectorListAppend(dl, dvect)
    return 0


def DVectorListToList(dl):
    lst = list()
    for i in range(dl[0].size):
        lst.append(vect.DVectorToList(dl[0].dvector[i].contents))
    return lst


class DVectorList(object):
    """
    Translate a list  into a libscientific double vector
    """
    def __init__(self, d_):
        self.dl = NewDVectorList(d_)

    def __del__(self):
        DelDVectorList(self.dl)
        del self.dl

    def __getitem__(self, key):
        return self.data_ptr()[key].contents
    
    def __setitem__(self, key, lst):
        for i in range(len(lst)):
                if i < self.data_ptr()[key].contents.size:
                    self.data_ptr()[key].contents.data[i] = lst[i]
                else:
                    break
        return 0

    def size(self):
        """
        return the size of the divector
        """
        return self.dl[0].size

    def data_ptr(self):
        """
        return the pointer to data
        """
        return self.dl[0].dvector

    def append(self, lst):
        """
        Append a value to the dvector
        """
        if type(lst) == vect.dvector:
            return DVectorAppend(self.dl, lst)
        else:
            lst_ = vect.NewDVector(lst)
            DVectorListAppend(self.dl, lst_)
            vect.DelDVector(lst_)
            return 0

    def tolist(self):
        """
        Convert the dvector list to a list of list
        """
        return DVectorListToList(self.dl)

    def fromlist(self, vlst_):
        """
        Convert a list of list to a dvector list
        """
        for lst in vlst_:
            DVectorListAppend(self.dl, lst)
        return 0

    def debug(self):
        """
        Debug the double list vector
        """
        for i in range(self.dl[0].size):
            vect.PrintDVector(self.dl[0].dvector[i].contents)


if __name__ in "__main__":
    from random import random
    a = [[random() for j in range(10)] for i in range(3)]
    d = DVectorList(a)
    d.debug()
    print("get value")
    vect.PrintDVector(d[1])
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
    DVectorListAppend(d.dl, lst)
    d.debug()

