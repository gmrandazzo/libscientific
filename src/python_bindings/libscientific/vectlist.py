# vectlist.py libscientific python binding
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
"""Lists of double vectors stored in the C library.

Used by the models for per-block or per-layer results (column averages,
explained variances, ...).

    >>> from libscientific.vectlist import DVectorList
    >>> DVectorList([[1.0, 2.0], [3.0]]).tolist()
    [[1.0, 2.0], [3.0]]
"""

import ctypes
from typing import Any

from libscientific import vector as vect
from libscientific.loadlibrary import load_libscientific_library

lsci = load_libscientific_library()


class DVECTLIST(ctypes.Structure):
    """C structure of a list of double vectors (low-level API).

    Attributes
    ----------
    dvector : POINTER(POINTER(DVECTOR))
        The vectors; ``dvector[i]`` is vector ``i``. The vectors can have
        different sizes.
    size : int
        Number of vectors.
    """

    _fields_ = [
        ("dvector", ctypes.POINTER(ctypes.POINTER(vect.DVECTOR))),
        ("size", ctypes.c_size_t),
    ]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


lsci.initDVectorList.argtypes = [ctypes.POINTER(ctypes.POINTER(DVECTLIST))]
lsci.initDVectorList.restype = None


def init_dvector_list():
    """Allocate an empty list of double vectors.

    Returns
    -------
    ctypes.POINTER(DVECTLIST)
        A new list; release it with :func:`del_dvector_list`.
    """
    dvl = ctypes.POINTER(DVECTLIST)()
    lsci.initDVectorList(ctypes.pointer(dvl))
    return dvl


def new_dvector_list(dvl_):
    """Create a list of double vectors from a list of lists.

    Parameters
    ----------
    dvl_ : list[list[float]]
        The vectors (they can have different lengths).

    Returns
    -------
    ctypes.POINTER(DVECTLIST)
        A new list; release it with :func:`del_dvector_list`.

    Examples
    --------
    >>> dvl = new_dvector_list([[1.0, 2.0], [3.0]])
    >>> dvector_list_tolist(dvl)
    [[1.0, 2.0], [3.0]]
    >>> del_dvector_list(dvl)
    """
    dvl = init_dvector_list()
    for v_lst in dvl_:
        dvector_list_append(dvl, v_lst)
    return dvl


lsci.DelDVectorList.argtypes = [ctypes.POINTER(ctypes.POINTER(DVECTLIST))]
lsci.DelDVectorList.restype = None


def del_dvector_list(dvl):
    """Free a list of double vectors and all its vectors.

    Parameters
    ----------
    dvl : ctypes.POINTER(DVECTLIST)
        The list to release. It must not be used afterwards.
    """
    lsci.DelDVectorList(ctypes.pointer(dvl))


lsci.DVectorListAppend.argtypes = [ctypes.POINTER(DVECTLIST), ctypes.POINTER(vect.DVECTOR)]
lsci.DVectorListAppend.restype = None


def dvector_list_append(dvl, dvect):
    """Append a copy of a vector to the list.

    Parameters
    ----------
    dvl : ctypes.POINTER(DVECTLIST)
        The list; its size grows by one.
    dvect : list[float], DVector or POINTER(DVECTOR)
        The vector to append.
    """
    if isinstance(dvect, vect.DVector):
        lsci.DVectorListAppend(dvl, dvect.dvect)
    elif isinstance(dvect, list):
        dvect_ = vect.new_dvector(dvect)
        lsci.DVectorListAppend(dvl, dvect_)
        vect.del_dvector(dvect_)
        del dvect_
    else:
        lsci.DVectorListAppend(dvl, dvect)
    return 0


def dvector_list_tolist(dvl: Any) -> list[list[float]]:
    """Convert a list of double vectors into Python lists.

    Parameters
    ----------
    dvl : ctypes.POINTER(DVECTLIST)
        The list.

    Returns
    -------
    list[list[float]]
    """
    lsts: list[list[float]] = []
    for i in range(dvl[0].size):
        lsts.append(vect.dvector_tolist(dvl[0].dvector[i].contents))
    return lsts


class DVectorList:
    """A list of double vectors stored in the C library.

    Parameters
    ----------
    dvl_ : list[list[float]] or None
        Initial vectors (copied), or ``None`` for an empty list.

    Attributes
    ----------
    dvl : ctypes.POINTER(DVECTLIST)
        The underlying C list, owned by this object.

    Examples
    --------
    >>> dvl = DVectorList([[1.0, 2.0], [3.0, 4.0, 5.0]])
    >>> dvl.append([6.0])
    >>> len(dvl), dvl.tolist()
    (3, [[1.0, 2.0], [3.0, 4.0, 5.0], [6.0]])
    """

    def __init__(self, dvl_=None):
        """Create the list (see the class documentation)."""
        if dvl_ is None:
            self.dvl = init_dvector_list()
        else:
            self.dvl = new_dvector_list(dvl_)

    def __del__(self):
        """Release the C list."""
        if getattr(self, "dvl", None) is not None:
            del_dvector_list(self.dvl)
            del self.dvl
            self.dvl = None

    def __repr__(self):
        return f"DVectorList({self.tolist()})"

    def __str__(self):
        return f"DVectorList({self.tolist()})"

    def __len__(self):
        """Return the number of vectors."""
        return self.size()

    def __iter__(self):
        """Iterate over the vectors, each as a Python list."""
        return iter(self.tolist())

    def __getitem__(self, key):
        """Return vector ``key`` as a C ``DVECTOR`` structure.

        Use ``dvl.tolist()[key]`` to get the values as a Python list.

        Raises
        ------
        IndexError
            If the index is out of range.
        """
        return self.data_ptr()[self._check_index(key)].contents

    def _check_index(self, key: int) -> int:
        """Validate an index; negative values count from the end."""
        size = self.size()
        idx = key + size if key < 0 else key
        if not 0 <= idx < size:
            raise IndexError(f"index {key} out of range for a list of {size} vectors")
        return idx

    def __setitem__(self, key, v_lst):
        """Overwrite the values of vector ``key`` (``dvl[key] = values``).

        Only the first ``len(vector)`` values are used; the vector size does not
        change.

        Raises
        ------
        IndexError
            If the index is out of range.
        """
        key = self._check_index(key)
        for i, val in enumerate(v_lst):
            if i < self.data_ptr()[key].contents.size:
                self.data_ptr()[key].contents.data[i] = val
            else:
                break

    def size(self):
        """Return the number of vectors.

        Returns
        -------
        int
        """
        return self.dvl[0].size

    def data_ptr(self):
        """Return the C array of vector pointers (low-level access).

        Returns
        -------
        POINTER(POINTER(DVECTOR))
        """
        return self.dvl[0].dvector

    def append(self, v_lst):
        """Append a vector.

        Parameters
        ----------
        v_lst : list[float]
            The values of the new vector.
        """
        dvector_list_append(self.dvl, v_lst)

    def tolist(self):
        """Return all the vectors as Python lists.

        Returns
        -------
        list[list[float]]
        """
        return dvector_list_tolist(self.dvl)

    def fromlist(self, v_lists):
        """Replace the content with a list of lists.

        Parameters
        ----------
        v_lists : list[list[float]]
            The new vectors.
        """
        for v_lst in v_lists:
            dvector_list_append(self.dvl, v_lst)
        return 0

    def debug(self):
        """Print the list to standard output."""
        for i in range(self.dvl[0].size):
            vect.print_dvector(self.dvl[0].dvector[i].contents)
