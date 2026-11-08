# vector.py libscientific python binding
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
"""Vectors stored in the C library.

- :class:`DVector` (double values) and :class:`UIVector` (unsigned
  integers): Pythonic vectors.
- Low-level functions for double (``dvector``), unsigned int (``uivector``)
  and int (``ivector``) C vectors.

    >>> from libscientific.vector import DVector
    >>> v = DVector([1.0, 2.0, 3.0])
    >>> len(v), v[-1]
    (3, 3.0)
"""

import ctypes
from typing import Any

from libscientific import misc
from libscientific.loadlibrary import load_libscientific_library

lsci = load_libscientific_library()


class STRVECTOR(ctypes.Structure):
    """C structure of a libscientific string vector (low-level API).

    Attributes
    ----------
    data : POINTER(c_char_p)
        The strings.
    size : int
        Number of strings.
    """

    _fields_ = [("data", ctypes.POINTER(ctypes.POINTER(ctypes.c_char))), ("size", ctypes.c_size_t)]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


class DVECTOR(ctypes.Structure):
    """C structure of a libscientific double vector (low-level API).

    Attributes
    ----------
    data : POINTER
        The values; ``data[i]`` is element ``i``.
    size : int
        Number of values.
    """

    _fields_ = [("data", ctypes.POINTER(ctypes.c_double)), ("size", ctypes.c_size_t)]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


class UIVECTOR(ctypes.Structure):
    """C structure of a libscientific unsigned int vector (low-level API).

    Attributes
    ----------
    data : POINTER
        The values; ``data[i]`` is element ``i``.
    size : int
        Number of values.
    """

    _fields_ = [("data", ctypes.POINTER(ctypes.c_size_t)), ("size", ctypes.c_size_t)]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


class IVECTOR(ctypes.Structure):
    """C structure of a libscientific int vector (low-level API).

    Attributes
    ----------
    data : POINTER
        The values; ``data[i]`` is element ``i``.
    size : int
        Number of values.
    """

    _fields_ = [("data", ctypes.POINTER(ctypes.c_int)), ("size", ctypes.c_size_t)]

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


def xvector_tolist(dvect: Any) -> list[Any]:
    """Convert any C vector (double, unsigned int or int) into a Python list.

    Parameters
    ----------
    dvect : POINTER(DVECTOR), POINTER(UIVECTOR) or POINTER(IVECTOR)
        The vector (a pointer or the structure itself).

    Returns
    -------
    list
    """
    d_lst: list[Any] = []
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
    """Allocate an empty (size 0) C double vector.

    Returns
    -------
    ctypes.POINTER(DVECTOR)
        A new vector; release it with :func:`del_dvector`.
    """
    d_vect = ctypes.POINTER(DVECTOR)()
    lsci.initDVector(ctypes.pointer(d_vect))
    return d_vect


lsci.NewDVector.argtypes = [ctypes.POINTER(ctypes.POINTER(DVECTOR)), ctypes.c_size_t]
lsci.NewDVector.restype = None


def new_dvector(v_lst):
    """Create a C double vector from a list.

    Parameters
    ----------
    v_lst : list[float]
        The values. Non-numeric values (e.g. ``"NA"``) are stored as the missing
        value (see :func:`libscientific.misc.missing_value`).

    Returns
    -------
    ctypes.POINTER(DVECTOR)
        A new vector; release it with :func:`del_dvector`.

    Examples
    --------
    >>> v = new_dvector([1.0, 2.0])
    >>> dvector_tolist(v)
    [1.0, 2.0]
    >>> del_dvector(v)
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
    """Free a C double vector created by this module.

    Parameters
    ----------
    dvect : ctypes.POINTER(DVECTOR)
        The vector to release. It must not be used afterwards.
    """
    lsci.DelDVector(ctypes.pointer(dvect))


lsci.DVectorResize.argtypes = [ctypes.POINTER(DVECTOR), ctypes.c_size_t]
lsci.DVectorResize.restype = None


def dvector_resize(dvect, size):
    """Resize a C double vector; all its values are set to zero.

    Parameters
    ----------
    dvect : ctypes.POINTER(DVECTOR)
        The vector.
    size : int
        The new size.
    """
    lsci.DVectorResize(dvect, size)


lsci.PrintDVector.argtypes = [ctypes.POINTER(DVECTOR)]
lsci.PrintDVector.restype = None


def print_dvector(dvect):
    """Print a C double vector to standard output.

    Parameters
    ----------
    dvect : ctypes.POINTER(DVECTOR)
        The vector to print.
    """
    lsci.PrintDVector(dvect)


lsci.DVectorHasValue.argtypes = [ctypes.POINTER(DVECTOR), ctypes.c_double]
lsci.DVectorHasValue.restype = ctypes.c_int


def dvector_has_value(dvect, val):
    """Check whether a C double vector contains a value.

    Parameters
    ----------
    dvect : ctypes.POINTER(DVECTOR)
        The vector to search.
    val : float
        The value to look for (compared within the library tolerance ``EPSILON``).

    Returns
    -------
    int
        **0 if the value is present, 1 if it is not** (C library convention).

    Examples
    --------
    >>> v = new_dvector([1.0, 2.0])
    >>> dvector_has_value(v, 2.0), dvector_has_value(v, 5.0)
    (0, 1)
    """
    return lsci.DVectorHasValue(dvect, val)


lsci.DVectorSet.argtypes = [ctypes.POINTER(DVECTOR), ctypes.c_double]
lsci.DVectorSet.restype = None


def dvector_set(dvect, val):
    """Set every element of a C double vector to the same value.

    Parameters
    ----------
    dvect : ctypes.POINTER(DVECTOR)
        The vector, modified in place.
    val : float
        The value to assign.
    """
    lsci.DVectorSet(dvect, val)


lsci.setDVectorValue.argtypes = [ctypes.POINTER(DVECTOR), ctypes.c_size_t, ctypes.c_double]
lsci.setDVectorValue.restype = None


def set_dvector_value(dvect, indx, val):
    """Set one element of a C double vector.

    Parameters
    ----------
    dvect : ctypes.POINTER(DVECTOR)
        The vector, modified in place.
    indx : int
        Index of the element (0-based, smaller than the size).
    val : float
        The new value.
    """
    lsci.setDVectorValue(dvect, indx, val)


lsci.getDVectorValue.argtypes = [ctypes.POINTER(DVECTOR), ctypes.c_size_t]
lsci.getDVectorValue.restype = ctypes.c_double


def get_dvector_value(dvect, indx):
    """Return one element of a C double vector.

    Parameters
    ----------
    dvect : ctypes.POINTER(DVECTOR)
        The vector.
    indx : int
        Index of the element (0-based, smaller than the size).

    Returns
    -------
    float
    """
    return lsci.getDVectorValue(dvect, indx)


lsci.DVectorAppend.argtypes = [ctypes.POINTER(DVECTOR), ctypes.c_double]
lsci.DVectorAppend.restype = None


def dvector_append(dvect, val):
    """Append a value at the end of a C double vector.

    Parameters
    ----------
    dvect : ctypes.POINTER(DVECTOR)
        The vector; its size grows by one.
    val : float
        The value to append.
    """
    return lsci.DVectorAppend(dvect, val)


def dvector_tolist(dvect: Any) -> list[float]:
    """Convert a C double vector into a Python list.

    Parameters
    ----------
    dvect : ctypes.POINTER(DVECTOR)
        The vector.

    Returns
    -------
    list[float]
    """
    return xvector_tolist(dvect)


lsci.DVectorRemoveAt.argtypes = [ctypes.POINTER(DVECTOR), ctypes.c_size_t]
lsci.DVectorRemoveAt.restype = None


def dvector_remove_at(dvect, indx):
    """Remove the element at an index; the following elements shift left.

    Parameters
    ----------
    dvect : ctypes.POINTER(DVECTOR)
        The vector; its size decreases by one.
    indx : int
        Index of the element to remove.
    """
    return lsci.DVectorRemoveAt(dvect, indx)


lsci.DVectorCopy.argtypes = [ctypes.POINTER(DVECTOR), ctypes.POINTER(DVECTOR)]
lsci.DVectorCopy.restype = None


def dvector_copy(dvect_src):
    """Return an independent copy of a C double vector.

    Parameters
    ----------
    dvect_src : ctypes.POINTER(DVECTOR)
        The vector to copy.

    Returns
    -------
    ctypes.POINTER(DVECTOR)
        A new vector; release it with :func:`del_dvector`.
    """
    dvect_dst = init_dvector()
    lsci.DVectorCopy(dvect_src, dvect_dst)
    return dvect_dst


lsci.initUIVector.argtypes = [ctypes.POINTER(ctypes.POINTER(UIVECTOR))]
lsci.initUIVector.restype = None


def init_uivector():
    """Allocate an empty (size 0) C unsigned int vector.

    Returns
    -------
    ctypes.POINTER(UIVECTOR)
        A new vector; release it with :func:`del_uivector`.
    """
    uivect = ctypes.POINTER(UIVECTOR)()
    lsci.initUIVector(ctypes.pointer(uivect))
    return uivect


lsci.NewUIVector.argtypes = [ctypes.POINTER(ctypes.POINTER(UIVECTOR)), ctypes.c_size_t]
lsci.NewUIVector.restype = None


def new_uivector(v_lst):
    """Create a C unsigned int vector from a list.

    Parameters
    ----------
    v_lst : list[int]
        The values. They must be non-negative.

    Returns
    -------
    ctypes.POINTER(UIVECTOR)
        A new vector; release it with :func:`del_uivector`.

    Examples
    --------
    >>> v = new_uivector([1, 2])
    >>> uivector_tolist(v)
    [1, 2]
    >>> del_uivector(v)
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
    """Free a C unsigned int vector created by this module.

    Parameters
    ----------
    uivect : ctypes.POINTER(UIVECTOR)
        The vector to release. It must not be used afterwards.
    """
    lsci.DelUIVector(ctypes.pointer(uivect))


lsci.UIVectorResize.argtypes = [ctypes.POINTER(UIVECTOR), ctypes.c_size_t]
lsci.UIVectorResize.restype = None


def uivector_resize(uivect, size):
    """Resize a C unsigned int vector; all its values are set to zero.

    Parameters
    ----------
    uivect : ctypes.POINTER(UIVECTOR)
        The vector.
    size : int
        The new size.
    """
    lsci.UIVectorResize(uivect, size)


lsci.PrintUIVector.argtypes = [ctypes.POINTER(UIVECTOR)]
lsci.PrintUIVector.restype = None


def print_uivector(uivect):
    """Print a C unsigned int vector to standard output.

    Parameters
    ----------
    uivect : ctypes.POINTER(UIVECTOR)
        The vector to print.
    """
    lsci.PrintUIVector(uivect)


lsci.setUIVectorValue.argtypes = [ctypes.POINTER(UIVECTOR), ctypes.c_size_t, ctypes.c_size_t]
lsci.setUIVectorValue.restype = None


def set_uivector_value(uivect, indx, val):
    """Set one element of a C unsigned int vector.

    Parameters
    ----------
    uivect : ctypes.POINTER(UIVECTOR)
        The vector, modified in place.
    indx : int
        Index of the element (0-based, smaller than the size).
    val : int
        The new value.
    """
    lsci.setUIVectorValue(uivect, indx, val)


lsci.getUIVectorValue.argtypes = [ctypes.POINTER(UIVECTOR), ctypes.c_size_t]
lsci.getUIVectorValue.restype = ctypes.c_size_t


def get_uivector_value(uivect, indx):
    """Return one element of a C unsigned int vector.

    Parameters
    ----------
    uivect : ctypes.POINTER(UIVECTOR)
        The vector.
    indx : int
        Index of the element (0-based, smaller than the size).

    Returns
    -------
    int
    """
    return lsci.getUIVectorValue(uivect, indx)


lsci.UIVectorAppend.argtypes = [ctypes.POINTER(UIVECTOR), ctypes.c_size_t]
lsci.UIVectorAppend.restype = None


def uivector_append(uivect, val):
    """Append a value at the end of a C unsigned int vector.

    Parameters
    ----------
    uivect : ctypes.POINTER(UIVECTOR)
        The vector; its size grows by one.
    val : int
        The value to append.
    """
    return lsci.UIVectorAppend(uivect, val)


def uivector_tolist(uivect: Any) -> list[int]:
    """Convert a C unsigned int vector into a Python list.

    Parameters
    ----------
    uivect : ctypes.POINTER(UIVECTOR)
        The vector.

    Returns
    -------
    list[int]
    """
    return xvector_tolist(uivect)


lsci.UIVectorRemoveAt.argtypes = [ctypes.POINTER(UIVECTOR), ctypes.c_size_t]
lsci.UIVectorRemoveAt.restype = None


def uivector_remove_at(uivect, indx):
    """Remove the element at an index; the following elements shift left.

    Parameters
    ----------
    uivect : ctypes.POINTER(UIVECTOR)
        The vector; its size decreases by one.
    indx : int
        Index of the element to remove.
    """
    return lsci.UIVectorRemoveAt(uivect, indx)


def uivector_copy(uivect_src):
    """Return an independent copy of a C unsigned int vector.

    Parameters
    ----------
    uivect_src : ctypes.POINTER(UIVECTOR)
        The vector to copy.

    Returns
    -------
    ctypes.POINTER(UIVECTOR)
        A new vector; release it with :func:`del_uivector`.
    """
    return new_uivector(xvector_tolist(uivect_src))


lsci.initIVector.argtypes = [ctypes.POINTER(ctypes.POINTER(IVECTOR))]
lsci.initIVector.restype = None


def init_ivector():
    """Allocate an empty (size 0) C int vector.

    Returns
    -------
    ctypes.POINTER(IVECTOR)
        A new vector; release it with :func:`del_ivector`.
    """
    ivect = ctypes.POINTER(IVECTOR)()
    lsci.initIVector(ctypes.pointer(ivect))
    return ivect


lsci.NewIVector.argtypes = [ctypes.POINTER(ctypes.POINTER(IVECTOR)), ctypes.c_size_t]
lsci.NewIVector.restype = None


def new_ivector(v_lst):
    """Create a C int vector from a list.

    Parameters
    ----------
    v_lst : list[int]
        The values.

    Returns
    -------
    ctypes.POINTER(IVECTOR)
        A new vector; release it with :func:`del_ivector`.

    Examples
    --------
    >>> v = new_ivector([-1, 2])
    >>> ivector_tolist(v)
    [-1, 2]
    >>> del_ivector(v)
    """
    size = len(v_lst)
    ivect = ctypes.POINTER(IVECTOR)()
    lsci.NewIVector(ctypes.pointer(ivect), size)
    for i in range(size):
        lsci.setIVectorValue(ivect, i, int(v_lst[i]))
    return ivect


lsci.DelIVector.argtypes = [ctypes.POINTER(ctypes.POINTER(IVECTOR))]
lsci.DelIVector.restype = None


def del_ivector(ivect):
    """Free a C int vector created by this module.

    Parameters
    ----------
    ivect : ctypes.POINTER(IVECTOR)
        The vector to release. It must not be used afterwards.
    """
    lsci.DelIVector(ctypes.pointer(ivect))


# lsci.PrintIVector.argtypes = [ctypes.POINTER(IVECTOR)]
# lsci.PrintIVector.restype = None


# def print_ivector(dvect):
#    """
#    PrintIVector: Print to video a libscientific int vector
#    """
#    lsci.PrintIVector(dvect)


lsci.setIVectorValue.argtypes = [ctypes.POINTER(IVECTOR), ctypes.c_size_t, ctypes.c_int]
lsci.setIVectorValue.restype = None


def set_ivector_value(ivect, indx, val):
    """Set one element of a C int vector.

    Parameters
    ----------
    ivect : ctypes.POINTER(IVECTOR)
        The vector, modified in place.
    indx : int
        Index of the element (0-based, smaller than the size).
    val : int
        The new value.
    """
    lsci.setIVectorValue(ivect, indx, val)


lsci.getIVectorValue.argtypes = [ctypes.POINTER(IVECTOR), ctypes.c_size_t]
lsci.getIVectorValue.restype = ctypes.c_int


def get_ivector_value(ivect, indx):
    """Return one element of a C int vector.

    Parameters
    ----------
    ivect : ctypes.POINTER(IVECTOR)
        The vector.
    indx : int
        Index of the element (0-based, smaller than the size).

    Returns
    -------
    int
    """
    return lsci.getIVectorValue(ivect, indx)


def ivector_tolist(ivect: Any) -> list[int]:
    """Convert a C int vector into a Python list.

    Parameters
    ----------
    ivect : ctypes.POINTER(IVECTOR)
        The vector.

    Returns
    -------
    list[int]
    """
    return xvector_tolist(ivect)


lsci.IVectorAppend.argtypes = [ctypes.POINTER(IVECTOR), ctypes.c_int]
lsci.IVectorAppend.restype = None


def ivector_append(ivect, val):
    """Append a value at the end of a C int vector.

    Parameters
    ----------
    ivect : ctypes.POINTER(IVECTOR)
        The vector; its size grows by one.
    val : int
        The value to append.
    """
    return lsci.IVectorAppend(ivect, val)


lsci.IVectorRemoveAt.argtypes = [ctypes.POINTER(IVECTOR), ctypes.c_size_t]
lsci.IVectorRemoveAt.restype = None


def ivector_remove_at(ivect, indx):
    """Remove the element at an index; the following elements shift left.

    Parameters
    ----------
    ivect : ctypes.POINTER(IVECTOR)
        The vector; its size decreases by one.
    indx : int
        Index of the element to remove.
    """
    return lsci.IVectorRemoveAt(ivect, indx)


def _check_index(key: int, size: int) -> int:
    """Return a valid 0-based index (negative indexes count from the end)."""
    if not isinstance(key, int):
        raise TypeError(f"vector indices must be integers, not {type(key).__name__}")
    idx = key + size if key < 0 else key
    if not 0 <= idx < size:
        raise IndexError(f"index {key} out of range for a vector of size {size}")
    return idx


class DVector:
    """A double vector stored in the C library, with a Pythonic interface.

    Parameters
    ----------
    dvect_ : list[float], DVector or None
        Initial values (copied), or ``None`` for an empty vector.

    Attributes
    ----------
    dvect : ctypes.POINTER(DVECTOR)
        The underlying C vector, owned by this object.

    Notes
    -----
    ``v[i]`` reads/writes one element (negative indexes count from the end and
    out-of-range indexes raise ``IndexError``), ``len(v)`` is the size and
    iterating yields the values.

    Examples
    --------
    >>> v = DVector([1.0, 2.0, 3.0])
    >>> v.append(4.0)
    >>> len(v), v[-1]
    (4, 4.0)
    >>> v.tolist()
    [1.0, 2.0, 3.0, 4.0]
    """

    def __init__(self, dvect_=None):
        """Create the vector (see the class documentation)."""
        if dvect_ is None:
            self.dvect = init_dvector()
        elif isinstance(dvect_, DVector):
            self.dvect = dvector_copy(dvect_.dvect)
        elif isinstance(dvect_, ctypes.POINTER(DVECTOR)):
            self.dvect = dvector_copy(dvect_)
        else:
            self.dvect = new_dvector(dvect_)

    def __del__(self):
        """Release the C vector."""
        if getattr(self, "dvect", None) is not None:
            del_dvector(self.dvect)
            del self.dvect
            self.dvect = None

    def __repr__(self):
        return f"DVector({self.tolist()})"

    def __str__(self):
        return f"DVector({self.tolist()})"

    def __len__(self):
        """Return the number of elements."""
        return self.size()

    def __iter__(self):
        """Iterate over the values."""
        return iter(self.tolist())

    def __getitem__(self, key):
        """Return the element at an index (``v[i]``).

        Raises
        ------
        IndexError
            If the index is out of range.
        """
        return self.data_ptr()[_check_index(key, self.size())]

    def __setitem__(self, key, value):
        """Set the element at an index (``v[i] = value``).

        Raises
        ------
        IndexError
            If the index is out of range.
        """
        set_dvector_value(self.dvect, _check_index(key, self.size()), value)

    def size(self):
        """Return the number of elements.

        Returns
        -------
        int
        """
        return self.dvect[0].size

    def data_ptr(self):
        """Return the C data pointer (low-level access).

        Returns
        -------
        ctypes pointer
            ``data_ptr()[i]`` is element ``i``; valid while the vector exists.
        """
        return self.dvect[0].data

    def append(self, value):
        """Append a value at the end.

        Parameters
        ----------
        value : float
            The value to append.
        """
        return dvector_append(self.dvect, value)

    def extend(self, v_lst):
        """Append all the values of a list.

        Parameters
        ----------
        v_lst : list[float]
            The values to append.
        """
        for val in v_lst:
            dvector_append(self.dvect, val)

    def tolist(self):
        """Return the values as a Python list.

        Returns
        -------
        list[float]
        """
        return xvector_tolist(self.dvect)

    def fromlist(self, vlst_):
        """Replace the content with the values of a list.

        Parameters
        ----------
        vlst_ : list[float]
            The new values.
        """
        del_dvector(self.dvect)
        del self.dvect
        self.dvect = new_dvector(vlst_)

    def debug(self):
        """Print the vector to standard output."""
        print_dvector(self.dvect)


class UIVector:
    """A unsigned int vector stored in the C library, with a Pythonic interface.

    Parameters
    ----------
    uivect_ : list[int], UIVector or None
        Initial values (copied), or ``None`` for an empty vector.

    Attributes
    ----------
    uivect : ctypes.POINTER(UIVECTOR)
        The underlying C vector, owned by this object.

    Notes
    -----
    ``v[i]`` reads/writes one element (negative indexes count from the end and
    out-of-range indexes raise ``IndexError``), ``len(v)`` is the size and
    iterating yields the values.

    Examples
    --------
    >>> v = UIVector([1, 2, 3])
    >>> v.append(4)
    >>> len(v), v[-1]
    (4, 4)
    >>> v.tolist()
    [1, 2, 3, 4]
    """

    def __init__(self, uivect_=None):
        """Create the vector (see the class documentation)."""
        if uivect_ is None:
            self.uivect = init_uivector()
        elif isinstance(uivect_, UIVector):
            self.uivect = uivector_copy(uivect_.uivect)
        elif isinstance(uivect_, ctypes.POINTER(UIVECTOR)):
            self.uivect = uivector_copy(uivect_)
        else:
            self.uivect = new_uivector(uivect_)

    def __del__(self):
        """Release the C vector."""
        if getattr(self, "uivect", None) is not None:
            del_uivector(self.uivect)
            del self.uivect
            self.uivect = None

    def __repr__(self):
        return f"UIVector({self.tolist()})"

    def __str__(self):
        return f"UIVector({self.tolist()})"

    def __len__(self):
        """Return the number of elements."""
        return self.size()

    def __iter__(self):
        """Iterate over the values."""
        return iter(self.tolist())

    def __getitem__(self, key):
        """Return the element at an index (``v[i]``).

        Raises
        ------
        IndexError
            If the index is out of range.
        """
        return self.data_ptr()[_check_index(key, self.size())]

    def __setitem__(self, key, value):
        """Set the element at an index (``v[i] = value``).

        Raises
        ------
        IndexError
            If the index is out of range.
        """
        set_uivector_value(self.uivect, _check_index(key, self.size()), value)

    def size(self):
        """Return the number of elements.

        Returns
        -------
        int
        """
        return self.uivect[0].size

    def data_ptr(self):
        """Return the C data pointer (low-level access).

        Returns
        -------
        ctypes pointer
            ``data_ptr()[i]`` is element ``i``; valid while the vector exists.
        """
        return self.uivect[0].data

    def append(self, value):
        """Append a value at the end.

        Parameters
        ----------
        value : int
            The value to append.
        """
        return uivector_append(self.uivect, value)

    def extend(self, v_lst):
        """Append all the values of a list.

        Parameters
        ----------
        v_lst : list[int]
            The values to append.
        """
        for val in v_lst:
            uivector_append(self.uivect, val)

    def tolist(self):
        """Return the values as a Python list.

        Returns
        -------
        list[int]
        """
        return xvector_tolist(self.uivect)

    def fromlist(self, vlst_):
        """Replace the content with the values of a list.

        Parameters
        ----------
        vlst_ : list[int]
            The new values.
        """
        del_uivector(self.uivect)
        del self.uivect
        self.uivect = new_uivector(vlst_)

    def debug(self):
        """Print the vector to standard output."""
        print_uivector(self.uivect)
