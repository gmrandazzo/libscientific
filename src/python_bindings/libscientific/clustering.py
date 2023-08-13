"""clustering.py libscientific python binding

Copyright (C) <2023>  Giuseppe Marco Randazzo

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
import os
import libscientific.matrix as mx
import libscientific.vector as vect
from libscientific.loadlibrary import load_libscientific_library


lsci = load_libscientific_library()

lsci.MDC.argtypes = [ctypes.POINTER(mx.MATRIX),
                     ctypes.c_size_t,
                     ctypes.c_int,
                     ctypes.POINTER(vect.UIVECTOR),
                     ctypes.c_size_t]
lsci.MDC.restype = None

def most_descriptive_compound(x_input, nobjects):
    """
    Most Descriptive Compound algorithm
    """
    if "Matrix" not in str(type(x_input)):
        x_input_ = mx.new_matrix(x_input)
        xalloc = True
    else:
        x_input_ = x_input

    obj_sel = vect.init_uivector()
    lsci.MDC(x_input_,
             nobjects,
             0,
             obj_sel,
             os.cpu_count())
    obj_sel_lst = vect.uivector_tolist(obj_sel)
    vect.del_uivector(obj_sel)
    if xalloc is True:
        mx.del_matrix(x_input_)
        del x_input_
    return obj_sel_lst


lsci.MaxDis_Fast.argtypes = [ctypes.POINTER(mx.MATRIX),
                             ctypes.c_size_t,
                             ctypes.c_int,
                             ctypes.POINTER(vect.UIVECTOR),
                             ctypes.c_size_t]
lsci.MaxDis_Fast.restype = None


def max_dissimilarity_selection(x_input, nobjects):
    """
    Max dissimilarity compound selection algorithm
    """
    if "Matrix" not in str(type(x_input)):
        x_input_ = mx.new_matrix(x_input)
        xalloc = True
    else:
        x_input_ = x_input

    obj_sel = vect.init_uivector()
    lsci.MaxDis_Fast(x_input_,
                     nobjects,
                     0,
                     obj_sel,
                     os.cpu_count())
    obj_sel_lst = vect.uivector_tolist(obj_sel)
    vect.del_uivector(obj_sel)
    if xalloc is True:
        mx.del_matrix(x_input_)
        del x_input_
    return obj_sel_lst


lsci.KMeans.argtypes = [ctypes.POINTER(mx.MATRIX),
                        ctypes.c_size_t,
                        ctypes.c_int,
                        ctypes.POINTER(vect.UIVECTOR),
                        ctypes.POINTER(mx.MATRIX),
                        ctypes.c_size_t]
lsci.KMeans.restype = None


def k_means_plus_plus(x_input, n_clusters):
    """
    K-Means++ clustering using david arthur initialization
    """
    if "Matrix" not in str(type(x_input)):
        x_input_ = mx.new_matrix(x_input)
        xalloc = True
    else:
        x_input_ = x_input

    labels = vect.init_uivector()
    lsci.KMeans(x_input_,
                n_clusters,
                1,
                labels,
                None,
                os.cpu_count())
    labels_lst = vect.uivector_tolist(labels)
    vect.del_uivector(labels)
    if xalloc is True:
        mx.del_matrix(x_input_)
        del x_input_
    return labels_lst
