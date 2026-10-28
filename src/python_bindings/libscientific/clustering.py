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

lsci.MDC.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.c_size_t,
    ctypes.c_int,
    ctypes.POINTER(vect.UIVECTOR),
    ctypes.c_size_t,
]
lsci.MDC.restype = None

lsci.MDC_Fast.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.c_size_t,
    ctypes.c_int,
    ctypes.POINTER(vect.UIVECTOR),
    ctypes.c_size_t,
]
lsci.MDC_Fast.restype = None


def most_descriptive_compound(x_input, nobjects, metric=0, nthreads=None):
    """
    Most Descriptive Compound Selection Algorithm.

    This function selects the most descriptive compounds/rows from
    an input matrix based on the specified number of objects to select.

    Parameters:
        x_input : List[List] or Matrix
            Input matrix.
        nobjects : int
            Number of objects to select.
        metric : int, optional
            Metric type: 0 euclidean (default), 1 manhattan, 2 cosine distance.
        nthreads : int, optional
            Number of CPU threads to use. Defaults to os.cpu_count().

    Returns:
        List[int] : A list of selected object/row IDs.

    Examples:
    >>> np.random.seed(12345)
    >>> x = np.random.rand(100, 2)
    >>> mdc_ids = most_descriptive_compound(x, 10)
    >>> mdc_ids
    [74, 97, 95, 7, 35, 25, 50, 10, 32, 59]
    """
    threads = nthreads if nthreads is not None else (os.cpu_count() or 1)
    mdcxalloc = False
    if isinstance(x_input, mx.Matrix):
        x_input_ = x_input.mtx
    elif isinstance(x_input, ctypes.POINTER(mx.MATRIX)):
        x_input_ = x_input
    elif "numpy" in str(type(x_input)):
        x_input_ = mx.matrix_from_numpy(x_input)
        mdcxalloc = True
    else:
        x_input_ = mx.new_matrix(x_input)
        mdcxalloc = True

    obj_sel = vect.init_uivector()
    lsci.MDC(x_input_, nobjects, metric, obj_sel, threads)
    obj_sel_lst = vect.uivector_tolist(obj_sel)
    vect.del_uivector(obj_sel)

    if mdcxalloc is True:
        mx.del_matrix(x_input_)
        del x_input_

    return obj_sel_lst


def most_descriptive_compound_fast(x_input, nobjects, metric=0, nthreads=None):
    """
    Most Descriptive Compound Selection Algorithm (Memory-Efficient Fast Implementation).
    """
    threads = nthreads if nthreads is not None else (os.cpu_count() or 1)
    mdcxalloc = False
    if isinstance(x_input, mx.Matrix):
        x_input_ = x_input.mtx
    elif isinstance(x_input, ctypes.POINTER(mx.MATRIX)):
        x_input_ = x_input
    elif "numpy" in str(type(x_input)):
        x_input_ = mx.matrix_from_numpy(x_input)
        mdcxalloc = True
    else:
        x_input_ = mx.new_matrix(x_input)
        mdcxalloc = True

    obj_sel = vect.init_uivector()
    lsci.MDC_Fast(x_input_, nobjects, metric, obj_sel, threads)
    obj_sel_lst = vect.uivector_tolist(obj_sel)
    vect.del_uivector(obj_sel)

    if mdcxalloc is True:
        mx.del_matrix(x_input_)
        del x_input_

    return obj_sel_lst


lsci.MaxDis_Fast.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.c_size_t,
    ctypes.c_int,
    ctypes.POINTER(vect.UIVECTOR),
    ctypes.c_size_t,
]
lsci.MaxDis_Fast.restype = None


def max_dissimilarity_selection(x_input, nobjects, metric=0, nthreads=None):
    """
    Max Dissimilarity Compound Selection Algorithm.

    This function selects compounds/rows from an input matrix in
    a way that maximizes their dissimilarity based on the specified
    number of objects to select.

    Parameters:
        x_input : List[List] or Matrix
            Input matrix.
        nobjects : int
            Number of objects to select.
        metric : int, optional
            Metric type: 0 euclidean (default), 1 manhattan, 2 cosine distance.
        nthreads : int, optional
            Number of CPU threads to use. Defaults to os.cpu_count().

    Returns:
        List[int] : A list of selected object/row IDs.

    Examples:
    >>> np.random.seed(12345)
    >>> x = np.random.rand(100, 2)
    >>> mdis_ids = max_dissimilarity_selection(x, 10)
    >>> mdis_ids
    [57, 89, 88, 6, 23, 94, 56, 61, 39, 24]
    """
    threads = nthreads if nthreads is not None else (os.cpu_count() or 1)
    mdisxalloc = False
    if isinstance(x_input, mx.Matrix):
        x_input_ = x_input.mtx
    elif isinstance(x_input, ctypes.POINTER(mx.MATRIX)):
        x_input_ = x_input
    elif "numpy" in str(type(x_input)):
        x_input_ = mx.matrix_from_numpy(x_input)
        mdisxalloc = True
    else:
        x_input_ = mx.new_matrix(x_input)
        mdisxalloc = True

    obj_sel = vect.init_uivector()
    lsci.MaxDis_Fast(x_input_, nobjects, metric, obj_sel, threads)
    obj_sel_lst = vect.uivector_tolist(obj_sel)
    vect.del_uivector(obj_sel)

    if mdisxalloc is True:
        mx.del_matrix(x_input_)
        del x_input_

    return obj_sel_lst


lsci.KMeansppCenters.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.c_size_t,
    ctypes.POINTER(vect.UIVECTOR),
    ctypes.c_int,
]
lsci.KMeansppCenters.restype = None


def k_means_pp_centers(x_input, n_centers, nthreads=None):
    """
    KMeans++ Centers Method: Select n_centers diverse seeds.
    """
    threads = nthreads if nthreads is not None else (os.cpu_count() or 1)
    alloc = False
    if isinstance(x_input, mx.Matrix):
        x_input_ = x_input.mtx
    elif isinstance(x_input, ctypes.POINTER(mx.MATRIX)):
        x_input_ = x_input
    elif "numpy" in str(type(x_input)):
        x_input_ = mx.matrix_from_numpy(x_input)
        alloc = True
    else:
        x_input_ = mx.new_matrix(x_input)
        alloc = True

    obj_sel = vect.init_uivector()
    lsci.KMeansppCenters(x_input_, n_centers, obj_sel, threads)
    obj_sel_lst = vect.uivector_tolist(obj_sel)
    vect.del_uivector(obj_sel)

    if alloc:
        mx.del_matrix(x_input_)
        del x_input_

    return obj_sel_lst


lsci.KMeans.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.c_size_t,
    ctypes.c_int,
    ctypes.POINTER(vect.UIVECTOR),
    ctypes.POINTER(mx.MATRIX),
    ctypes.c_size_t,
]
lsci.KMeans.restype = None


def k_means_plus_plus(x_input, n_clusters, nthreads=None):
    """
    K-Means++ Clustering Algorithm (kmeans + David Arthur initialization).

    This function performs K-Means++ clustering on the given input matrix
    to create the specified number of clusters.

    Parameters:
        x_input : List[List] or Matrix
            Input matrix.
        n_clusters : int
            Number of clusters.
        nthreads : int, optional
            Number of threads to use. Defaults to os.cpu_count().

    Returns:
        List[int] : A list containing the cluster assignments for each input point.

    Examples:
    >>> np.random.seed(12345)
    >>> x = np.random.rand(100, 2)
    >>> clusters = k_means_plus_plus(x, 3)
    >>> clusters
    [1, 2, 0, 0, 0, 0, 2, 2, 0, 0, 0, 1, ..., 0, 1, 1, 2, 2]
    """
    threads = nthreads if nthreads is not None else (os.cpu_count() or 1)
    kmppxalloc = False
    if isinstance(x_input, mx.Matrix):
        x_input_ = x_input.mtx
    elif isinstance(x_input, ctypes.POINTER(mx.MATRIX)):
        x_input_ = x_input
    elif "numpy" in str(type(x_input)):
        x_input_ = mx.matrix_from_numpy(x_input)
        kmppxalloc = True
    else:
        x_input_ = mx.new_matrix(x_input)
        kmppxalloc = True

    labels = vect.init_uivector()
    lsci.KMeans(x_input_, n_clusters, 1, labels, None, threads)
    labels_lst = vect.uivector_tolist(labels)
    vect.del_uivector(labels)
    if kmppxalloc is True:
        mx.del_matrix(x_input_)
        del x_input_
    return labels_lst
