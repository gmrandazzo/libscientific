# clustering.py libscientific python binding
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
"""Object selection and clustering.

- :func:`most_descriptive_compound` and :func:`max_dissimilarity_selection`
  pick a representative subset of objects (rows), e.g. to build a diverse
  training set;
- :func:`k_means_pp_centers` picks well spread starting centers;
- :func:`k_means_plus_plus` assigns every object to one of ``k`` clusters.

All functions accept a list of rows, a numpy array or a
:class:`~libscientific.matrix.Matrix`, and work on multiple threads
(``nthreads`` defaults to the number of CPUs).

Distance metrics (``metric``): 0 = euclidean, 1 = manhattan, 2 = cosine
(``1 - cosine similarity``).

    >>> from libscientific.clustering import max_dissimilarity_selection
    >>> X = [[0.0, 0.0], [0.1, 0.0], [5.0, 5.0], [5.1, 5.0], [0.0, 5.0], [10.0, 0.0]]
    >>> sorted(max_dissimilarity_selection(X, 3))   # 3 far apart objects
    [3, 4, 5]
"""

from __future__ import annotations

import ctypes
import os

import libscientific.matrix as mx
import libscientific.vector as vect
from libscientific.loadlibrary import load_libscientific_library
from libscientific.matrix import MatrixLike

lsci = load_libscientific_library()

_SELECT_ARGTYPES = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.c_size_t,
    ctypes.c_int,
    ctypes.POINTER(vect.UIVECTOR),
    ctypes.c_size_t,
]
for _name in ("MDC", "MDC_Fast", "MaxDis_Fast"):
    getattr(lsci, _name).argtypes = _SELECT_ARGTYPES
    getattr(lsci, _name).restype = None


def _threads(nthreads: int | None) -> int:
    return nthreads if nthreads is not None else (os.cpu_count() or 1)


def _select(
    func_name: str, x_input: MatrixLike, nobjects: int, metric: int, nthreads: int | None
) -> list[int]:
    x_input_, owned = mx.as_c_matrix(x_input)
    obj_sel = vect.init_uivector()
    try:
        getattr(lsci, func_name)(x_input_, nobjects, metric, obj_sel, _threads(nthreads))
        return vect.uivector_tolist(obj_sel)
    finally:
        vect.del_uivector(obj_sel)
        if owned:
            mx.del_matrix(x_input_)


def most_descriptive_compound(
    x_input: MatrixLike, nobjects: int, metric: int = 0, nthreads: int | None = None
) -> list[int]:
    """
    Select the most descriptive objects (MDC algorithm).

    The first selected object is the one closest to all the others (the most
    representative); each following object is the one that best describes
    the objects not yet represented by the selection.

    Parameters
    ----------
    x_input : MatrixLike
        Data, shape ``(n_objects, n_features)``.
    nobjects : int
        Number of objects to select (at most ``n_objects``).
    metric : int, default 0
        Distance: 0 = euclidean, 1 = manhattan, 2 = cosine.
    nthreads : int, optional
        Number of threads (default: number of CPUs).

    Returns
    -------
    list[int]
        Row indexes of the selected objects, in selection order.

    See Also
    --------
    most_descriptive_compound_fast : same selection, lower memory use.

    Examples
    --------
    >>> X = [[0.0, 0.0], [0.1, 0.0], [0.2, 0.1], [5.0, 5.0], [5.1, 5.0]]
    >>> ids = most_descriptive_compound(X, 2)
    >>> len(ids), len(set(ids))
    (2, 2)
    """
    return _select("MDC", x_input, nobjects, metric, nthreads)


def most_descriptive_compound_fast(
    x_input: MatrixLike, nobjects: int, metric: int = 0, nthreads: int | None = None
) -> list[int]:
    """
    Memory-efficient version of :func:`most_descriptive_compound`.

    The distances are computed on the fly instead of being stored, so it
    works on data sets that do not fit an ``n_objects x n_objects`` distance
    matrix in memory. It selects the same objects as
    :func:`most_descriptive_compound`.

    Parameters
    ----------
    x_input : MatrixLike
        Data, shape ``(n_objects, n_features)``.
    nobjects : int
        Number of objects to select.
    metric : int, default 0
        Distance: 0 = euclidean, 1 = manhattan, 2 = cosine.
    nthreads : int, optional
        Number of threads (default: number of CPUs).

    Returns
    -------
    list[int]
        Row indexes of the selected objects, in selection order.
    """
    return _select("MDC_Fast", x_input, nobjects, metric, nthreads)


def max_dissimilarity_selection(
    x_input: MatrixLike, nobjects: int, metric: int = 0, nthreads: int | None = None
) -> list[int]:
    """
    Select objects as different as possible from each other (MaxMin).

    Each new object is the one whose distance to the closest already
    selected object is the largest, which spreads the selection over the
    whole data space (including outliers).

    Parameters
    ----------
    x_input : MatrixLike
        Data, shape ``(n_objects, n_features)``.
    nobjects : int
        Number of objects to select.
    metric : int, default 0
        Distance: 0 = euclidean, 1 = manhattan, 2 = cosine.
    nthreads : int, optional
        Number of threads (default: number of CPUs).

    Returns
    -------
    list[int]
        Row indexes of the selected objects, in selection order.

    Examples
    --------
    >>> X = [[0.0, 0.0], [0.1, 0.0], [5.0, 5.0], [5.1, 5.0], [0.0, 5.0], [10.0, 0.0]]
    >>> sorted(max_dissimilarity_selection(X, 3))
    [3, 4, 5]
    """
    return _select("MaxDis_Fast", x_input, nobjects, metric, nthreads)


lsci.KMeansppCenters.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.c_size_t,
    ctypes.POINTER(vect.UIVECTOR),
    ctypes.c_int,
]
lsci.KMeansppCenters.restype = None


def k_means_pp_centers(
    x_input: MatrixLike, n_centers: int, nthreads: int | None = None
) -> list[int]:
    """
    Choose well spread starting centers with the k-means++ seeding.

    The first center is chosen at random; each following one is drawn with a
    probability proportional to its squared distance from the nearest center
    already chosen (D. Arthur and S. Vassilvitskii, 2007).

    Parameters
    ----------
    x_input : MatrixLike
        Data, shape ``(n_objects, n_features)``.
    n_centers : int
        Number of centers.
    nthreads : int, optional
        Number of threads (default: number of CPUs).

    Returns
    -------
    list[int]
        Row indexes of the chosen objects. The choice is random: it can
        change between runs.

    Examples
    --------
    >>> X = [[0.0, 0.0], [0.1, 0.0], [5.0, 5.0], [5.1, 5.0], [0.0, 5.0], [10.0, 0.0]]
    >>> centers = k_means_pp_centers(X, 3)
    >>> len(centers), all(0 <= c < 6 for c in centers)
    (3, True)
    """
    x_input_, owned = mx.as_c_matrix(x_input)
    obj_sel = vect.init_uivector()
    try:
        lsci.KMeansppCenters(x_input_, n_centers, obj_sel, _threads(nthreads))
        return vect.uivector_tolist(obj_sel)
    finally:
        vect.del_uivector(obj_sel)
        if owned:
            mx.del_matrix(x_input_)


lsci.KMeans.argtypes = [
    ctypes.POINTER(mx.MATRIX),
    ctypes.c_size_t,
    ctypes.c_int,
    ctypes.POINTER(vect.UIVECTOR),
    ctypes.POINTER(mx.MATRIX),
    ctypes.c_size_t,
]
lsci.KMeans.restype = None


def k_means_plus_plus(
    x_input: MatrixLike, n_clusters: int, nthreads: int | None = None
) -> list[int]:
    """
    Cluster objects with k-means, initialized with k-means++.

    Parameters
    ----------
    x_input : MatrixLike
        Data, shape ``(n_objects, n_features)``.
    n_clusters : int
        Number of clusters ``k``.
    nthreads : int, optional
        Number of threads (default: number of CPUs).

    Returns
    -------
    list[int]
        Cluster label of every object, from 0 to ``k - 1``. The
        initialization is random: the numbering of the clusters (and, for
        poorly separated data, the clusters themselves) can change between
        runs.

    Examples
    --------
    Two well separated groups always end up in two different clusters:

    >>> X = [[0.0, 0.0], [0.1, 0.2], [0.2, 0.1], [9.0, 9.0], [9.1, 9.2], [9.2, 9.1]]
    >>> labels = k_means_plus_plus(X, 2)
    >>> len(labels), labels[0] == labels[1] == labels[2], labels[0] != labels[3]
    (6, True, True)
    """
    x_input_, owned = mx.as_c_matrix(x_input)
    labels = vect.init_uivector()
    try:
        lsci.KMeans(x_input_, n_clusters, 1, labels, None, _threads(nthreads))
        return vect.uivector_tolist(labels)
    finally:
        vect.del_uivector(labels)
        if owned:
            mx.del_matrix(x_input_)
