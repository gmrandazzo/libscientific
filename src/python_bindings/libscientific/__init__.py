# libscientific python binding
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
"""libscientific: multivariate analysis and chemometrics in C, with Python bindings.

The package wraps the libscientific C library through ``ctypes``; the
compiled library is bundled in the wheel.

Models (all with ``fit``/``predict``, ``save``/``load`` and documented getters):

- :class:`PCA` - Principal Component Analysis (:mod:`libscientific.pca`)
- :class:`PLS` - Partial Least Squares regression (:mod:`libscientific.pls`)
- :class:`CPCA` - Consensus (multi-block) PCA (:mod:`libscientific.cpca`)
- :class:`MLR` - Multiple Linear Regression (:mod:`libscientific.mlr`)
- :class:`LDA` - Linear Discriminant Analysis (:mod:`libscientific.lda`)
- :class:`UPCA` - Multi-way PCA (:mod:`libscientific.upca`)
- :class:`UPLS` - Multi-way PLS (:mod:`libscientific.upls`)

Data structures and utilities: :mod:`~libscientific.matrix`,
:mod:`~libscientific.vector`, :mod:`~libscientific.tensor`,
:mod:`~libscientific.vectlist`, :mod:`~libscientific.clustering`,
:mod:`~libscientific.interpolate`, :mod:`~libscientific.model_io`,
:mod:`~libscientific.misc`, :mod:`~libscientific.info`.

Inputs can be Python lists of rows, numpy arrays or the library's own
:class:`Matrix`; results are returned as Python lists.

Examples
--------
>>> import libscientific
>>> model = libscientific.PCA(scaling=1, npc=2)
>>> model.fit([[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [3.0, 4.0, 2.0], [4.0, 3.0, 1.0]])
>>> len(model.get_scores()), len(model.get_scores()[0])
(4, 2)
"""

from . import (
    clustering,
    cpca,
    info,
    interpolate,
    lda,
    matrix,
    misc,
    mlr,
    model_io,
    pca,
    pls,
    tensor,
    upca,
    upls,
    vectlist,
    vector,
)
from .cpca import CPCA
from .interpolate import CubicSplineInterpolation
from .lda import LDA
from .matrix import Matrix
from .mlr import MLR
from .pca import PCA
from .pls import PLS
from .tensor import Tensor
from .upca import UPCA
from .upls import UPLS
from .vectlist import DVectorList
from .vector import DVector, UIVector

__all__ = [
    "info",
    "misc",
    "vector",
    "matrix",
    "tensor",
    "pca",
    "pls",
    "cpca",
    "mlr",
    "lda",
    "upca",
    "upls",
    "interpolate",
    "vectlist",
    "clustering",
    "model_io",
    "Matrix",
    "DVector",
    "UIVector",
    "Tensor",
    "DVectorList",
    "PCA",
    "PLS",
    "CPCA",
    "MLR",
    "LDA",
    "UPCA",
    "UPLS",
    "CubicSplineInterpolation",
]


__version__ = "1.8.0"
