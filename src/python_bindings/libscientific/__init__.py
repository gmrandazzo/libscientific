"""
init.py

Library Initialization for LibScientific

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

Submodules
----------
- info : Module providing information about the LibScientific library.
- misc : Module for miscellaneous library utilities and missing value definitions.
- vector : Module for handling vector-related operations.
- matrix : Module for handling matrix-related operations.
- tensor : Module for handling tensor-related operations.
- pca : Module for Principal Component Analysis (PCA).
- pls : Module for Partial Least Squares (PLS).
- cpca : Module for Consensus Principal Component Analysis (CPCA).
- interpolate : Module for cubic spline interpolation.
- vectlist : Module for handling lists of double vectors.
- clustering : Module for clustering algorithms.
- model_io : Module for high-performance zero-dependency model serialization.

Attributes
----------
__version__ : str
    The version of the LibScientific library.

Examples
--------
The following example demonstrates how to import and use the LibScientific library:

>>> import libscientific
>>> print(libscientific.__version__)
'1.5.1'
>>> x = [[1.0, 2.0], [3.0, 4.0]]
>>> vec = libscientific.vector.new_dvector(x)
>>> print(libscientific.vector.dvector_tolist(vec))
[[1.0, 2.0], [3.0, 4.0]]

Note
----
This module serves as the main entry point for using the LibScientific library.
Users can import submodules and access various functionalities provided by the library.
"""

from . import (
    clustering,
    cpca,
    info,
    interpolate,
    matrix,
    misc,
    model_io,
    pca,
    pls,
    tensor,
    vectlist,
    vector,
)
from .cpca import CPCA
from .interpolate import CubicSplineInterpolation
from .matrix import Matrix
from .pca import PCA
from .pls import PLS
from .tensor import Tensor
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
    "CubicSplineInterpolation",
]


# __version__ = info.Scientific().__version__
__version__ = "1.8.0"
