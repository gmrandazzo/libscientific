"""io.py libscientific python binding

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
from libscientific.loadlibrary import load_libscientific_library
import libscientific.matrix as mx
import libscientific.vector as vect
from .pca import PCAMODEL
from .cpca import CPCAMODEL
from .pls import PLSMODEL

lsci = load_libscientific_library()

lsci.WritePCA.argtypes = [
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(PCAMODEL)
]
lsci.WritePCA.restype = None


def write_pca(dbpath : str, pca : PCAMODEL):
    """
    write_pca(dbpath, pca):
    
    Writes a PCA (Principal Component Analysis) model to a SQLite database.

    Parameters:
        dbpath (str): The path to the SQLite database file where the PCA 
            model will be stored.
        pca (PCAMODEL): The libscientific data structure representing the PCA 
            model.

    Returns:
        None
    """
    lsci.WritePCA(dbpath, pca)


lsci.ReadPCA.argtypes = [
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(PCAMODEL)
]
lsci.ReadPCA.restype = None

def read_pca(dbpath : str, pca : PCAMODEL):
    """
    read_pca(dbpath, pca):
    
    Reads a PCA (Principal Component Analysis) model from a SQLite database.

    Parameters:
        dbpath (str): The path to the SQLite database file from which the PCA
            model will be read.
        pca (PCAMODEL): The libscientific data structure where the PCA
            model will be loaded.

    Returns:
        None
    """
    lsci.ReadPCA(dbpath, pca)


lsci.WriteCPCA.argtypes = [
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(CPCAMODEL)
]
lsci.WriteCPCA.restype = None

def write_cpca(dbpath : str, pca : CPCAMODEL):
    """
    write_pca(dbpath, cpca):
    
    Writes a CPCA (Consensus Principal Component Analysis) model to
    a SQLite database.

    Parameters:
        dbpath (str): The path to the SQLite database file where the CPCA 
            model will be stored.
        cpca (CPCAMODEL): The libscientific data structure representing the CPCA 
            model.

    Returns:
        None
    """
    lsci.WriteCPCA(dbpath, pca)


lsci.ReadCPCA.argtypes = [
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(CPCAMODEL)
]
lsci.ReadCPCA.restype = None

def read_cpca(dbpath : str, pca : CPCAMODEL):
    """
    read_cpca(dbpath, cpca):
    
    Reads a CPCA (Consensus Principal Component Analysis) model from
    a SQLite database.

    Parameters:
        dbpath (str): The path to the SQLite database file from which the CPCA
            model will be read.
        cpca (CPCAMODEL): The libscientific data structure where the CPCA
            model will be loaded.

    Returns:
        None
    """
    lsci.ReadCPCA(dbpath, pca)


lsci.WritePLS.argtypes = [
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(PLSMODEL)
]
lsci.WritePLS.restype = None


def write_pls(dbpath : str, pca : PLSMODEL):
    """
    write_pls(dbpath, pls):
    
    Writes a PLS (Partial Least Squares) model to
    a SQLite database.

    Parameters:
        dbpath (str): The path to the SQLite database file where the PLS 
            model will be stored.
        pls (PLSMODEL): The libscientific data structure representing the PLS 
            model.

    Returns:
        None
    """
    lsci.WritePLS(dbpath, pca)


lsci.ReadPLS.argtypes = [
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(PLSMODEL)
]
lsci.ReadPLS.restype = None

def read_pls(dbpath : str, pca : PLSMODEL):
    """
    read_pls(dbpath, pls):
    
    Reads a PLS (Partial Least Squares) model from a SQLite database.

    Parameters:
        dbpath (str): The path to the SQLite database file from which the PLS
            model will be read.
        pls (PLSMODEL): The libscientific data structure where the PLS
            model will be loaded.

    Returns:
        None
    """
    lsci.ReadPLS(dbpath, pca)


