# loadlibrary libscientific python binding
#
# Copyright (C) <2019>  Giuseppe Marco Randazzo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import pathlib
import ctypes
from ctypes.util import find_library


def get_posix_library():
    paths = ["/usr/lib/libscientific.so",
             "/usr/lib64/libscientific.so",
             "/usr/local/lib/libscientific.so",
             "/usr/local/lib64/libscientific.so"]
    
    for p in paths:
        if pathlib.Path(p).is_file() == True:
            return p
        else:
            continue
    return None


def load_libscientific_library():
    """Load the libscientific library
    """
    LIBSCIENTIFIC_LIB_DIR = os.environ.get("LIBSCIENTIFIC_LIB_DIR")
    lsci = None
    if os.name == "nt":
        if LIBSCIENTIFIC_LIB_DIR:
            lib_path = str(pathlib.Path(LIBSCIENTIFIC_LIB_DIR))+'libscientific.dll'
        else:
            lib_path = find_library("scientific")

        try:
            lsci = ctypes.CDLL(lib_path, winmode=0)
        except TypeError:
            lsci = ctypes.CDLL(lib_path)

    elif os.name == "posix":
        libname = None
        if os.uname()[0] == "Darwin":
            libname = "libscientific.dylib"
        else:
            libname = "libscientific.so"

        if LIBSCIENTIFIC_LIB_DIR:
            lib_path = str(pathlib.Path(LIBSCIENTIFIC_LIB_DIR))+f'/{libname}'
        else:
            lib_path = get_posix_library()
            if lib_path is None:
                lib_path = find_library("scientific")
            else:
                message = "Please install libscientific. "
                message += "Go to https://github.com/gmrandazzo/libscientific"
                message += "or if it is installed please check/specify the location in LD_LIBRARY_PATH"
                message += "or add this environment variable export LIBSCIENTIFIC_LIB_DIR=<path libscientific>"
                raise RuntimeError(message) 
        lsci = ctypes.cdll.LoadLibrary(lib_path)
    else:
        raise RuntimeError(f"Don't know how to load library on {os.name}")
    return lsci


if __name__ == '__main__':
    a = load_libscientific_library()
