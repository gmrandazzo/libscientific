"""loadlibrary.py method to load libscientific library

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

import os
import pathlib
import ctypes
from ctypes.util import find_library


def get_posix_library():
    """Get libscientific posix library location
    """
    paths = ["/usr/lib/libscientific.so",
             "/usr/lib64/libscientific.so",
             "/usr/local/lib/libscientific.so",
             "/usr/local/lib64/libscientific.so"] 
    for path in paths:
        if pathlib.Path(path).is_file() is True:
            return path
    return None


def load_libscientific_library():
    """Load the libscientific library
    """
    lsci_lib_dir = os.environ.get("LIBSCIENTIFIC_LIB_DIR")
    lsci = None
    if os.name == "nt":
        if lsci_lib_dir:
            lib_path = str(pathlib.Path(lsci_lib_dir))+'libscientific.dll'
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

        if lsci_lib_dir:
            lib_path = str(pathlib.Path(lsci_lib_dir))+f'/{libname}'
        else:
            lib_path = get_posix_library()
            if lib_path is None:
                lib_path = find_library("scientific")
            else:
                msg = "Please install libscientific. Go to "
                msg += "https://github.com/gmrandazzo/libscientific"
                msg += "or if it is installed please "
                msg += "check/specify the location in LD_LIBRARY_PATH "
                msg += "or add this environment variable "
                msg += "export LIBSCIENTIFIC_LIB_DIR=<path libscientific>"
                raise RuntimeError(msg)
        lsci = ctypes.cdll.LoadLibrary(lib_path)
    else:
        raise RuntimeError(f"Don't know how to load library on {os.name}")
    return lsci


if __name__ == '__main__':
    a = load_libscientific_library()
