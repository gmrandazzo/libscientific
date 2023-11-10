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

def load_library_for_nt():
    """Load the libscientific library for NT systems
    """
    try:
        lib_path = f'{pathlib.Path(__file__).parent}'
        return ctypes.WinDLL(f'{lib_path}\\libscientific.dll')
    except TypeError:
        return None

def get_posix_library():
    """Find the path to the libscientific POSIX library.

    Returns
    -------
    str or None
        The path to the libscientific library file if found, or None if not found.
    """
    paths = ["/usr/lib/libscientific.so",
             "/usr/lib64/libscientific.so",
             "/usr/lib/x86_64-linux-gnu/libscientific.so",
             "/usr/local/lib/libscientific.so",
             "/usr/local/lib64/libscientific.so",
             "/opt/homebrew/lib/libscientific.dylib"]
    for path in paths:
        if pathlib.Path(path).is_file() is True:
            return path
    return None

def load_library_for_posix():
    """Load the libscientific library for posix systems
    """
    lib_path = get_posix_library()
    if lib_path is not None:
        return ctypes.cdll.LoadLibrary(lib_path)
    try:
        lib_path = str(pathlib.Path(__file__).parent.absolute())
        if os.uname()[0] == "Darwin":
            return ctypes.CDLL(f'{lib_path}/libscientific.dylib')
        # else load library for linux systems
        for i in range(3, 5):
            if pathlib.Path(f'{lib_path}/libgfortran.so.{i}').is_file():
                _ = ctypes.CDLL(f'{lib_path}/libgfortran.so.{i}')
                break
        _ = ctypes.CDLL(f'{lib_path}/libquadmath.so.0.0.0')
        _ = ctypes.CDLL(f'{lib_path}/libblas.so')
        _ = ctypes.CDLL(f'{lib_path}/liblapack.so')
        return ctypes.CDLL(f'{lib_path}/libscientific.so')
    except OSError as err:
        msg = f"Don't know how to load libscientific on {os.name}\n"
        msg += "Please install the library directly "
        msg += "from its source code "
        msg += "https://github.com/gmrandazzo/libscientific"
        msg += "\nIf the library is installed please "
        msg += "check/specify the location in LD_LIBRARY_PATH "
        msg += "with export LD_LIBRARY_PATH=<path libscientific>"
        raise RuntimeError(msg) from err

def load_libscientific_library():
    """Load the libscientific library dynamically.

    Returns
    -------
    CDLL
        A ctypes dynamic link library object representing the loaded libscientific library.

    Raises
    ------
    RuntimeError
        If the library cannot be loaded or the platform is not supported.
    """
    if os.name == "nt":
        return load_library_for_nt()
    return load_library_for_posix()
