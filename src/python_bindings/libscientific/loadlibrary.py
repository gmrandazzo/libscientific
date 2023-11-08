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
    lsci = None
    if os.name == "nt":
        lib_path = f'{pathlib.Path(__file__).parent}'
        try:
            _ = ctypes.CDLL(f'{lib_path}/libsqlite3-0.dll', winmode=0)
            _ = ctypes.CDLL(f'{lib_path}/libwinpthread-1.dll', winmode=0)
            lsci = ctypes.CDLL(f'{lib_path}/libscientific.dll', winmode=0)
        except TypeError:
            lsci = ctypes.CDLL(f'{lib_path}/libscientific.dll', winmode=0)

    elif os.name == "posix":
        lib_path = get_posix_library()
        if lib_path is not None:
            lsci = ctypes.cdll.LoadLibrary(lib_path)
        else:
            try:
                lib_path = str(pathlib.Path(__file__).parent.absolute())
                if os.uname()[0] == "Darwin":
                    lsci = ctypes.CDLL(f'{lib_path}/libscientific.dylib')
                else:
                    _ = ctypes.CDLL(f'{lib_path}/libgfortran.so.3')
                    _ = ctypes.CDLL(f'{lib_path}/libquadmath.so.0')
                    _ = ctypes.CDLL(f'{lib_path}/libblas.so.3')
                    _ = ctypes.CDLL(f'{lib_path}/liblapack.so.3')
                    lsci = ctypes.CDLL(f'{lib_path}/libscientific.so')
            except OSError as err:
                msg = "Please install sqlite3 and lapack library "
                msg += "if you want to use the binary library in this package.\n"
                msg += "Alternativelly you can install the library directly "
                msg += "from its source code from here "
                msg += "https://github.com/gmrandazzo/libscientific"
                msg += "\nIf the library is installed please "
                msg += "check/specify the location in LD_LIBRARY_PATH "
                msg += "with export LD_LIBRARY_PATH=<path libscientific>"
                raise RuntimeError(msg) from err
    else:
        raise RuntimeError(f"Don't know how to load library on {os.name}")
    return lsci
