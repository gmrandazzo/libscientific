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
    lsci_lib_dir = os.environ.get("LIBSCIENTIFIC_LIB_DIR")
    lsci = None
    if os.name == "nt":
        if lsci_lib_dir:
            lib_path = f'{pathlib.Path(lsci_lib_dir)}/libscientific.dll'
        else:
            lib_path = f'{pathlib.Path(__file__).parent}/libscientific.dll'
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
                curr_path = str(pathlib.Path(__file__).parent.absolute())
                os_path = os.environ.get('PATH', '')
                os.environ['PATH'] = f"{os_path};{curr_path}"
                lib_path = f'{curr_path}/{libname}'
        try:
            lsci = ctypes.cdll.LoadLibrary(lib_path)
        except OSError as err:
            msg = "Please install sqlite3 and lapack library "
            msg += "if you want to use the binary library in this package.\n"
            msg += "Alternativelly you can install the library directly "
            msg += "from its source code from here "
            msg += "https://github.com/gmrandazzo/libscientific"
            msg += "\nIf the library is installed please "
            msg += "check/specify the location in LD_LIBRARY_PATH "
            msg += "or add this environment variable "
            msg += "export LIBSCIENTIFIC_LIB_DIR=<path libscientific>"
            raise RuntimeError(msg) from err
    else:
        raise RuntimeError(f"Don't know how to load library on {os.name}")
    return lsci
