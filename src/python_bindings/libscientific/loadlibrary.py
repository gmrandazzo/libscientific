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

from ctypes import cdll
from ctypes.util import find_library

def load_libscientific_library():
    """
    Load the libscientific library
    """
    lsci = None
    library = find_library("scientific")
    if library is None:
        message = "Please install libscientific.\n"
        message += "Go to https://github.com/gmrandazzo/libscientific\n"
        message += "or if it is installed please check/specify the location in LD_LIBRARY_PATH"
        print(message)
    else:
        lsci = cdll.LoadLibrary(library)
    return lsci


if __name__ == '__main__':
    a = load_libscientific_library()
