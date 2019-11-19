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
#
import ctypes
from ctypes import util


def LoadLibrary():
    """
    Load the libscientific library
    """
    library = util.find_library("scientific")
    if library is None:
        print("Please install libscientific. see https://github.com/gmrandazzo/libscientific")
        return 0
    else:
        lsci = ctypes.cdll.LoadLibrary(library)
        return lsci


if __name__ == '__main__':
    a = LoadLibrary()
