"""
info libscientific python binding

Copyright (C) <2019>  Giuseppe Marco Randazzo

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
from libscientific.loadlibrary import LoadLibrary

lsci = LoadLibrary()

lsci.GetScientificVersion.argtypes = (None)
lsci.GetScientificVersion.restype = ctypes.c_char_p


class Scientific():
    """
    Libscientific info
    """
    def __init__(self):
        self.__version__ = str(lsci.GetScientificVersion().decode('UTF-8'))

    def print_version(self):
        """
        Print to video libscientific version
        """
        print(self.__version__)

    def version(self):
        """
        Get libscientific version as string
        """
        return self.__version__

if __name__ == '__main__':
    libscientific = Scientific()
    print(libscientific.__version__)
