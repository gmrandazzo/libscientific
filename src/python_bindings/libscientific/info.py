"""info.py libscientific python binding

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

lsci = load_libscientific_library()

lsci.GetScientificVersion.argtypes = None
lsci.GetScientificVersion.restype = ctypes.c_char_p


class Scientific():
    """
    Libscientific Information

    This class provides methods to retrieve and print information
    about the libscientific library.

    Attributes:
        __version__ (str): The version of the libscientific library.

    Methods:
        print_version(): Prints the libscientific version to the console.
        version(): Returns the libscientific version as a string.

    Examples:
    >>> sci = Scientific()
    >>> sci.print_version()
    1.0.0
    >>> ver = sci.version()
    >>> print(ver)
    1.0.0
    """
    def __init__(self):
        """
        Initialize the Scientific class and retrieve the libscientific version.
        """
        self.__version__ = str(lsci.GetScientificVersion().decode('UTF-8'))

    def print_version(self):
        """
        Print the libscientific version to the console.
        """
        print(self.__version__)

    def version(self):
        """
        Get the libscientific version as a string.

        Returns:
            str: The libscientific version.
        """
        return self.__version__


if __name__ == '__main__':
    libscientific = Scientific()
    print(libscientific.__version__)
