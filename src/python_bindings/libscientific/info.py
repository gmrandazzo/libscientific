# info.py libscientific python binding
#
# Copyright (C) <2023>  Giuseppe Marco Randazzo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""Information about the libscientific C library (version)."""

import ctypes

from libscientific.loadlibrary import load_libscientific_library

lsci = load_libscientific_library()

lsci.GetScientificVersion.argtypes = []
lsci.GetScientificVersion.restype = ctypes.c_char_p


class Scientific:
    """Information about the loaded libscientific C library.

    Attributes
    ----------
    __version__ : str
        Version of the C library, e.g. ``"1.8.0"``.

    Examples
    --------
    >>> from libscientific.info import Scientific
    >>> isinstance(Scientific().version(), str)
    True
    >>> Scientific().version()  # doctest: +SKIP
    '1.8.0'
    """

    def __init__(self):
        """Read the version from the C library."""
        self.__version__ = str(lsci.GetScientificVersion().decode("UTF-8"))

    def print_version(self):
        """Print the version of the C library to standard output."""
        print(self.__version__)

    def version(self):
        """Return the version of the C library.

        Returns
        -------
        str
            For example ``"1.8.0"``. It normally matches ``libscientific.__version__``
            (the version of the Python package).
        """
        return self.__version__


if __name__ == "__main__":
    libscientific = Scientific()
    print(libscientific.__version__)
