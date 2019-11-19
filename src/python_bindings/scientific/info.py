# info libscientific python binding
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
from scientific.loadlibrary import LoadLibrary

lsci = LoadLibrary()

lsci.GetScientificVersion.argtypes = (None)
lsci.initMatrix.restype = ctypes.c_char_p


class Scientific(object):
    def __init__(self):
        # ctypes.cast(lsci.GetScientificVersion(), ctypes.c_char_p)
        self.__version__ = "1.2.0"


if __name__ == '__main__':
    libscientific = Scientific()
    print(libscientific.__version__)
