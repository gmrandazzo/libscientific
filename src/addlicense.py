#!/usr/bin/env python


import string
import sys

fi = open(sys.argv[1], "r")
fo = open(sys.argv[2], "w")

fo.write("/* %s\n" % (sys.argv[1]))
fo.write("*\n* Copyright (C) <2016>  Giuseppe Marco Randazzo\n*\n* This program is free software: you can redistribute it and/or modify\n* it under the terms of the GNU General Public License as published by\n* the Free Software Foundation, either version 3 of the License, or\n* (at your option) any later version.\n*\n* This program is distributed in the hope that it will be useful,\n* but WITHOUT ANY WARRANTY; without even the implied warranty of\n* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n* GNU General Public License for more details.\n*\n* You should have received a copy of the GNU General Public License\n* along with this program.  If not, see <http://www.gnu.org/licenses/>.\n*/\n\n")

for line in fi:
    fo.write(line)
fi.close()
fo.close()
