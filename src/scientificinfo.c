/* scientificinfo.c
*
* Copyright (C) <2016>  Giuseppe Marco Randazzo
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "scientificinfo.h"

#define major_ 0
#define minor_ 7
#define patch_ 5

void ScientificVersion()
{
  printf("Scientific Library was writen by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>\nVersion: %d.%d.%d\n", major_, minor_, patch_);
}

const char *GetScientificVersion()
{
  static char c[10];
  sprintf(c, "%d.%d.%d", (int)major_, (int)minor_, (int)patch_);
  return c;
}
