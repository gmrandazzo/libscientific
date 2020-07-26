/* io.h
 * 
 * Copyright (C) <2020>  Giuseppe Marco Randazzo
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

#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stdlib.h>
#include <matrix.h>
#include <vector.h>
#include "pls.h"

/*
 * Serialization/Deserialization methods
 */
void WriteTensor(tensor *, char **);
void ReadTensor(char *, tensor **);

void WriteMatrix(matrix *, char **);
void ReadMatrix(char *, matrix **);

void WriteDVector(dvector *, char **);
void ReadDVector(char *, dvector **);

/*
 * Save a PLSMODEL to file
 */
void SavePLSModel(PLSMODEL *, char *path);

/*
 * Load a PLSMODEL from file
 */
void LoadPLSModel(char *path,PLSMODEL *);

