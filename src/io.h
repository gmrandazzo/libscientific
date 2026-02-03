/* io.h
*
* Copyright (C) <2023>  Giuseppe Marco Randazzo
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
#include "pca.h"
#include "pls.h"
#include "cpca.h"
#include "upca.h"
#include "upls.h"
#include "scientificinfo.h"

/**
 * Save a PCA model to sqlite session
 * 
 * @param [in] dbname output path of th sqlite3 database
 * @param [out] pca PCA model to save
 */
void WritePCA(char *dbpath, PCAMODEL *pca);

/**
 * Read a PLS model from an sqlite session
 * 
 * @param [in] dbname input path of th sqlite3 database
 * @param [out] pca PCA model to fill
 */
void ReadPCA(char *dbpath, PCAMODEL *pca);


/**
 * Save a CPCA model to sqlite session
 * 
 * @param [in] dbname output path of th sqlite3 database
 * @param [out] cpca CPCA model to save
 */
void WriteCPCA(char *dbpath, CPCAMODEL *cpca);

/**
 * Read a CPCA model from an sqlite session
 * 
 * @param [in] dbname input path of th sqlite3 database
 * @param [out] cpca CPCA model to fill
 */
void ReadCPCA(char *dbpath, CPCAMODEL *cpca);

void WritePLS(char *dbpath, PLSMODEL *pls);
void ReadPLS(char *dbpath, PLSMODEL *pls);

/**
 * Save a matrix to a CSV file
 * 
 * @param [in] path output path of the CSV file
 * @param [in] m matrix to save
 */
void WriteMatrixCSV(char *path, matrix *m);

#endif
