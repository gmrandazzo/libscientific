/* Handles input/output operations and file formats.
 * Copyright (C) 2023-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef IO_H
#define IO_H
#include <stdio.h>
#include "pca.h"
#include "pls.h"
#include "cpca.h"
#include "mlr.h"
#include "lda.h"
#include "ica.h"
#include "upca.h"
#include "upls.h"

void WritePCA(char *dbpath, PCAMODEL *pca);
void ReadPCA(char *dbpath, PCAMODEL *pca);

void WritePLS(char *dbpath, PLSMODEL *pls);
void ReadPLS(char *dbpath, PLSMODEL *pls);

void WriteCPCA(char *dbpath, CPCAMODEL *cpca);
void ReadCPCA(char *dbpath, CPCAMODEL *cpca);

void WriteMLR(char *dbpath, MLRMODEL *mlr);
void ReadMLR(char *dbpath, MLRMODEL *mlr);

void WriteLDA(char *dbpath, LDAMODEL *lda);
void ReadLDA(char *dbpath, LDAMODEL *lda);

void WriteICA(char *dbpath, ICAMODEL *ica);
void ReadICA(char *dbpath, ICAMODEL *ica);

void WriteUPCA(char *dbpath, UPCAMODEL *upca);
void ReadUPCA(char *dbpath, UPCAMODEL *upca);

void WriteUPLS(char *dbpath, UPLSMODEL *upls);
void ReadUPLS(char *dbpath, UPLSMODEL *upls);

/**
 * Save a matrix to a CSV file
 * 
 * @param [in] path output path of the CSV file
 * @param [in] m matrix to save
 */
void WriteMatrixCSV(char *path, matrix *m);

#endif
