/* upca.h
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

#ifndef UPCA_H
#define UPCA_H
#include "list.h"
#include "tensor.h"
#include "matrix.h"
#include "vector.h"
#include "preprocessing.h"
#include "scientificinfo.h"

typedef struct{
  matrix *scores;
  tensor *loadings;
  dvector *varexp;
  dvectorlist *colaverage;
  dvectorlist *colscaling;
} UPCAMODEL;

#define UPCACONVERGENCE 1e-10

void NewUPCAModel(UPCAMODEL **m);
void DelUPCAModel(UPCAMODEL **m);

int CheckTensor(tensor *X_);
void UPCA(tensor *X_, size_t npc, size_t autoscaling, UPCAMODEL *m, ssignal *s);

void UPCAScorePredictor(tensor *X_,
                        UPCAMODEL *model,
                        size_t npc,
                        matrix *pscores);

void UPCAIndVarPredictor(matrix *T,
                         tensor *P,
                         dvectorlist *colaverage,
                         dvectorlist *colscaling,
                         size_t npc,
                         tensor *X);

void PrintUPCAModel(UPCAMODEL *m);

#endif
