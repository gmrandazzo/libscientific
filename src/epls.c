/* epls.c
*
* Copyright (C) <2018>  Giuseppe Marco Randazzo
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

#include "epls.h"

#include "memwrapper.h"
#include "pls.h"
#include "pca.h" /*Using: MatrixAutoScaling(); and calcVarExpressed(); */
#include "numeric.h" /* Using:  if(FLOAT_EQ(NumOne, NumTwo));*/
#include "metricspace.h"
#include <math.h>

void NewEPLSModel(EPLSMODEL** m)
{
  (*m) = xmalloc(sizeof(EPLSMODEL));
  (*m)->models = NULL;
  (*m)->n_models = 0;
}

void DelEPLSModel(EPLSMODEL** m)
{
  size_t i;
  for(i = 0; i < (*m)->n_models; i++){
    DelPLSModel(&(*m)->models[i]);
  }
  xfree((*m)->models);
  xfree((*m));
}

void EPLS(matrix *mx, matrix *my, size_t nlv, size_t xautoscaling, size_t yautoscaling, EPLSMODEL *model, ELearningMethod agtype, ssignal *s){}
void EPLSGetSXScore(EPLSMODEL *m, CombinationRule crule, matrix *sxscores){}
void EPLSGetSXLoadings(EPLSMODEL *m, CombinationRule crule, matrix *sxloadings){}
void EPLSGetSYScore(EPLSMODEL *m, CombinationRule crule, matrix *syscores){}
void EPLSGetSYLoadings(EPLSMODEL *m, CombinationRule crule, matrix *syloadings){}
void EPLSGetSWeights(EPLSMODEL *m, CombinationRule crule, matrix *sweights){}
void EPLSGetSBetaCoefficients(EPLSMODEL *m, CombinationRule crule, matrix *sbetas){}
void EPLSYPrediction(matrix *mx, EPLSMODEL *m, CombinationRule crule, matrix *py){}
