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
#include "modelvalidation.h"
#include <math.h>

void NewEPLSModel(EPLSMODEL** m)
{
  (*m) = xmalloc(sizeof(EPLSMODEL));
  (*m)->models = NULL;
  (*m)->q2 = NULL;
  (*m)->sdep = NULL;
  (*m)->n_models = 0;
  (*m)->nlv = 0;
  (*m)->ny = 0;
}

void DelEPLSModel(EPLSMODEL** m)
{
  size_t i;
  for(i = 0; i < (*m)->n_models; i++){
    DelPLSModel(&(*m)->models[i]);
    DelMatrix(&(*m)->q2[i]);
    DelMatrix(&(*m)->sdep[i]);
  }
  xfree((*m)->sdep);
  xfree((*m)->q2);
  xfree((*m)->models);
  xfree((*m));
}

void EPLS(matrix *mx, matrix *my, size_t nlv, size_t xautoscaling, size_t yautoscaling, size_t n_models, double testsize, EPLSMODEL *m, ELearningMethod agtype, ssignal *s)
{
  if(agtype == Bagging){
    matrix *x_train, *y_train, *x_test, *y_test, *p_y_test;
    uivector *testids;
    size_t i;

    m->models = xmalloc(sizeof(PLSMODEL)*n_models);
    m->q2 = xmalloc(sizeof(matrix)*n_models);
    m->sdep = xmalloc(sizeof(matrix)*n_models);
    m->n_models = n_models;
    m->nlv = nlv;
    m->ny = my->col;
    unsigned int srand_init = mx->row * testsize + xautoscaling + yautoscaling + mx->col + my->col;
    for(i = 0; i < n_models; i++){
      initMatrix(&x_train);
      initMatrix(&y_train);
      initMatrix(&x_test);
      initMatrix(&y_test);
      initUIVector(&testids);
      train_test_split(mx, my, testsize, &x_train, &y_train, &x_test, &y_test, &testids, &srand_init);
      srand_init++;
      NewPLSModel(&m->models[i]);
      PLS(x_train, y_train, nlv, xautoscaling, yautoscaling, m->models[i], s);
      initMatrix(&p_y_test);
      PLSYPredictorAllLV(x_test, m->models[i], nlv, &p_y_test);
      initMatrix(&m->q2[i]);
      initMatrix(&m->sdep[i]);
      PLSRegressionStatistics(y_test, p_y_test, &m->q2[i], &m->sdep[i], NULL);
      DelMatrix(&p_y_test);
      DelMatrix(&x_train);
      DelMatrix(&y_train);
      DelMatrix(&x_test);
      DelMatrix(&y_test);
    }
  }
  else{
    //Random Subspace Method
  }
}

void EPLSGetSXScore(EPLSMODEL *m, CombinationRule crule, matrix *sxscores){}
void EPLSGetSXLoadings(EPLSMODEL *m, CombinationRule crule, matrix *sxloadings){}
void EPLSGetSYScore(EPLSMODEL *m, CombinationRule crule, matrix *syscores){}
void EPLSGetSYLoadings(EPLSMODEL *m, CombinationRule crule, matrix *syloadings){}
void EPLSGetSWeights(EPLSMODEL *m, CombinationRule crule, matrix *sweights){}
void EPLSGetSBetaCoefficients(EPLSMODEL *m, CombinationRule crule, matrix *sbetas){}
void EPLSYPrediction(matrix *mx, EPLSMODEL *m, CombinationRule crule, matrix *py){}

void PrintEPLSModel(EPLSMODEL *m)
{
  size_t i, j, k;
  for(j = 0; j < m->n_models; j++){
    PrintMatrix(m->q2[j]);
    i = k = 1;
  }


  for(k = 0; k < m->ny; k++){
    printf("Y No. %lu\n", k);
    puts("Q2 external");
    for(i = 0; i < m->nlv; i++){
      printf("%lu;", i+1);
      for(j = 0; j < m->n_models-1; j++){
        printf("%f;", m->q2[j]->data[i][k]);
      }
      printf("%f\n", m->q2[m->n_models-1]->data[i][k]);
    }
    puts("SDEP external");
    for(i = 0; i < m->nlv; i++){
      printf("%lu;", i+1);
      for(j = 0; j < m->n_models-1; j++){
        printf("%f;", m->sdep[j]->data[i][k]);
      }
      printf("%f\n", m->sdep[m->n_models-1]->data[i][k]);
    }
    printf(">>> END\n");
  }
}
