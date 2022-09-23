/* cpca.c
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

#include "cpca.h"

#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "numeric.h"
#include "array.h"
#include "vector.h"


void NewCPCAModel(CPCAMODEL **m){
  initArray(&(*m)->b_scores);
  initArray(&(*m)->b_loadings);
  initMatrix(&(*m)->sscores);
  initMatrix(&(*m)->sweights);
  initDVector(&(*m)->b_scaling);
  initDVector(&(*m)->expvar);
  initDVectorList(&(*m)->colaverage);
  initDVectorList(&(*m)->colscaling);
}


void DelCPCAModel(CPCAMODEL **m)
{
  DelArray(&(*m)->b_scores);
  DelArray(&(*m)->b_loadings);
  DelMatrix(&(*m)->sscores);
  DelMatrix(&(*m)->sweights);
  DelDVector(&(*m)->expvar);
  DelDVectorList(&(*m)->colaverage);
  DelDVectorList(&(*m)->colscaling);
}


void CPCA(array *x, size_t npc, size_t scaling, CPCAMODEL *model)
{
  size_t i, j, pc;
  dvector *t;
  dvector *p;
  dvector *colvar;
  dvector *eval; /* t't */

  double min, max, mod_p, mod_t_old, mod_t_new, ss;

  array *E; /* data matrix of autoscaled / mean centred object */
  NewArray(&E, x->order);
  for(i = 0; i < E->order; i++){
    NewArrayMatrix(&E, i, x->m[i]->row, x->m[i]->row);

  }

  /*CENTERING */
  for(k = 0; k < x->order; k++){
    dvector *colavg;
    initDVector(&colavg);
    MatrixColAverage(x->m[k], &colavg);
    for(j = 0; j < x->m[i]->col; j++){
      for(i = 0; i < x->m[i]->row; i++){
        E->m[k]->m->data[i][j] = x->m[k]->data[i][j] - colavg->data[j];
      }
    }
    /* Append column average to dvectorlist */
    DVectorListAppend(&model->colaverage, colavg);
    DelDVector(&colavg);
  }
  if(scaling > 0){
    for(k = 0; k < x->order; k++){
      dvector *colscaling;
      initDVector(&colscaling);
      if(scaling == 1){
        MatrixColSDEV(x->m[k], &colscaling);
      }
      else if(scaling == 2){
        MatrixColSDEV(x->m[k], &colscaling);
      }
      else if(scaling == 3){ /* PARETO Autoscaling */
        MatrixColSDEV(x->m[k], &colscaling);
        for(i = 0; i < colscaling->size; i++){
          colscaling->data[i] = sqrt(colscaling->data[i])
        }
      }
      else if(scaling == 4){ /* Range Scaling */
        for(i = 0; i < mx->col; i++){
          MatrixColumnMinMax(x->m[k], i, &min, &max);
          DVectorAppend(colscaling, (max - min));
        }
      }
      else if(scaling == 5){ /* Level Scaling  */
        DVectorCopy(model->colaverage[k]->d, &model->colscaling);
      }
      else{
        for(int i = 0; i < model->colaverage[k]->d->size; i++){
          DVectorAppend(colscaling, 1.0);
        }
      }

      for(j = 0; j < E->col; j++){
        if(FLOAT_EQ(colscaling->data[j], 0.f, EPSILON)){
          for(i = 0; i< E->row; i++){
            E->m[k]->data[i][j] = 0.f;
          }
        }
        else{
          for(i = 0; i < E->row; i++){
            E->m[k]->data[i][j] /= model->colscaling[k]->data[j];
          }
        }
      }
      DVectorListAppend(model->colscaling, colscaling);
    }
  }

   /* if the number of principal component selected is major of the permitted */
  if(npc > E->col)
    npc = E->col;    /* set the value to the max value */

  ss = 0.f;
  for(i = 0; i < E->row; i++){
    for(j = 0; j < E->col; j++)
      ss += square(E->data[i][j]);
  }

  NewDVector(&t, E->row);
  NewDVector(&p, E->col);
  NewDVector(&eval, npc);


  ResizeMatrix(&(model->scores), E->row, npc);
  ResizeMatrix(&(model->loadings), E->col, npc);

  for(i = 0; i < npc; i++){
    for(j = 0; j < )
  }

  DelDVector(&t);
  DelDVector(&p);
  DelDVector(&eval);
}
