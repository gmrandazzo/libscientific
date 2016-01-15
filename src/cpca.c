#include "cpca.h"

#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "array.h"
#include "dvector.h"
#include "uivector.h"


NewCPCAModel(CPCAMODEL **m){
  initArray(&(*m)->b_scores);
  initArray(&(*m)->b_loadings);
  initMatrix(&(*m)->sscores);
  initMatrix(&(*m)->sweights);
  initDVector(&(*m)->b_scaling);
  initDVector(&(*m)->expvar);
  initDVector(&(*m)->colaverage);
  initDVector(&(*m)->colscaling);
}

 
DelCPCAModel(CPCAMODEL **m)
{
  DelArray(&(*m)->b_scores);
  DelArray(&(*m)->b_loadings);
  DelMatrix(&(*m)->sscores);
  DelMatrix(&(*m)->sweights);
  DelDVector(&(*m)->b_scaling);
  DelDVector(&(*m)->expvar);
  DelDVector(&(*m)->colaverage);
  DelDVector(&(*m)->colscaling);
}
 
 
void CPCA(array *x, size_t npc, size_t scaling, CPCAMODEL *m)
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
  MatrixColAverage(mx, &(model->colaverage));
  for(j = 0; j < mx->col; j++){
    for(i = 0; i < mx->row; i++){
      setMatrixValue(E, i, j, getMatrixValue(mx, i, j) - getDVectorValue(model->colaverage, j));
    }
  }
  
  if(scaling > 0){
    if(scaling == 1){
      MatrixColSDEV(mx, &(model->colscaling));
    }
    else if(scaling == 2){
      MatrixColRMS(mx, &(model->colscaling));
    }
    else if(scaling == 3){ /* PARETO Autoscaling */
      MatrixColSDEV(mx, &(model->colscaling));
      for(i = 0; i < model->colscaling->size; i++){
        setDVectorValue(model->colscaling, i, sqrt(getDVectorValue(model->colscaling, i)));
      }
    }
    else if(scaling == 4){ /* Range Scaling */
      for(i = 0; i < mx->col; i++){
        MatrixColumnMinMax(mx, i, &min, &max);
        DVectorAppend(&model->colscaling, (max - min));
      }
    }
    else if(scaling == 5){ /* Level Scaling  */
      DVectorCopy(model->colaverage, &model->colscaling);
    }
    else{
      for(int i = 0; i < model->colaverage->size; i++){
        DVectorAppend(&model->colscaling, 1.0);
      }
    }
    
    for(j = 0; j < E->col; j++){
      if(FLOAT_EQ(getDVectorValue(model->colscaling, j), 0, EPSILON)){
        for(i = 0; i< E->row; i++){
          E->data[i][j] = 0.f;
        }
      }
      else{
        for(i = 0; i < E->row; i++){
          E->data[i][j] /= model->colscaling->data[j];
        }
      }
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
  
}