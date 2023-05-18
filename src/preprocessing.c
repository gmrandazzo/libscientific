#include "preprocessing.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "tensor.h"
#include <math.h>

void MatrixPreprocess(matrix *orig,
                      int type,
                      dvector *colaverage,
                      dvector *colscaling,
                      matrix *trans)
{
  size_t i, j;
  double min, max;
  /* check if matrix have nan or infinite and convert them to MISSING */
  MatrixCheck(orig);

  if(type >= 0){
    /* CENTERING */
    MatrixColAverage(orig, colaverage);
    for(j = 0; j < orig->col; j++){
      for(i = 0; i < orig->row; i++){
        if(FLOAT_EQ(orig->data[i][j], MISSING, 1e-1)){
          continue;
        }
        else{
          trans->data[i][j] = orig->data[i][j] - colaverage->data[j];
        }
      }
    }

    if(type == 1){
      MatrixColSDEV(orig, colscaling);
    }
    else if(type == 2){
      MatrixColRMS(orig, colscaling);
    }
    else if(type == 3){ /* PARETO Autoscaling */
      MatrixColSDEV(orig, colscaling);
      for(i = 0; i < colscaling->size; i++){
        colscaling->data[i] = sqrt(colscaling->data[i]);
      }
    }
    else if(type == 4){ /* Range Scaling */
      for(i = 0; i < orig->col; i++){
        MatrixColumnMinMax(orig, i, &min, &max);
        DVectorAppend(colscaling, (max - min));
      }
    }
    else if(type == 5){ /* Level Scaling  */
      DVectorCopy(colaverage, colscaling);
    }
    else{
      for(i = 0; i < colaverage->size; i++){
        DVectorAppend(colscaling, 1.0);
      }
    }

    for(j = 0; j < trans->col; j++){
      if(FLOAT_EQ(getDVectorValue(colscaling, j), 0, EPSILON)){
        for(i = 0; i< trans->row; i++){
          trans->data[i][j] = 0.f;
        }
      }
      else{
        for(i = 0; i < trans->row; i++){
          if(FLOAT_EQ(trans->data[i][j], MISSING, 1e-1)){
            continue;
          }
          else{
            trans->data[i][j] /= colscaling->data[j];
          }
        }
      }
    }
  }
  else{
    /* Just copy AS IS*/
    for(j = 0; j < orig->col; j++){
      for(i = 0; i < orig->row; i++){
        trans->data[i][j] = orig->data[i][j];
      }
    }
  }
}

void MatrixWhitening(matrix *X, matrix *X_whiten)
{
  size_t i;
  matrix *Xcov, *D, *Dinv, *evect, *evectT, *eX, *DinveX;
  dvector *eval;
  NewMatrix(&Xcov, X->row, X->row);
  MatrixCovariance(X, Xcov);
  initDVector(&eval);
  initMatrix(&evect);
  EVectEval(Xcov, eval, evect);
  NewMatrix(&D, eval->size, eval->size);
  for(i = 0; i < eval->size; i++)
    D->data[i][i] = eval->data[i];

  NewMatrix(&Dinv, D->row, D->col);
  MatrixInversion(D, Dinv);

  NewMatrix(&evectT, evect->col, evect->row);
  MatrixTranspose(evect, evectT);

  /* here to predict we need evect/evectT and Dinv*/
  NewMatrix(&eX, evectT->row, X->col);
  MatrixDotProduct(evectT, X, eX);

  NewMatrix(&DinveX, Dinv->row, eX->col);
  MatrixDotProduct(Dinv, eX, DinveX);

  DelMatrix(&eX);
  ResizeMatrix(X_whiten, evect->row, DinveX->col);
  MatrixDotProduct(evect, DinveX, X_whiten);

  DelMatrix(&DinveX);
  DelMatrix(&evect);
  DelMatrix(&Dinv);
  DelMatrix(&evectT);
  DelMatrix(&D);
  DelDVector(&eval);
  DelMatrix(&Xcov);
}

void TensorPreprocess(tensor *orig,
                      int type,
                      dvectorlist *colaverages,
                      dvectorlist *colscalings,
                      tensor *trans)
{
  size_t k;
  dvector *colaverage;
  dvector *colscaling;

  for(k = 0; k < orig->order; k++){
    initDVector(&colaverage);
    initDVector(&colscaling);

    MatrixPreprocess(orig->m[k], type, colaverage, colscaling, trans->m[k]);

    DVectorListAppend(colaverages, colaverage);
    DVectorListAppend(colscalings, colscaling);

    DelDVector(&colaverage);
    DelDVector(&colscaling);
  }
}

