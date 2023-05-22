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

  if(colaverage->size  == 0 && colscaling->size == 0){
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
  else{
    if(trans->row == 0){
      ResizeMatrix(trans, orig->row, orig->col);
    }

    if(colaverage->size > 0){
      for(j = 0; j < trans->col; j++){
        for(i = 0; i < trans->row; i++){
          trans->data[i][j] = orig->data[i][j] - colaverage->data[j];
        }
      }
    }
    else{
      for(j = 0; j < orig->col; j++){
        for(i = 0; i < orig->row; i++){
          trans->data[i][j] = orig->data[i][j];
        }
      }
    }

    if(colscaling->size > 0){
      for(j = 0; j < trans->col; j++){
        if(FLOAT_EQ(colscaling->data[j], 0.f, 1e-2)){
          for(i = 0; i < trans->row; i++){
            trans->data[i][j] = 0.f;
          }
        }
        else{
          for(i = 0; i < trans->row; i++){
            trans->data[i][j] /= colscaling->data[j];
          }
        }
      }
    }
  }
}

void MatrixWhitening(matrix *X,
                     matrix *whitening_matrix,
                     matrix *X_whiten)
{
  size_t i;
  matrix *Xcov;
  matrix *D;
  matrix *u;
  matrix *u_T;
  matrix *vt;
  matrix *s;
  matrix *a;
  matrix *X_whiten_;
  matrix *Xt;

  if(whitening_matrix->row == 0){
    NewMatrix(&Xt, X->col, X->row);
    MatrixTranspose(X, Xt);
    initMatrix(&Xcov);
    MatrixCovariance(X, Xcov);
    DelMatrix(&Xt);

    initMatrix(&u);
    initMatrix(&s);
    initMatrix(&vt);
    SVDlapack(Xcov, u, s, vt);

    NewMatrix(&D, s->col, s->col);
    for(i = 0; i < s->col; i++)
      D->data[i][i] = 1/sqrt(s->data[i][i]);
    
    NewMatrix(&u_T, u->col, u->row);
    MatrixTranspose(u, u_T);
    NewMatrix(&a, D->row, u_T->col);
    MatrixDotProduct(D, u_T, a);
    DelMatrix(&u_T);
    
    ResizeMatrix(whitening_matrix, u->row, a->col);
    MatrixDotProduct(u, a, whitening_matrix);

    DelMatrix(&a);
    DelMatrix(&u);
    DelMatrix(&s);
    DelMatrix(&vt);

    NewMatrix(&X_whiten_, X->row, X->col);
    MatrixDotProduct(whitening_matrix, X, X_whiten_);
    MatrixCopy(X_whiten_, &X_whiten);
    
    DelMatrix(&X_whiten_);
    DelMatrix(&D);
    DelMatrix(&Xcov);
  }
  else{
    NewMatrix(&X_whiten_, X->row, X->col);
    MatrixDotProduct(X, whitening_matrix, X_whiten_);
    MatrixCopy(X_whiten_, &X_whiten);
    DelMatrix(&X_whiten_);
  }
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

