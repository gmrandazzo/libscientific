#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "memwrapper.h"
#include "upca.h"
#include "array.h"
#include "matrix.h"
#include "vector.h"
#include "pca.h"
#include "numeric.h"


void NewUPCAModel(UPCAMODEL** m)
{
  (*m) = xmalloc(sizeof(UPCAMODEL));
  initMatrix(&(*m)->scores);
  initArray(&(*m)->loadings);
  initDVector(&(*m)->varexp);
  initMatrix(&(*m)->colaverage);
  initMatrix(&(*m)->colscaling);
}

void DelUPCAModel(UPCAMODEL** m)
{
  DelMatrix(&(*m)->scores);
  DelArray(&(*m)->loadings);
  DelDVector(&(*m)->varexp);
  DelMatrix(&(*m)->colaverage);
  DelMatrix(&(*m)->colscaling);
  xfree((*m));
}

int CheckArray(array *X_){
  size_t i, j;
  for(i = 0; i < X_->order; i++){
    for(j = i+1; j < X_->order; j++){
      if(X_->m[i]->row == X_->m[j]->row
        && X_->m[i]->col == X_->m[j]->col)
        continue;
      else
        return -1;
    }
  }
  return 0;
}


/*
 * PCA Nipals algorithm for Multi-Way Principal Component Analysis 
 * Journal of Chemometrics, Vol 1, 41-56 (1987)
 * 
 * Svante Woold, Paul Geladi and Kim Esbensen
 * 
 * X_ = input Array
 * X = array MeanCentered and Autoscaled
 * P = Matrix 
 * 
 * step 1 : Preprocessing X_ data: X_ is mean centered and after autoscaled.
 * 
 * step 2: Initialize the component index a = 0, then for each compoent do:
 * 
 * step 3: a = a +1
 * 
 * step 4: t-start = column in X with largest Variance
 *  
 * 
 * step 5: P = t'X    ||P|| = 1 -> Loadings
 * 
 * step 5.1:  additional NIPALS loop that decomposes P into u_1r_1'+ u_2r_2'. This loop should the be run through twice.
 *     a) h = 0; u-start is the "largest" column in P;
 *     b) r' = u'P/u'u
 *     c) norm to ||r|| = 1
 *     d) u = Pr
 *     e) h = h + 1; if h < 2 then back to "b"
 *     f) set P = ur' and procede with "step 6"
 *     
 * step 6: t = X..P/||P||^2
 * 
 * step 7: check convergente: 
 *     if d = (t_new - t_old)'(t_new-t_old)/N*t_new't_new) > 10^-10
 *        back to "step 5"
 * 
 * step 8: After convergence, residuals: E = X - t <scalar_product> 
 * Then set X = E and extimante the next component by returning to step 3
 * 
 */

void UPCA(array *X_, size_t npc, size_t autoscaling, UPCAMODEL *m, ssignal *s)
{
  if(CheckArray(X_) == 0){
    size_t pc, i, j, k, order, column;
    double a, ss, min, max;
    array *X;
    dvector *colvar;
    dvector *t_old; /* row vector with size X->m->row */
    dvector *t_new;
    dvector *t_diff;
    dvector *eval; /* t't */
    
    matrix *P;
    dvector *r;
    dvector *u;
    
    dvector *tmpv;
    
    /* step 1 */    
    NewArray(&X, X_->order);

    for(k = 0; k < X_->order; k++){
      NewArrayMatrix(&X, k, X_->m[k]->row, X_->m[k]->col);
      
    /* For Wold Geladi test enable this 
    for(i = 0; i < X_->m[k]->row; i++)
      for(j = 0; j < X_->m[k]->col; j++)
        setArrayValue(X, k, i, j, getArrayValue(X_, k, i, j));
      */
    }
        
    /* For Wold Geladi test disable this  */
     
    for(k = 0; k < X_->order; k++){
      initDVector(&tmpv);
      MatrixColAverage(X_->m[k], &tmpv);
      MatrixAppendCol(&(m->colaverage), tmpv);
      DelDVector(&tmpv);
      
      for(j = 0; j < X_->m[k]->col; j++){
        if(X_->m[k]->row > 1){
          for(i = 0; i < X_->m[k]->row; i++){
            setArrayValue(X, k, i, j, getArrayValue(X_, k, i, j) - getMatrixValue(m->colaverage, j, k));
          }
        }
        else{
          for(i = 0; i < X_->m[k]->row; i++){
            setArrayValue(X, k, i, j, getArrayValue(X_, k, i, j));
          }
        }
      }

    }
    
    /* AUTOSCALING */
    if(autoscaling > 0){
      if(autoscaling == 1){
        for(k = 0; k < X_->order; k++){
          initDVector(&tmpv);
          MatrixColSDEV(X_->m[k], &tmpv);
          MatrixAppendCol(&(m->colscaling), tmpv);
          DelDVector(&tmpv);
        }
      }
      else if(autoscaling == 2){
        for(k = 0; k < X_->order; k++){
          initDVector(&tmpv);
          MatrixColRMS(X_->m[k], &tmpv);
          MatrixAppendCol(&(m->colscaling), tmpv);
          DelDVector(&tmpv);
        }
      }
      else if(autoscaling == 3){ /* PARETO Autoscaling */
        for(k = 0; k < X_->order; k++){
          initDVector(&tmpv);
          MatrixColSDEV(X_->m[k], &tmpv);
          for(i = 0; i < tmpv->size; i++){
            setDVectorValue(tmpv, i, sqrt(getDVectorValue(tmpv, i)));
          }
          MatrixAppendCol(&(m->colscaling), tmpv);
          DelDVector(&tmpv);
        }
      }
      else if(autoscaling == 4){ /* Range Scaling */
        for(k = 0; k < X_->order; k++){
          initDVector(&tmpv);
          for(i = 0; i < X_->m[k]->col; i++){
            MatrixColumnMinMax(X_->m[k], i, &min, &max);
            DVectorAppend(&tmpv, (max-min));
          }
          MatrixAppendCol(&(m->colscaling), tmpv);
          DelDVector(&tmpv);
        }
      }
      else if(autoscaling == 5){ /* Level Scaling  */
        MatrixCopy(m->colaverage, &m->colscaling);
      }
      else{ /*no autoscaling so divide all for 1 */
        for(k = 0; k < X_->order; k++){
          initDVector(&tmpv);
          for(i = 0; i < X_->m[k]->col; i++){
            DVectorAppend(&tmpv, 1.0);
          }
          MatrixAppendCol(&(m->colscaling), tmpv);
          DelDVector(&tmpv);
        }
      }

      for(k = 0; k < X->order; k++){
        for(j = 0; j < X->m[k]->col; j++){
          if(getMatrixValue(m->colscaling, j, k) == 0){
            for(i = 0; i< X->m[k]->row; i++){
              setArrayValue(X, k, i, j, 0.f);
            }
          }
          else{
            for(i = 0; i < X->m[k]->row; i++){
              setArrayValue(X, k, i, j, getArrayValue(X, k, i, j)/getMatrixValue(m->colscaling, j, k));
            }
          }
        }
      }
    }
    
   /* END Disable */
   
    if(npc > X->m[0]->col)
      npc = X->m[0]->col;
    
    ss = 0.f;
    for(k = 0; k < X->order; k++)
      for(i = 0; i < X->m[k]->row; i++)
        for(j = 0; j < X->m[k]->col; j++)
          ss += square(getArrayValue(X, k, i, j));
    
    #ifdef DEBUG
    puts("Processed Array");
    PrintArray(X);
    #endif
    
    NewDVector(&t_old, X->m[0]->row); /* the object number is the same for all the matrix in array */
    NewDVector(&t_new, X->m[0]->row); /* the object number is the same for all the matrix in array */
    NewMatrix(&P, X->m[0]->col, X->order); /*Check that all the matrix of array have the same row and column */
    NewDVector(&eval, npc);

    /* step 2 and 3 */
    for(pc = 0; pc < npc; pc++){
      if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
        break;
      }
      else{
        /* step 4 select the column vector t with the largest column variance */
        order = column = 0;
        for(i = 0; i < X->order; i++){
          initDVector(&colvar);
          MatrixColVar(X->m[i], &colvar);
          
          #ifdef DEBUG      
          puts("Column Variance for chosing the best t-start column vector");
          PrintDVector(colvar);
          #endif
          
          /* check the column in array with the largest column variance */
          k = 0;
          for(j = 1; j < X->m[i]->col; j++){
            if(getDVectorValue(colvar, j) > getDVectorValue(colvar, k)){
              k = j;
            }
            else{
              continue;
            }
          }
          
          if(i > 0){
            if(getDVectorValue(colvar, k) > a){
              a = getDVectorValue(colvar, k);
              order = i;
              column = k;
            }
            else{
              DelDVector(&colvar);
              continue;
            }
          }
          else{
            a = getDVectorValue(colvar, k);
            order = i;
            column = k;
          }
          
          DelDVector(&colvar);
        }
      
        
        /* copy the best column vector into t */
        for(i = 0; i < X->m[order]->row; i++){
          setDVectorValue(t_old, i, getArrayValue(X, order, i, column));
        }
        
        #ifdef DEBUG
        puts("The best t-start vector");
        PrintDVector(t_old);
        #endif
        
        while(1){
          /*step 5*/
          MatrixSet(P, 0.f);
          DvectorArrayDotProduct(X, t_old, P);
          
          MatrixNorm(P, P);
          
          #ifdef DEBUG
          puts("P Loadings Matrix");
          PrintMatrix(P);
          #endif
          
          /*step 5.1*/
          #ifdef DEBUG
          puts("Nipals SHUNT");
          #endif
          NewDVector(&r, P->col);
          NewDVector(&u, P->row);
          
          /*u-start is the larget column in P*/
          initDVector(&colvar);
          MatrixColAverage(P, &colvar);
          j = 0;
          for(i = 1; i < colvar->size; i++){
            if(getDVectorValue(colvar, i) > getDVectorValue(colvar, j))
              j = i;
            else
              continue;
          }
          DelDVector(&colvar);
          
          /*copy the best column vector from P to u*/
          for(i = 0; i < P->row; i++)
            setDVectorValue(u, i, getMatrixValue(P, i, j));
          
          #ifdef DEBUG
          puts("U-start");
          PrintDVector(u);
          #endif
          
          for(i = 0; i < 2; i++){ /* for the simplest case */
            
            /* 2.  p' := (t'X)/(t't)        Project the matrix E onto t in order to find the corresponding loading p
              * 3.  p' := p'/|p'|            Normalize the loading vector p to length 1
              * 4.  t_old := t
              *     t := (Xp)/(p'p)          Store the score vector t into uold and project the matrix E onto p in order to find corresponding score vector t
              */           

            
            /* r' = u'P/u'u */
            DVectorSet(r, 0.f);
            DVectorMatrixDotProduct(P, u, r);
            a = DVectorDVectorDotProd(u, u);
            
            for(j = 0; j < r->size; j++)
              setDVectorValue(r, j, getDVectorValue(r, j) / a);
            
            /*||r|| = 1; */
            DVectNorm(r, r);
            
            #ifdef DEBUG
            puts("r vector");
            PrintDVector(r);
            #endif
            
            /*u = Pr */
            DVectorSet(u, 0.f); /* reset the u vector*/
            MatrixDVectorDotProduct(P, r, u);
            
            /*
            a = DVectorDVectorDotProd(r, r);
            
            for(j = 0; j < u->size; j++)
              setDVectorValue(u, j, getDVectorValue(u, j) / a);
            
            */
            #ifdef DEBUG
            puts("u vector");
            PrintDVector(u);
            #endif
          }
          
          /*step 5= f  set P = ur' */
          MatrixSet(P, 0.f);
          RowColOuterProduct(u, r, P);
          
          #ifdef DEBUG
          puts("Recomputed P Loadings Matrix from ur'");
          PrintMatrix(P);
          #endif
          
          
          DelDVector(&r);
          DelDVector(&u);
          
          /* step 7 
          #ifdef DEBUG
          puts("old t score vector");
          PrintDVector(t_old);
          #endif
          */
          
          DVectorSet(t_new, 0.f);
          ArrayMatrixDotProduct(X, P, t_new);
          
          a = Matrixnorm(P);
          a = square(a);
          
          for(i = 0; i < t_new->size; i++)
            setDVectorValue(t_new, i, getDVectorValue(t_new, i)/a);

          
          #ifdef DEBUG
          puts("new t score vector");
          PrintDVector(t_new);
          #endif
          
          /*step 8 */
          initDVector(&t_diff);
          DVectorDVectorDiff(t_new, t_old, &t_diff);
          
          a = DVectorDVectorDotProd(t_diff, t_diff);
          
          #ifdef DEBUG
          puts("t_old - t_new vector");
          PrintDVector(t_diff);
          #endif
          
          DelDVector(&t_diff);
          
          
          /*Check Distance */
  /*         if(a/(t_new->size*DVectorDVectorDotProd(t_new, t_new)) < 1e-10){*/
          if(a/(t_new->size*DVectorDVectorDotProd(t_new, t_new)) < EPSILON){
            /* storing results */
            MatrixAppendCol(&(m->scores), t_new);
            ArrayAppendMatrix(&(m->loadings), P);
            setDVectorValue(eval, pc, DVectorDVectorDotProd(t_new, t_new));
            break;
          }
          else if(_isnan_(DVectorDVectorDotProd(t_new, t_new))){
            fprintf(stderr, "UPCA Error! The Solver Engine was Unable to Converge! Please Check your data.\n");
            fflush(stderr);
            /* abort(); */
          }
          else{
            for(i = 0; i < t_new->size; i++)
              setDVectorValue(t_old, i, getDVectorValue(t_new, i));
            continue;
          }
        }
        
        
        /* step 8 remove computed component */   
        for(k = 0; k < X->order; k++){
          for(i = 0; i < X->m[k]->row; i++){
            for(j = 0; j < X->m[k]->col; j++){
              setArrayValue(X, k, i, j, getArrayValue(X, k, i, j) - (getDVectorValue(t_new, i)*getMatrixValue(P, j, k)));
            }
          }    
        }

        /*Reset all*/
        MatrixSet(P, 0.f);
        DVectorSet(t_new, 0.f);
        DVectorSet(t_old, 0.f);
      }
    }
    
    calcVarExpressed(ss, eval, &(m->varexp));
    
    DelDVector(&eval);
    DelMatrix(&P);
    DelDVector(&t_old);
    DelDVector(&t_new);
    DelArray(&X);
  }
  else{
    fprintf(stderr, "Error!! Unable to compute Multi Way PCA. The (row,col) size of array submatrix are different!\n");
    /* abort(); */
  }
}

void UPCAScorePredictor(array *X_, UPCAMODEL *model, size_t npc, matrix **pscores)
{
  size_t pc, i, j , k;
  array *X;
  dvector *t;
  
  double a;

  NewArray(&X, X_->order);
  for(k = 0; k < X_->order; k++){
    NewArrayMatrix(&X, k, X_->m[k]->row, X_->m[k]->col);
  }
  
  NewDVector(&t, X->m[0]->row);

  for(k = 0; k < X_->order; k++){
    for(j = 0; j < X_->m[k]->col; j++){
      /* Mean Centering */
      if(model->colaverage->row > 0){
        for(i = 0; i < X_->m[k]->row; i++){
          setArrayValue(X, k, i, j, getArrayValue(X_, k, i, j) -  getMatrixValue(model->colaverage, j, k));
        }
      }
      else{
        for(i = 0; i < X_->m[k]->row; i++){
          setArrayValue(X, k, i, j, getArrayValue(X_, k, i, j));
        }
      }
      
      if(model->colscaling->row > 0){
        /* Scaling to Column SDEV */
        for(i = 0; i < X->m[k]->row; i++){
          setArrayValue(X, k, i, j, getArrayValue(X, k, i, j) / getMatrixValue(model->colscaling, j,  k));
        }
      }
    }
  } 
  
  if(model->loadings->order < npc)
    npc = model->loadings->order;
  
  for(pc = 0; pc < npc; pc++){
    DVectorSet(t, 0.f);
    ArrayMatrixDotProduct(X, model->loadings->m[pc], t);
    
    
     
    a = Matrixnorm(model->loadings->m[pc]);
    a = square(a);
        
    for(i = 0; i < t->size; i++)
      setDVectorValue(t, i, getDVectorValue(t, i)/a);

    MatrixAppendCol(pscores, t);
    

    for(k = 0; k < X->order; k++){
      for(i = 0; i < X->m[k]->row; i++){
        for(j = 0; j < X->m[k]->col; j++){
          setArrayValue(X, k, i, j, getArrayValue(X, k, i, j) - (getDVectorValue(t, i)*getMatrixValue(model->loadings->m[pc], j, k)));
        }
      }    
    }
  }
 
  DelDVector(&t);
  DelArray(&X);
}


/*
 * X = TP' + ....
 * 
 */
void UPCAIndVarPredictor(matrix *T, array *P, matrix *colaverage, matrix *colscaling,  size_t npc, array **X)
{
  size_t pc, i, j, k;
  
  /*Allocate the output array*/
  for(i = 0; i < P->m[0]->col; i++)
    AddArrayMatrix(X, T->row, P->m[0]->row);
    /*AddArrayMatrix(X, T->row, P->m[0]->col);*/
  
  for(pc = 0; pc < npc; pc++){ 
    for(k = 0; k < (*X)->order; k++){
      for(i = 0; i < (*X)->m[k]->row; i++){
        for(j = 0; j < (*X)->m[k]->col; j++){
          setArrayValue((*X), k, i, j, getArrayValue((*X), k, i, j) + (getMatrixValue(T, i, pc) * getArrayValue(P, pc, j, k)));
        }
      }    
    }
  }
  
  /*WARNING da testare se tornano i valori predetti*/
  
 /*
  * xcolaverage->row and xcolscaling->row are the number of variable for each submatrix in X_, so these are equal to X_->m[i]->col
  * 
  * xcolaverage->col and xcolscaling->col are the number of order present in X_, so these are equal to X_->order
  */
  
  if(colaverage->row > 0 && colaverage->col > 0){
    for(k = 0; k < (*X)->order; k++){
      for(i = 0; i < (*X)->m[k]->row; i++){
        for(j = 0; j < (*X)->m[k]->col; j++){
          if(colscaling->row > 0 && colscaling->col > 0){
            setArrayValue((*X), k, i, j, (getArrayValue((*X), k, i, j) * getMatrixValue(colscaling, j, k)) + getMatrixValue(colaverage, j, k));
          }
          else{
            setArrayValue((*X), k, i, j, getArrayValue((*X), k, i, j) + getMatrixValue(colaverage, j, k));
          }
        }
      }
    }
  }
  
}


void PrintUPCAModel(UPCAMODEL* m)
{
  puts("Scores");
  PrintMatrix(m->scores);
  
  puts("Loadings");
  PrintArray(m->loadings);
  
  puts("Variance Explained");
  PrintDVector(m->varexp);
}

