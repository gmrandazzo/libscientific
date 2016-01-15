#include "statistic.h"
#include "numeric.h"
#include "memwrapper.h"
#include <math.h>

/*
 * TP = True PositivePredictedValue
 * FN = False Negative
 * 
 * Sensitivity = TP / TP + FN
 */
void Sensitivity(dvector *dtp, double thmin, double thmax, double thstep, matrix** s)
{
  size_t i, row, tp, fn;
  double dx;
  
  ResizeMatrix(s, (size_t) ceil((thmax-thmin)/thstep), 2); /* two column: x and y */
  
  row = 0;
  dx = thstep;
  
  for(row = 0; row < (*s)->row; row++){
    tp = fn = 0;
    
    for(i = 0; i < dtp->size; i++){
      if(getDVectorValue(dtp, i) < dx)
        tp++;
      else if(FLOAT_EQ(getDVectorValue(dtp, i), dx, EPSILON))
        tp++;
      else
        fn++;
    }
    
    if(tp == 0){
      setMatrixValue((*s), row, 0, dx); setMatrixValue((*s), row, 1, 0);
    }
    else{
      setMatrixValue((*s), row, 0, dx); setMatrixValue((*s), row, 1, (double)tp/(double)(tp+fn));
    }
    
    dx+=thstep;
  }

}

/*
 * TP = True Positive
 * FP = False Positive
 * 
 * PPV = TP / TP + FP
 */
void PositivePredictedValue(dvector* dtp, dvector* dtn, double thmin, double thmax, double thstep, matrix** p)
{
  size_t i, row, tp, fp;
  double dx;
  
  ResizeMatrix(p, (size_t) ceil((thmax-thmin)/thstep), 2); /* two column: x and y */
  
  row = 0;
  dx = thstep;

  for(row = 0; row < (*p)->row; row++){
    tp = fp = 0;
    
    for(i = 0; i < dtp->size; i++){
      if(getDVectorValue(dtp, i) < dx)
        tp++;
      else if(FLOAT_EQ(getDVectorValue(dtp, i), dx, EPSILON))
        tp++;
    }
    
    for(i = 0; i < dtn->size; i++){
      if(getDVectorValue(dtn, i) < dx)
        fp++;
      else if(FLOAT_EQ(getDVectorValue(dtn, i), dx, EPSILON))
        fp++;
    }
    
    if(tp == 0){
      setMatrixValue((*p), row, 0, dx); setMatrixValue((*p), row, 1, 0);
    }
    else{
      setMatrixValue((*p), row, 0, dx); setMatrixValue((*p), row, 1, (double)tp/(double)(tp+fp));
    }

    dx+=thstep;
  }
}


/*
 * This function Code the matrix by using the 0, 1 scale
 * according to the function
 * 
 * x = (val - mid) / step
 */
void MatrixCode(matrix* inmx, matrix* outmx)
{
  size_t i, j;
  double min, nmin, max, nmax, mid, step;
  ResizeMatrix(&outmx, inmx->row, inmx->col);
  
  for(j = 0; j < inmx->col; j++){
    min = max = nmin = nmax = getMatrixValue(inmx, 0, j);
    for(i = 0; i < inmx->row; i++){
      
      if(getMatrixValue(inmx, i, j) < min){
        min = getMatrixValue(inmx, i, j);
      }
      
      if(getMatrixValue(inmx, i, j) > max){
        max = getMatrixValue(inmx, i, j);
      }
    }
    
    mid = min + ((max - min) / 2);
    
    for(i = 0; i < inmx->row; i++){
      if(getMatrixValue(inmx, i, j) > min && getMatrixValue(inmx, i, j) < mid && getMatrixValue(inmx, i, j) < nmin){
        nmin = getMatrixValue(inmx, i, j);
      }
      
      if(getMatrixValue(inmx, i, j) < max && getMatrixValue(inmx, i, j) > mid &&  getMatrixValue(inmx, i, j) > nmax){
        nmax = getMatrixValue(inmx, i, j);
      }
    }
    
    if(nmax != nmin){
      step = (nmax- nmin) / 2;
      mid = nmin + step;
    }
    else{
      step = (max - min) / 2;
      mid = min + step;
    }
    
    for(i = 0; i < inmx->row; i++){
      setMatrixValue(outmx, i, j, (getMatrixValue(inmx, i, j) - mid) / step);
    }
  }
}

/*
 * This function expand the matrix to Bifactorial cross product
 * a, b, c
 * a*a, b*b, c*c,
 * a*b, a*c, b*c
 */
void BifactorialMatrixExpansion(matrix* inmx, matrix* outmx)
{
  size_t i, j, k;
  size_t n = inmx->col, n_k, comb = 1;
  
  if(inmx->col-2 != 0){
    /*n! / (n-k)!*k!*/
    n = Factorial(inmx->col);
    n_k = Factorial(inmx->col-2);
    comb = (n / (n_k *2));
  }
  
  ResizeMatrix(&outmx, inmx->row, 2*inmx->col + comb + 1);
  
  
  /* Generate the matrix */
  for(i = 0; i < inmx->row; i++){
    setMatrixValue(outmx, i, 0, 1);
  }
  
  for(i = 0; i < inmx->row; i++){
    for(j = 0; j < inmx->col; j++){
      setMatrixValue(outmx, i, j+1, getMatrixValue(inmx, i, j));
      setMatrixValue(outmx, i, inmx->col+1+j+comb, square(getMatrixValue(inmx, i, j)));
    }
  }
  
  for(i = 0; i < inmx->row; i++){
    for(j = 0; j < inmx->col; j++){
      for(k = j+1; k < inmx->col; k++){
        setMatrixValue(outmx, i, inmx->col+j+k, getMatrixValue(inmx, i, j)*getMatrixValue(inmx, i, k));
      }
    }
  }

}

void YatesVarEffect(matrix* mx, dvector* veff)
{

}
