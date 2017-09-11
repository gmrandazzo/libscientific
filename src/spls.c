#include "memwrapper.h"
#include "spls.h"
#include "pls.h" /* Use LVGet to calculate the latent variable according NIPALS */
#include "pca.h" /*Using: MatrixAutoScaling(); and calcVarExpressed(); */
#include "numeric.h" /* Using:  if(FLOAT_EQ(NumOne, NumTwo));*/
#include "metricspace.h"
#include <math.h>
#include <pthread.h>

void SetNBlocks(BLOCKS **b, size_t n)
{
  size_t i;
  (*b) = xmalloc(sizeof(BLOCKS)*n);
  for(i = 0; i < n; i++){
    (*b)[i].from = 0;
    (*b)[i].to = 0;
  }
}

void DelBlocks( BLOCKS **b)
{
  xfree((*b));
}

void NewSPLSModel(SPLSMODEL** m)
{
  (*m) = xmalloc(sizeof(PLSMODEL));
  initMatrix(&(*m)->xscores);
  initMatrix(&(*m)->xloadings);
  initMatrix(&(*m)->xweights);
  initMatrix(&(*m)->yscores);
  initMatrix(&(*m)->yloadings);
  initDVector(&(*m)->b);
  initDVector(&(*m)->xvarexp);
  /*initDVector(&(*m)->yvarexp);*/
  initDVector(&(*m)->xcolaverage);
  initDVector(&(*m)->xcolscaling);
  initDVector(&(*m)->ycolaverage);
  initDVector(&(*m)->ycolscaling);
  initMatrix(&(*m)->r2y_model);
  initMatrix(&(*m)->r2y_validation);
  initMatrix(&(*m)->recalc_residuals);
  initMatrix(&(*m)->q2y);
  initMatrix(&(*m)->sdep);
  initMatrix(&(*m)->sdec);
  initMatrix(&(*m)->bias);
  initMatrix(&(*m)->recalculated_y);
  initMatrix(&(*m)->predicted_y);
  initMatrix(&(*m)->pred_residuals);
  initMatrix(&(*m)->r2q2scrambling);
  initMatrix(&(*m)->q2_sample_validation);
  initMatrix(&(*m)->sdep_sample_validation);
  initMatrix(&(*m)->q2_sample_validation_surface);
  initMatrix(&(*m)->sdep_sample_validation_surface);
}

void DelSPLSModel(SPLSMODEL** m)
{
  DelMatrix(&(*m)->xscores);
  DelMatrix(&(*m)->xloadings);
  DelMatrix(&(*m)->xweights);
  DelMatrix(&(*m)->yscores);
  DelMatrix(&(*m)->yloadings);
  DelDVector(&(*m)->b);
  DelDVector(&(*m)->xvarexp);
  /*DelDVector(&(*m)->yvarexp); */
  DelDVector(&(*m)->xcolaverage);
  DelDVector(&(*m)->xcolscaling);
  DelDVector(&(*m)->ycolaverage);
  DelDVector(&(*m)->ycolscaling);
  DelMatrix(&(*m)->r2y_model);
  DelMatrix(&(*m)->r2y_validation);
  DelMatrix(&(*m)->recalc_residuals);
  DelMatrix(&(*m)->q2y);
  DelMatrix(&(*m)->sdep);
  DelMatrix(&(*m)->sdec);
  DelMatrix(&(*m)->bias);
  DelMatrix(&(*m)->recalculated_y);
  DelMatrix(&(*m)->predicted_y);
  DelMatrix(&(*m)->pred_residuals);
  DelMatrix(&(*m)->r2q2scrambling);
  DelMatrix(&(*m)->q2_sample_validation);
  DelMatrix(&(*m)->sdep_sample_validation);
  DelMatrix(&(*m)->q2_sample_validation_surface);
  DelMatrix(&(*m)->sdep_sample_validation_surface);
  xfree((*m));
}


/* Algorithm S-PLS
 * A serial extension of multiblock PLS
 * Anders Berglund, Svante Wold
 * May 1999
 * DOI: 10.1002/(SICI)1099-128X(199905/08)13:3/4<461::AID-CEM555>3.0.CO;2-B
 *
 *
 * X1 = T1 x P1' + E1
 * X2 = T2 x P2' + E2
 * Y = T1 x C1' + T2 x C2' + F
 *
 * Algorithm steps for two blocks:
 *
 * 1) Set F2 = Y
 * 2) Calculate the first PLS with the X and F2
 * 3) Calculate the residuals F1 = Y - T1 x C1'
 * 4) Calculate the second PLS model with X2 and F1
 * 5) Calculate the residuals F2 = Y - T2 x C2'
 * 6) Check for residuals congergence; if not converged, return to step 2
 *
 * Giuseppe Marco Randazzo <gmrandazzo@gmail.com> 8 Mar 2017
 */
void SPLS(matrix *mx, matrix *my, BLOCKS *b, size_t nlv, size_t xautoscaling, size_t yautoscaling, SPLSMODEL *m, ssignal *s)
{
  size_t i, j, k, n;
  dvector *t, *p, *w, *u, *q, *xecal, *yeval;

  /*Create the n X matrix according the b blocks */
  n = sizeof(b)/sizeof(BLOCKS);
  matrix **x = xmalloc(sizeof(matrix*)*n);
  matrix **y = xmalloc(sizeof(matrix*)*n);
  matrix *f;
  NewMatrix(&f, my->row, my->col);

  /*initialize the n matrix by splitting in  scaling and centering */
  for(i = 0; i < n; i++){
    NewMatrix(&x[i], mx->row, (b[i].to - b[i].from));
    for(j = 0; j < mx->row; j++)
      for(k = b[i].from; k < b[i].to; k++)
        x[i]->data[j][k-b[i].from] = mx->data[j][k];
  }

  NewDVector(&t, X->row);
  NewDVector(&p, X->col);
  NewDVector(&w, X->col);

  NewDVector(&u, Y->row);
  NewDVector(&q, Y->col);

  NewDVector(&xeval, nlv);
  NewDVector(&yeval, nlv);



  while(1){

    for(i = 0; i < n; i++){
      double bcoef = 0.f;
      LVCalc(&x[i], &f[i], &t, &u, &p, &q, &w, &bcoef);
      
    }
  }



  /* free the memory */
  for(i = 0; i < n; i++){
    DelMatrix(&x[i]);
    DelMatrix(&y[i]);
  }
  xfree(&y);
  xfree(&x);
  DelDVector(&t);
  DelDVector(&p);
  DelDVector(&w);
  DelDVector(&u);
  DelDVector(&q);
  DelDVector(&xeval);
  DelDVector(&yeval);
}
