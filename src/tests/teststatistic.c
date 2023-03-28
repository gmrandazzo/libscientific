#include <stdio.h>
#include <math.h>
#include "numeric.h"
#include "statistic.h"
#include <time.h>

/*
 * 
 double R2(dvector *ytrue, dvector *ypred);
 double MAE(dvector *ytrue, dvector *ypred);
 double MSE(dvector *ytrue, dvector *ypred);
 double RMSE(dvector *ytrue, dvector *ypred);
 double BIAS(dvector *ytrue, dvector *ypred);$
 void Sensitivity(dvector *dtp,
                 double thmin,
                 double thmax,
                 double thstep,
                 matrix *s);
void PositivePredictedValue(dvector *dtp,
        dvector *dtn,
        double thmin,
        double thmax,
        double thstep,
        matrix *p);
void MatrixCode(matrix *inmx, matrix *outmx);
void BifactorialMatrixExpansion(matrix* inmx, matrix* outmx);
void YatesVarEffect(matrix *mx, dvector *veff);

*/
void test1()
{
  puts("Test1: ROC and Precision-Recall test.");
  /* ROC curve test */
  dvector *y_score, *y_true;
  matrix *roc, *pr;
  double auc, ap;
  NewDVector(&y_score, 13);
  NewDVector(&y_true, 13);

  y_score->data[0] = 0.1;
  y_score->data[1] = 0.2;
  y_score->data[2] = 0.3;
  y_score->data[3] = 0.4;
  y_score->data[4] = 0.5;
  y_score->data[5] = 0.6;
  y_score->data[6] = 0.7;
  y_score->data[7] = 0.8;
  y_score->data[8] = 0.9;
  y_score->data[9] = 1.0;
  y_score->data[10] = 2.0;
  y_score->data[11] = 3.0;
  y_score->data[12] = 3.4;


  y_true->data[0] = 0;
  y_true->data[1] = 0;
  y_true->data[2] = 0;
  y_true->data[3] = 0;
  y_true->data[4] = 0;
  y_true->data[5] = 1;
  y_true->data[6] = 1;
  y_true->data[7] = 1;
  y_true->data[8] = 1;
  y_true->data[9] = 0;
  y_true->data[10] = 1;
  y_true->data[11] = 1;
  y_true->data[12] = 1;

  initMatrix(&roc);
  ROC(y_true, y_score,  roc, &auc);
  PrintMatrix(roc);
  if(FLOAT_EQ(auc, 0.904762, 1e-6)){
    printf("AUC OK!\n");
  }else{
    printf("AUC ERROR!\n");
  }
  printf("AUC: %f\n", auc);

  initMatrix(&pr);
  PrecisionRecall(y_true, y_score, pr, &ap);
  PrintMatrix(pr);

  if(FLOAT_EQ(ap, 0.900425, 1e-6)){
    printf("AVERAGE PRECISION-RECALL OK!\n");
  }else{
    printf("AVERAGE PRECISION-RECALL ERROR!\n");
  }
  printf("AVERAGE PRECISION-RECALL: %f\n", ap);

  DelMatrix(&pr);
  DelMatrix(&roc);
  DelDVector(&y_true);
  DelDVector(&y_score);
}


int main(void)
{
  test1();
}
