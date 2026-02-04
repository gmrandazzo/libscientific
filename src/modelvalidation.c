/* Provides tools for model validation and error estimation.
 * Copyright (C) 2016-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include "modelvalidation.h"
#include "memwrapper.h"
#include "mlr.h"
#include "pls.h"
#include "pca.h" /*Using: MatrixAutoScaling(); and calcVarExpressed(); */
#include "epls.h"
#include "lda.h"
#include "numeric.h" /* Using:  if(FLOAT_EQ(NumOne, NumTwo));*/
#include "metricspace.h"
#include "statistic.h"
#include <math.h>
#include <pthread.h>
#include <stdarg.h>


ValidationArg initValidationArg()
{
  return (ValidationArg) {.vtype = BootstrapRGCV,
                          .rgcv_group = 5,
                          .rgcv_iterations = 5};
}

MODELINPUT initModelInput()
{
  return (MODELINPUT) {.mx = NULL,
                       .my = NULL,
                       .nlv = 0,
                       .xautoscaling = 0,
                       .yautoscaling = 0
                      };
}
void random_kfold_group_generator(matrix *gid,
                                  size_t ngroups,
                                  size_t nobj,
                                  unsigned int *srand_init)
{
  /*equal distribution of objects*/
  ResizeMatrix(gid, ngroups, (size_t)ceil(nobj/(double)ngroups));
  MatrixSet(gid, -1); /*default is -1*/
  size_t i;
  size_t j;
  size_t k = 0;
  size_t n;
  srand_((uint32_t)(*srand_init));
  for(i = 0; i <  gid->row; i++){
    for(j = 0; j <  gid->col; j++){
      do{
        n = randInt(0, nobj);
      } while(ValInMatrix(gid, n) == 1 && k < nobj);
      if(k < nobj){
        gid->data[i][j] = n;
        k++;
      }
      else
        continue;
    }
  }
}

void kfold_group_train_test_split(matrix *x,
                                  matrix *y,
                                  matrix *gid,
                                  size_t group_id,
                                  matrix *x_train,
                                  matrix *y_train,
                                  matrix *x_test,
                                  matrix *y_test)
{
  /* Estimate how many objects are utilised to build the model and how many to predict */
  size_t i, j, k, l, n;
  size_t m_obj = 0;
  size_t p_obj = 0;
  for(i = 0; i < gid->row; i++){
    if(i != group_id){
      for(j = 0; j < gid->col; j++){
        if((int)gid->data[i][j] != -1)
          m_obj++;
        else
          continue;
      }
    }
    else{
      for(j = 0; j < gid->col; j++){
        if((int)gid->data[i][j] != -1)
          p_obj++;
        else
          continue;
      }
    }
  }

  /*Allocate the train matrix */
  ResizeMatrix(x_train, m_obj, x->col);
  ResizeMatrix(y_train, m_obj, y->col);

  /*Allocate the test matrix */
  ResizeMatrix(x_test, p_obj, x->col);
  ResizeMatrix(y_test, p_obj, y->col);
  /* copy the train and test values */

  for(i = 0, k = 0, l = 0; i < gid->row; i++){
    if(i != group_id){
      for(j = 0; j < gid->col; j++){
        int a =  (int)gid->data[i][j]; /* get the row index */
        if(a != -1){
          for(n = 0; n < x->col; n++){
            x_train->data[k][n] = x->data[a][n];
          }
          for(n = 0; n < y->col; n++){
            y_train->data[k][n] = y->data[a][n];
          }
          k++;
        }
        else{
          continue;
        }
      }
    }
    else{
      for(j = 0; j < gid->col; j++){
        int a = (int)gid->data[i][j];
        if(a != -1){
          for(n = 0; n < x->col; n++){
            x_test->data[l][n] = x->data[a][n];
          }
          for(n = 0; n < y->col; n++){
            y_test->data[l][n] = y->data[a][n];
          }
          l++;
        }
        else{
          continue;
        }
      }
    }
  }
}

void train_test_split(matrix *x,
                      matrix *y,
                      double testsize,
                      matrix *x_train,
                      matrix *y_train,
                      matrix *x_test,
                      matrix *y_test,
                      uivector *testids,
                      unsigned int *srand_init)
{
  size_t i, j, n, testsize_;
  if(testsize > 1.0 || testsize < 0.){
    testsize_ = (int)ceil(0.2*x->row); /* default value */
  }
  else{
    testsize_ = (int)ceil(testsize*x->row);
  }

  ResizeMatrix(x_train, (x->row-testsize_), x->col);
  ResizeMatrix(y_train, (y->row-testsize_), y->col);
  ResizeMatrix(x_test, testsize_, x->col);
  ResizeMatrix(y_test, testsize_, y->col);

  /*Fill test set*/
  srand_((uint32_t)(*srand_init));
  for(i = 0; i < testsize_; i++){
    do{
      n = randInt(0, x->row);
    } while(UIVectorHasValue(testids, n) == 0);
    UIVectorAppend(testids, n);

    for(j = 0; j < x->col; j++){
      x_test->data[i][j] = x->data[n][j];
    }
    for(j = 0; j < y->col; j++){
      y_test->data[i][j] = y->data[n][j];
    }
  }
  /* Fill training set */
  for(i = 0, n = 0; i < x->row; i++){
    if(UIVectorHasValue(testids, i) == 1){
      for(j = 0; j < x->col; j++){
        x_train->data[n][j] = x->data[i][j];
      }
      for(j = 0; j < y->col; j++){
        y_train->data[n][j] = y->data[i][j];
      }
      n++;
    }
    else{
      continue;
    }
  }
}

typedef struct{
  ELearningParameters eparm;
  CombinationRule crule;
  matrix *mx, *my; /*INPUT*/
  matrix *predicted_y;  /*OUPUT*/
  uivector *predictioncounter; /*OUPUT*/
  size_t xautoscaling, yautoscaling, nlv, group; /*INPUT*/
  unsigned int srand_init;
} rgcv_th_arg;

void *PLSRandomGroupCVModel(void *arg_)
{
  size_t j, k, n, g;
  rgcv_th_arg *arg;
  matrix *gid; /* randomization and storing id for each random group into a matrix */

  /* Matrix for compute the PLS models for groups */
  matrix *x_train;
  matrix *y_train;
  PLSMODEL *subm;

  /*matrix to predict*/
  matrix *x_test;
  matrix *y_test;
  matrix *y_test_predicted;

  arg = (rgcv_th_arg*) arg_;

  /* step 1 generate the random groups */
  initMatrix(&gid);
  random_kfold_group_generator(gid, arg->group, arg->mx->row, &arg->srand_init);

  /*
  puts("Gid Matrix");
  PrintMatrix(gid);
  */

  /* step 2 */
  for(g = 0; g < gid->row; g++){ /*For aeach group */
    /* Estimate how many objects are inside the group "g" to predict and outside the group "g" to build the model  */
    initMatrix(&x_train);
    initMatrix(&y_train);
    initMatrix(&x_test);
    initMatrix(&y_test); /* unused variable here .... */

    kfold_group_train_test_split(arg->mx,
                                 arg->my,
                                 gid,
                                 g,
                                 x_train,
                                 y_train,
                                 x_test,
                                 y_test);

    NewPLSModel(&subm);

    PLS(x_train, y_train, arg->nlv, arg->xautoscaling, arg->yautoscaling, subm, NULL);
    initMatrix(&y_test_predicted);
    PLSYPredictorAllLV(x_test, subm, NULL, y_test_predicted);

    for(j = 0, k = 0; j < gid->col; j++){
      int a = (int)gid->data[g][j]; /*object id*/
      if(a != -1){
        arg->predictioncounter->data[a] += 1; /* this object was visited */
        /* updating y */
        for(n = 0; n < y_test_predicted->col; n++){
          arg->predicted_y->data[a][n] +=  y_test_predicted->data[k][n];
        }
        k++;
      }
      else{
        continue;
      }
    }

    DelMatrix(&y_test_predicted);
    DelPLSModel(&subm);
    DelMatrix(&x_train);
    DelMatrix(&y_train);
    DelMatrix(&x_test);
    DelMatrix(&y_test);
  }

  DelMatrix(&gid);
  return 0;
}

void *MLRRandomGroupCVModel(void *arg_)
{
  size_t j, k, n, g;
  rgcv_th_arg *arg;
  matrix *gid; /* randomization and storing id for each random group into a matrix */

  /* Matrix for compute the PLS models for groups */
  matrix *x_train;
  matrix *y_train;
  MLRMODEL *subm;

  /*matrix to predict*/
  matrix *x_test;
  matrix *y_test;
  matrix *y_test_predicted;

  arg = (rgcv_th_arg*) arg_;

  /* step 1 generate the random groups */
  initMatrix(&gid);
  random_kfold_group_generator(gid, arg->group, arg->mx->row, &arg->srand_init);

  /*
  puts("Gid Matrix");
  PrintMatrix(gid);
  */

  /* step 2 */
  for(g = 0; g < gid->row; g++){ /*For aeach group */
    /* Estimate how many objects are inside the group "g" to predict and outside the group "g" to build the model  */
    initMatrix(&x_train);
    initMatrix(&y_train);
    initMatrix(&x_test);
    initMatrix(&y_test); /* unused variable here .... */

    kfold_group_train_test_split(arg->mx,
                                 arg->my,
                                 gid,
                                 g,
                                 x_train,
                                 y_train,
                                 x_test,
                                 y_test);

    NewMLRModel(&subm);

    MLR(x_train, y_train, subm, NULL);

    initMatrix(&y_test_predicted);
    MLRPredictY(x_test, NULL, subm, y_test_predicted, NULL, NULL, NULL);

    for(j = 0, k = 0; j < gid->col; j++){
      int a = (int)gid->data[g][j]; /*object id*/
      if(a != -1){
        arg->predictioncounter->data[a] += 1; /* this object was visited */
        /* updating y */
        for(n = 0; n < y_test_predicted->col; n++){
          arg->predicted_y->data[a][n] +=  y_test_predicted->data[k][n];
        }
        k++;
      }
      else{
        continue;
      }
    }

    DelMatrix(&y_test_predicted);
    DelMLRModel(&subm);
    DelMatrix(&x_train);
    DelMatrix(&y_train);
    DelMatrix(&x_test);
    DelMatrix(&y_test);
  }

  DelMatrix(&gid);
  return 0;
}

void *EPLSRandomGroupCVModel(void *arg_)
{
  size_t j, k, n, g;
  rgcv_th_arg *arg;
  matrix *gid; /* randomization and storing id for each random group into a matrix */

  /* Matrix for compute the PLS models for groups */
  matrix *x_train;
  matrix *y_train;
  EPLSMODEL *subm;

  /*matrix to predict*/
  matrix *x_test;
  matrix *y_test;
  matrix *y_test_predicted;

  arg = (rgcv_th_arg*) arg_;

  /* step 1 generate the random groups */
  initMatrix(&gid);
  random_kfold_group_generator(gid,
                               arg->group,
                               arg->mx->row,
                               &arg->srand_init);
  /*
  puts("Gid Matrix");
  PrintMatrix(gid);
  */

  /* step 2 */
  for(g = 0; g < gid->row; g++){ /*For aeach group */
    /* Estimate how many objects are inside the group "g" to predict and outside the group "g" to build the model  */
    initMatrix(&x_train);
    initMatrix(&y_train);
    initMatrix(&x_test);
    initMatrix(&y_test); /* unused variable here .... */
    kfold_group_train_test_split(arg->mx,
                                 arg->my,
                                 gid,
                                 g,
                                 x_train,
                                 y_train,
                                 x_test,
                                 y_test);

    NewEPLSModel(&subm);

    EPLS(x_train,
         y_train,
         arg->nlv,
         arg->xautoscaling,
         arg->yautoscaling,
         subm,
         arg->eparm,
         NULL);

    initMatrix(&y_test_predicted);
    EPLSYPRedictorAllLV(x_test, subm, arg->crule, NULL, &y_test_predicted);

    for(j = 0, k = 0; j < gid->col; j++){
      int a = (int)gid->data[g][j]; /*object id*/
      if(a != -1){
        arg->predictioncounter->data[a] += 1; /* this object was visited */
        /* updating y */
        for(n = 0; n < y_test_predicted->col; n++){
          arg->predicted_y->data[a][n] +=  y_test_predicted->data[k][n];
        }
        k++;
      }
      else{
        continue;
      }
    }

    DelMatrix(&y_test_predicted);
    DelEPLSModel(&subm);
    DelMatrix(&x_train);
    DelMatrix(&y_train);
    DelMatrix(&x_test);
    DelMatrix(&y_test);
  }

  DelMatrix(&gid);
  return 0;
}

void *LDARandomGroupCVModel(void *arg_)
{
  size_t j, k, n, g;
  rgcv_th_arg *arg;
  matrix *gid; /* randomization and storing id for each random group into a matrix */

  /* Matrix for compute the PLS models for groups */
  matrix *x_train;
  matrix *y_train;
  LDAMODEL *subm;

  /*matrix to predict*/
  matrix *x_test;
  matrix *y_test;
  matrix *y_test_predicted;

  matrix *pfeatures;
  matrix *probability;
  matrix *mnpdf;

  arg = (rgcv_th_arg*) arg_;

  /* step 1 generate the random groups */
  initMatrix(&gid);
  random_kfold_group_generator(gid, arg->group, arg->mx->row, &arg->srand_init);

  /*
  puts("Gid Matrix");
  PrintMatrix(gid);
  */

  /* step 2 */
  for(g = 0; g < gid->row; g++){ /*For aeach group */
    /* Estimate how many objects are inside the group "g" to predict and outside the group "g" to build the model  */
    initMatrix(&x_train);
    initMatrix(&y_train);
    initMatrix(&x_test);
    initMatrix(&y_test); /* unused variable here .... */

    kfold_group_train_test_split(arg->mx,
                                 arg->my,
                                 gid,
                                 g,
                                 x_train,
                                 y_train,
                                 x_test,
                                 y_test);

    NewLDAModel(&subm);

    LDA(x_train, y_train, subm);

    initMatrix(&y_test_predicted);
    initMatrix(&pfeatures);
    initMatrix(&probability);
    initMatrix(&mnpdf);

    LDAPrediction(x_test, subm, pfeatures, probability, mnpdf, y_test_predicted);

    for(j = 0, k = 0; j < gid->col; j++){
      int a = (int)gid->data[g][j]; /*object id*/
      if(a != -1){
        arg->predictioncounter->data[a] += 1; /* this object was visited */
        /* updating y */
        for(n = 0; n < y_test_predicted->col; n++){
          arg->predicted_y->data[a][n] +=  y_test_predicted->data[k][n];
        }
        k++;
      }
      else{
        continue;
      }
    }

    DelMatrix(&mnpdf);
    DelMatrix(&probability);
    DelMatrix(&pfeatures);
    DelMatrix(&y_test_predicted);
    DelLDAModel(&subm);
    DelMatrix(&x_train);
    DelMatrix(&y_train);
    DelMatrix(&x_test);
    DelMatrix(&y_test);
  }

  DelMatrix(&gid);
  return 0;
}

void BootstrapRandomGroupsCV(MODELINPUT *input,
                             size_t group,
                             size_t iterations,
                             AlgorithmType algo,
                             matrix *predicted_y,
                             matrix *pred_residuals,
                             size_t nthreads,
                             ssignal *s,
                             int num_arg,
                             ...)
{
  matrix *mx = input->mx;
  matrix *my = input->my;
  size_t i, j;
  size_t nlv;
  size_t xautoscaling;
  size_t yautoscaling;
  ELearningParameters eparm;
  CombinationRule crule = Averaging;

  if(algo == _PLS_ || algo == _PLS_DA_ || algo == _EPLS_ || algo == _EPLS_DA_){
    if(input->nlv > mx->col){
      nlv = mx->col;
    }
    else{
      nlv = input->nlv;
    }

    xautoscaling = input->xautoscaling;
    yautoscaling = input->yautoscaling;
  }
  else{
    nlv = 0;
    xautoscaling = 0;
    yautoscaling = 0;
  }

  if(num_arg > 0){
    va_list valist;
    va_start(valist, num_arg);
    for (i = 0; i < num_arg; i++){
      if(i == 0)
        eparm = va_arg(valist, ELearningParameters);
      else if(i == 1)
        crule = va_arg(valist, CombinationRule);
      else
        continue;
    }
    /* clean memory reserved for valist */
    va_end(valist);
  }

  if(mx->row == my->row && group > 0 && iterations > 0){
    size_t th, iterations_;
    pthread_t *threads;
    uivector *predictcounter;
    matrix *sum_ypredictions;

    size_t scol;
    if(nlv > 0){
      if(nlv > mx->col){
        nlv = mx->col;
      }
      scol = my->col*nlv;  /* each component have my->col ypsilon */
    }
    else{
      nlv = 1;
      scol = my->col;
    }

    NewMatrix(&sum_ypredictions, my->row, scol);
    NewUIVector(&predictcounter, my->row);

    /* each thread have its argument type */
    threads = xmalloc(sizeof(pthread_t)*nthreads);
    rgcv_th_arg *arg;
    arg = xmalloc(sizeof(rgcv_th_arg)*nthreads);

    /*
    iterations_ = 0;
    while(iterations_ <  iterations){*/
    for(iterations_ = 0; iterations_ < iterations; iterations_ += nthreads){
      if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
        break;
      }
      else{
        /* Create independent threads each of which will execute function */
        for(th = 0; th < nthreads; th++){
          arg[th].eparm = eparm;
          arg[th].crule = crule;
          arg[th].mx = mx;
          arg[th].my = my;
          arg[th].group = group;
          arg[th].nlv = nlv;
          arg[th].xautoscaling = xautoscaling;
          arg[th].yautoscaling = yautoscaling;
          arg[th].srand_init = (size_t) group + mx->row + my->col + iterations + th + iterations_;
          NewMatrix(&arg[th].predicted_y, my->row, scol);
          NewUIVector(&arg[th].predictioncounter, my->row);
          if(algo == _PLS_ || algo == _PLS_DA_){
            pthread_create(&threads[th], NULL, PLSRandomGroupCVModel, (void*) &arg[th]);
          }
          else if(algo == _MLR_){
            pthread_create(&threads[th], NULL, MLRRandomGroupCVModel, (void*) &arg[th]);
          }
          else if(algo == _EPLS_ || algo == _EPLS_DA_){
            pthread_create(&threads[th], NULL, EPLSRandomGroupCVModel, (void*) &arg[th]);
          }
          else if(algo == _LDA_){
            pthread_create(&threads[th], NULL, LDARandomGroupCVModel, (void*) &arg[th]);
          }
          else{
            continue;
          }
        }

        /* Wait till threads are complete before main continues. Unless we  */
        /* wait we run the risk of executing an exit which will terminate   */
        /* the process and all threads before the threads have completed.   */
        for(th = 0; th < nthreads; th++){
          pthread_join(threads[th], NULL);
        }

        /* finalize thread outputs and free the memory.....*/
        for(th = 0; th < nthreads; th++){
          for(i = 0; i < arg[th].predicted_y->row; i++){
            for(j = 0; j < arg[th].predicted_y->col; j++){
              sum_ypredictions->data[i][j] += arg[th].predicted_y->data[i][j];
            }

            predictcounter->data[i] += arg[th].predictioncounter->data[i];
          }

          DelUIVector(&arg[th].predictioncounter);
          DelMatrix(&arg[th].predicted_y);
        }
      }
    }

    /*Finalize the output by dividing for the number of times that the object was predicted*/

    if(predicted_y != NULL){
      ResizeMatrix(predicted_y, sum_ypredictions->row, sum_ypredictions->col);
    }

    if(pred_residuals != NULL){
      ResizeMatrix(pred_residuals, my->row, my->col*nlv); /* each component have my->col ypsilon */
    }

    for(i = 0; i < sum_ypredictions->row; i++){
      for(j = 0; j < sum_ypredictions->col; j++){
        sum_ypredictions->data[i][j] /= (double)predictcounter->data[i];

        if(predicted_y != NULL)
          predicted_y->data[i][j] = sum_ypredictions->data[i][j];

        if(pred_residuals != NULL)
          pred_residuals->data[i][j] = sum_ypredictions->data[i][j] - my->data[i][(size_t)floor(j/nlv)];
      }
    }

    DelMatrix(&sum_ypredictions);
    DelUIVector(&predictcounter);
    xfree(threads);
    xfree(arg);
  }
  else{
    char *algo_;
    if(algo == _PLS_)
      algo_ = "PLS";
    else if(algo == _PLS_DA_)
      algo_ = "PLS-DA";
    else if(algo == _MLR_)
      algo_ = "MLR";
    else
      algo_ = "LDA";

    fprintf(stderr, "Error!! Unable to compute Random Group Cross Validation for %s\n", algo_);
  }
}

/* Leave One Out Cross Validation
 *
 * 1) remove one object
 * 2) calculate the model
 * 3) predict the removed object and so the r2 and q2
 */

typedef struct{
  ELearningParameters eparm;
  CombinationRule crule;
  matrix *x_train, *y_train, *x_test, *y_test, *y_test_predicted;
  size_t nlv, xautoscaling, yautoscaling;
} loocv_th_arg;

void *MLRLOOModel_(void *arg_)
{
  loocv_th_arg *arg;
  arg = (loocv_th_arg*) arg_;

  MLRMODEL *subm;
  NewMLRModel(&subm);

  MLR(arg->x_train, arg->y_train, subm, NULL);

  MLRPredictY(arg->x_test, NULL, subm, arg->y_test_predicted, NULL, NULL, NULL);

  DelMLRModel(&subm);
  return 0;
}

void *PLSLOOModel_(void *arg_)
{
  loocv_th_arg *arg;
  arg = (loocv_th_arg*) arg_;

  PLSMODEL *subm;
  NewPLSModel(&subm);
  PLS(arg->x_train, arg->y_train, arg->nlv, arg->xautoscaling, arg->yautoscaling, subm, NULL);
  /* Predict Y for each latent variable */
  PLSYPredictorAllLV(arg->x_test, subm, NULL, arg->y_test_predicted);
  DelPLSModel(&subm);
  return 0;
}

void *EPLSLOOModel_(void *arg_)
{
  loocv_th_arg *arg;
  arg = (loocv_th_arg*) arg_;
  EPLSMODEL *subm;
  NewEPLSModel(&subm);
  EPLS(arg->x_train, arg->y_train, arg->nlv, arg->xautoscaling, arg->yautoscaling, subm, arg->eparm, NULL);
  //Predict Y for each latent variable
  EPLSYPRedictorAllLV(arg->x_test, subm, arg->crule, NULL, &arg->y_test_predicted);
  DelEPLSModel(&subm);
  return 0;
}

void *LDALOOModel_(void *arg_)
{
  loocv_th_arg *arg;
  arg = (loocv_th_arg*) arg_;
  matrix *probability;
  matrix *pfeatures;
  matrix *mnpdf;
  LDAMODEL *subm;
  NewLDAModel(&subm);

  LDA(arg->x_train, arg->y_train, subm);

  initMatrix(&pfeatures);
  initMatrix(&probability);
  initMatrix(&mnpdf);

  LDAPrediction(arg->x_test, subm, pfeatures, probability, mnpdf, arg->y_test_predicted);

  DelMatrix(&mnpdf);
  DelMatrix(&probability);
  DelMatrix(&pfeatures);
  DelLDAModel(&subm);
  return 0;
}

void LeaveOneOut(MODELINPUT *input,
                 AlgorithmType algo,
                 matrix *predicted_y,
                 matrix *pred_residuals,
                 size_t nthreads,
                 ssignal *s,
                 int num_arg,
                 ...)
{
  size_t i, j, k, l, th, model;
  va_list valist;
  pthread_t *threads;
  loocv_th_arg *loo_arg;
  size_t scol;
  matrix *loopredictedy;
  matrix *mx = input->mx;
  matrix *my = input->my;
  size_t nlv;
  size_t xautoscaling;
  size_t yautoscaling;
  ELearningParameters eparm;
  CombinationRule crule = Averaging;

  if(algo == _PLS_ || algo == _PLS_DA_ || algo == _EPLS_ || algo == _EPLS_DA_){
    if(input->nlv > mx->col){
      nlv = mx->col;
    }
    else{
      nlv = input->nlv;
    }

    xautoscaling = input->xautoscaling;
    yautoscaling = input->yautoscaling;
  }
  else{
    nlv = 0;
    xautoscaling = 0;
    yautoscaling = 0;
  }

  if(num_arg > 0){
    va_start(valist, num_arg);
    for(i = 0; i < num_arg; i++){
      if(i == 0)
        eparm = va_arg(valist, ELearningParameters);
      else if(i == 1)
        crule = va_arg(valist, CombinationRule);
      else
        continue;
    }
    /* clean memory reserved for valist */
    va_end(valist);
  }

  if(mx->row == my->row){

    if(nlv > 0){
      if(nlv > mx->col){
        nlv = mx->col;
      }
      scol = my->col*nlv;  /* each component have my->col ypsilon */
    }
    else{
      nlv = 1;
      scol = my->col;
    }


    NewMatrix(&loopredictedy, my->row, scol);

    threads = xmalloc(sizeof(pthread_t)*nthreads);
    loo_arg = xmalloc(sizeof(loocv_th_arg)*nthreads);

    /* initialize threads arguments.. */
    for(th = 0; th < nthreads; th++){
      loo_arg[th].eparm = eparm;
      loo_arg[th].crule = crule;
      loo_arg[th].nlv = nlv;
      loo_arg[th].xautoscaling = xautoscaling;
      loo_arg[th].yautoscaling = yautoscaling;
      NewMatrix(&loo_arg[th].x_train, mx->row-1, mx->col);
      NewMatrix(&loo_arg[th].y_train, my->row-1, my->col);
      NewMatrix(&loo_arg[th].x_test, 1, mx->col);
      NewMatrix(&loo_arg[th].y_test, 1, my->col);
      NewMatrix(&loo_arg[th].y_test_predicted, 1, my->col*nlv);
    }

    for(model = 0; model < mx->row; model += nthreads){ /* we compute mx->row models  */
      if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
        break;
      }
      else{
        /* copy data into subX, subY, predictX and predictY for each thread argument
         * and run the thread
         */
        for(th = 0; th < nthreads; th++){
          if(th+model < mx->row){
            l = 0;
            for(j = 0; j < mx->row; j++){
              if(j != model+th){
                for(k = 0; k  < mx->col; k++){
                  /*setMatrixValue(arg[th].submx, l, k, getMatrixValue(mx, j, k));*/
                  loo_arg[th].x_train->data[l][k] = mx->data[j][k];
                }
                for(k = 0; k < my->col; k++){
                  /*setMatrixValue(arg[th].submy, l, k, getMatrixValue(my, j, k));*/
                  loo_arg[th].y_train->data[l][k] = my->data[j][k];
                }
                l++;
              }
              else{
                for(k = 0; k < mx->col; k++){
                  loo_arg[th].x_test->data[0][k] = mx->data[j][k];
                }
                for(k = 0; k < my->col; k++){
                  loo_arg[th].y_test->data[0][k] = my->data[j][k];
                }
              }
            }

            if(algo == _PLS_ || algo == _PLS_DA_)
              pthread_create(&threads[th], NULL, PLSLOOModel_, (void*) &loo_arg[th]);
            else if(algo == _MLR_)
              pthread_create(&threads[th], NULL, MLRLOOModel_, (void*) &loo_arg[th]);
            else if(algo == _EPLS_ || algo == _EPLS_DA_)
              pthread_create(&threads[th], NULL, EPLSLOOModel_, (void*) &loo_arg[th]);
            else if(algo == _LDA_)
                pthread_create(&threads[th], NULL, LDALOOModel_, (void*) &loo_arg[th]);
            else
              continue;
          }
          else{
            continue;
          }
        }

        /* Wait till threads are complete before main continues. Unless we  */
        /* wait we run the risk of executing an exit which will terminate   */
        /* the process and all threads before the threads have completed.   */
        for(th = 0; th < nthreads; th++){
          if(th+model < mx->row){
            pthread_join(threads[th], NULL);
          }
          else{
            continue;
          }
        }

        /*Collapse the threads output*/
        for(th = 0; th < nthreads; th++){
          if(th+model < mx->row){
            for(j = 0; j < loo_arg[th].y_test_predicted->col; j++){
              loopredictedy->data[model+th][j] = loo_arg[th].y_test_predicted->data[0][j];
            }
          }
        }
      }
    }

    /*Delete thread arguments*/

    for(th = 0; th < nthreads; th++){
      DelMatrix(&loo_arg[th].x_train);
      DelMatrix(&loo_arg[th].y_train);
      DelMatrix(&loo_arg[th].x_test);
      DelMatrix(&loo_arg[th].y_test);
      DelMatrix(&loo_arg[th].y_test_predicted);
    }

    /*Finalize the output by dividing for the number of models*/
    if(predicted_y != NULL){
      ResizeMatrix(predicted_y, loopredictedy->row, loopredictedy->col);
      MatrixCopy(loopredictedy, &predicted_y);
    }

    if(pred_residuals != NULL){
      ResizeMatrix(pred_residuals, my->row, my->col*nlv); /* each component have my->col ypsilon */
      for(i = 0; i < loopredictedy->row; i++){
        for(j = 0; j < loopredictedy->col; j++){
          pred_residuals->data[i][j] = loopredictedy->data[i][j] - my->data[i][(size_t)floor(j/nlv)];
        }
      }
    }

    DelMatrix(&loopredictedy);
    xfree(loo_arg);
    xfree(threads);
  }
  else{
    fprintf(stderr, "Error!! Unable to compute PLS Leave One Out Validation!!\n");
  }
}


void KFoldCV(MODELINPUT *input,
             uivector *groups,
             AlgorithmType algo,
             matrix *predicted_y,
             matrix *pred_residuals,
             size_t nthreads,
             ssignal *s,
             int arg,
             ...)
{
  matrix *mx = input->mx;
  matrix *my = input->my;
  size_t i, j;
  size_t nlv;
  size_t xautoscaling;
  size_t yautoscaling;
  ELearningParameters eparm;
  CombinationRule crule = Averaging;
  matrix *y_predicted;


  if(algo == _PLS_ || algo == _PLS_DA_ || algo == _EPLS_ || algo == _EPLS_DA_){
    if(input->nlv > mx->col){
      nlv = mx->col;
    }
    else{
      nlv = input->nlv;
    }

    xautoscaling = input->xautoscaling;
    yautoscaling = input->yautoscaling;
  }
  else{
    nlv = 0;
    xautoscaling = 0;
    yautoscaling = 0;
  }

  if(arg > 0){
    va_list valist;
    va_start(valist, arg);
    for(i = 0; i < arg; i++){
      if(i == 0)
        eparm = va_arg(valist, ELearningParameters);
      else if(i == 1)
        crule = va_arg(valist, CombinationRule);
      else
        continue;
    }
    /* clean memory reserved for valist */
    va_end(valist);
  }

  if(mx->row == my->row && groups->size > 0){
    size_t th, it;
    pthread_t *threads;
    loocv_th_arg *kcv_arg;
    matrix *gid;

    size_t scol;
    if(nlv > 0){
      if(nlv > mx->col){
        nlv = mx->col;
      }
      scol = my->col*nlv;  /* each component have my->col ypsilon */
    }
    else{
      nlv = 1;
      scol = my->col;
    }

    /* Create groups matrix */
    size_t gmax = groups->data[0];
    for(i = 1; i < groups->size; i++){
      if(groups->data[i] > gmax)
        gmax = groups->data[i];
      else
        continue;
    }
    gmax++;

    size_t objgmax = 0, objgmax_tmp;
    for(i = 0; i < gmax; i++){
      objgmax_tmp = 0;
      for(j = 0; j < groups->size; j++){
        if(groups->data[j] == i){
          objgmax_tmp++;
        }
        else{
          continue;
        }
      }
      if(objgmax_tmp > objgmax){
        objgmax = objgmax_tmp;
      }
      else{
        continue;
      }
    }

    NewMatrix(&gid, gmax, objgmax);
    MatrixSet(gid, -1);
    uivector *indx;
    NewUIVector(&indx, gmax);
    for(j = 0; j < groups->size; j++){
      size_t g = groups->data[j];
      gid->data[g][indx->data[g]] = j;
      indx->data[g] += 1;
    }
    DelUIVector(&indx);


    /* each thread have its argument type */
    threads = xmalloc(sizeof(pthread_t)*nthreads);
    kcv_arg = xmalloc(sizeof(loocv_th_arg)*nthreads);
    NewMatrix(&y_predicted, my->row, scol);

    for(it = 0; it < gmax; it += nthreads){
      if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
        break;
      }
      else{
        /* Create independent threads each of which will execute function */
        for(th = 0; th < nthreads; th++){
          if(it+th >= gmax){
            break;
          }
          else{
            kcv_arg[th].eparm = eparm;
            kcv_arg[th].crule = crule;
            kcv_arg[th].nlv = nlv;
            kcv_arg[th].xautoscaling = xautoscaling;
            kcv_arg[th].yautoscaling = yautoscaling;
            initMatrix(&kcv_arg[th].x_train);
            initMatrix(&kcv_arg[th].y_train);
            initMatrix(&kcv_arg[th].x_test);
            initMatrix(&kcv_arg[th].y_test);
            initMatrix(&kcv_arg[th].y_test_predicted);
            kfold_group_train_test_split(mx,
                                         my,
                                         gid,
                                         it+th,
                                         kcv_arg[th].x_train,
                                         kcv_arg[th].y_train,
                                         kcv_arg[th].x_test,
                                         kcv_arg[th].y_test);

            if(algo == _PLS_ || algo == _PLS_DA_){
              pthread_create(&threads[th], NULL, PLSLOOModel_, (void*) &kcv_arg[th]);
            }
            else if(algo == _MLR_){
              pthread_create(&threads[th], NULL, MLRLOOModel_, (void*) &kcv_arg[th]);
            }
            else if(algo == _EPLS_ || algo == _EPLS_DA_){
              pthread_create(&threads[th], NULL, EPLSLOOModel_, (void*) &kcv_arg[th]);
            }
            else{
              continue;
            }
          }
        }

        /* Wait till threads are complete before main continues. Unless we  */
        /* wait we run the risk of executing an exit which will terminate   */
        /* the process and all threads before the threads have completed.   */
        for(th = 0; th < nthreads; th++){
          if(it+th >= gmax){
            break;
          }
          else{
            pthread_join(threads[th], NULL);
          }
        }


        /* finalize thread outputs and free the memory.....*/
        for(th = 0; th < nthreads; th++){
          if(it+th >= gmax){
            break;
          }
          else{
            for(i = 0; i < kcv_arg[th].y_test_predicted->row; i++){
              int id = (int)gid->data[th+it][i];
              for(j = 0; j < kcv_arg[th].y_test_predicted->col; j++){
                y_predicted->data[id][j] = kcv_arg[th].y_test_predicted->data[i][j];
                //setMatrixValue(y_predicted, id, j, kcv_arg[th].y_test_predicted->data[i][j]);
              }
            }

            DelMatrix(&kcv_arg[th].x_train);
            DelMatrix(&kcv_arg[th].y_train);
            DelMatrix(&kcv_arg[th].x_test);
            DelMatrix(&kcv_arg[th].y_test);
            DelMatrix(&kcv_arg[th].y_test_predicted);
          }
        }
      }
    }

    /*Finalize the output by dividing for the number of times that the object was predicted*/

    if(predicted_y != NULL){
      MatrixCopy(y_predicted, &predicted_y);
    }

    if(pred_residuals != NULL){
      ResizeMatrix(pred_residuals, my->row, my->col*nlv); /* each component have my->col ypsilon */

      for(i = 0; i < y_predicted->row; i++){
        for(j = 0; j < y_predicted->col; j++){
          pred_residuals->data[i][j] = y_predicted->data[i][j] - my->data[i][(size_t)floor(j/nlv)];
        }
      }
    }


    DelMatrix(&y_predicted);
    DelMatrix(&gid);
    xfree(threads);
    xfree(kcv_arg);
  }
  else{
    char *algo_;
    if(algo == _PLS_)
      algo_ = "PLS";
    else if(algo == _PLS_DA_)
      algo_ = "PLS-DA";
    else if(algo == _MLR_)
      algo_ = "MLR";
    else
      algo_ = "LDA";

    fprintf(stderr, "Error!! Unable to compute Random Group Cross Validation for %s\n", algo_);
  }
}

/* Get the best r2 and q2 for PLS Model */
void PLSRegressionYScramblingPipeline(matrix *mx,
                                      matrix *my,
                                      size_t xautoscaling,
                                      size_t yautoscaling,
                                      size_t nlv,
                                      ValidationArg varg,
                                      size_t nthreads,
                                      dvector *r2,
                                      dvector *q2)
{
  size_t i;
  PLSMODEL *tmpmod;
  matrix *py;
  matrix *pres;
  matrix *yrec;
  matrix *tmpq2;

  MODELINPUT minpt;
  minpt.mx = mx;
  minpt.my = my;
  minpt.nlv = nlv;
  minpt.xautoscaling = xautoscaling;
  minpt.yautoscaling = yautoscaling;

  /*compude q2y*/
  initMatrix(&py);
  initMatrix(&pres);
  if(varg.vtype == LOO){
    LeaveOneOut(&minpt, _PLS_, py, pres, nthreads, NULL, 0);
  }
  else{
    BootstrapRandomGroupsCV(&minpt, 3, 100, _PLS_, py, pres, 4, NULL, 0);
  }

  initMatrix(&tmpq2);
  PLSRegressionStatistics(my, py, tmpq2, NULL, NULL);

  /*Computer r2y*/
  NewPLSModel(&tmpmod);
  PLS(mx, my, nlv, xautoscaling, yautoscaling, tmpmod, NULL);

  initMatrix(&yrec);
  PLSYPredictorAllLV(mx, tmpmod, NULL, yrec);
  PLSRegressionStatistics(my, yrec, tmpmod->r2y_recalculated, NULL, NULL);

  /* Calculate y real vs yscrambled and add other r2 q2 */
  size_t r2cutoff = GetLVCCutoff(tmpmod->r2y_recalculated);
  size_t q2cutoff = GetLVCCutoff(tmpq2);

  for(i = 0; i < my->col; i++){
    r2->data[i] = tmpmod->r2y_recalculated->data[r2cutoff][i];
    q2->data[i] = tmpq2->data[q2cutoff][i];
  }

  DelMatrix(&yrec);
  DelPLSModel(&tmpmod);
  DelMatrix(&tmpq2);
  DelMatrix(&py);
  DelMatrix(&pres);
}

void PLSDiscriminantAnalysisYScramblingPipeline(matrix *mx,
                                                matrix *my,
                                                size_t xautoscaling,
                                                size_t yautoscaling,
                                                size_t nlv,
                                                ValidationArg varg,
                                                size_t nthreads,
                                                dvector *auc_recalc_,
                                                dvector *auc_validation_)
{
  size_t i;
  PLSMODEL *tmpmod;
  matrix *py;
  matrix *pres;
  matrix *yrec;
  matrix *auc_recalc;
  matrix *auc_validation;


  MODELINPUT minpt;
  minpt.mx = mx;
  minpt.my = my;
  minpt.nlv = nlv;
  minpt.xautoscaling = xautoscaling;
  minpt.yautoscaling = yautoscaling;

  /*compude q2y*/
  initMatrix(&py);
  initMatrix(&pres);
  if(varg.vtype == LOO){
    LeaveOneOut(&minpt, _PLS_, py, pres, nthreads, NULL, 0);
  }
  else{
    BootstrapRandomGroupsCV(&minpt, 3, 100, _PLS_, py, pres, 4, NULL, 0);
  }


  /*Computer the model in recalculation */
  NewPLSModel(&tmpmod);
  PLS(mx, my, nlv, xautoscaling, yautoscaling, tmpmod, NULL);
  initMatrix(&yrec);
  initMatrix(&auc_recalc);
  PLSYPredictorAllLV(mx, tmpmod, NULL, yrec);
  PLSDiscriminantAnalysisStatistics(my, yrec, NULL, auc_recalc, NULL, NULL);

  /*compute the model invalidation */
  initMatrix(&auc_validation);
  PLSDiscriminantAnalysisStatistics(my, py, NULL, auc_validation, NULL, NULL);

  /* Calculate y real vs yscrambled and add other r2 q2 */
  size_t auc_cutoff = GetLVCCutoff(auc_validation);

  for(i = 0; i < my->col; i++){
    auc_recalc_->data[i] = auc_recalc->data[auc_cutoff][i];
    auc_validation_->data[i] = auc_validation->data[auc_cutoff][i];
  }

  DelMatrix(&yrec);
  DelPLSModel(&tmpmod);
  DelMatrix(&auc_validation);
  DelMatrix(&auc_recalc);
  DelMatrix(&py);
  DelMatrix(&pres);
}

/* Get the best r2 and q2 for PLS Model */
void MLRYScramblingPipeline(matrix *mx,
                            matrix *my,
                            ValidationArg varg,
                            size_t nthreads,
                            dvector *r2,
                            dvector *q2)
{
  size_t i;
  MLRMODEL *tmpmod;
  matrix *py;
  matrix *pres;
  matrix *yrec;
  dvector *tmpq2;

  MODELINPUT minpt;
  minpt.mx = mx;
  minpt.my = my;

  /*compude q2y*/
  initMatrix(&py);
  initMatrix(&pres);
  if(varg.vtype == LOO){
    LeaveOneOut(&minpt, _MLR_, py, pres, nthreads, NULL, 0);
  }
  else{
    BootstrapRandomGroupsCV(&minpt, 3, 100, _MLR_, py, pres, 4, NULL, 0);
  }

  initDVector(&tmpq2);
  MLRRegressionStatistics(my, py, tmpq2, NULL, NULL);

  /*Computer r2y*/
  NewMLRModel(&tmpmod);
  MLR(mx, my, tmpmod, NULL);

  initMatrix(&yrec);
  MLRPredictY(mx, my, tmpmod, yrec,  NULL, tmpmod->r2y_model, NULL);

  for(i = 0; i < my->col; i++){
    r2->data[i] = tmpmod->r2y_model->data[i];
    q2->data[i] = tmpq2->data[i];
  }

  DelMatrix(&yrec);
  DelMLRModel(&tmpmod);
  DelDVector(&tmpq2);
  DelMatrix(&py);
  DelMatrix(&pres);
}


void LDAYScramblingPipeline(matrix *mx,
                            matrix *my,
                            ValidationArg varg,
                            size_t nthreads,
                            dvector *auc_recalc_,
                            dvector *auc_validation_)
{
  size_t i;
  LDAMODEL *tmpmod;
  matrix *yrec;
  matrix *py;
  matrix *probability;
  matrix *pfeatures;
  matrix *mnpdf;
  dvector *roc_auc;
  dvector *pr_auc;

  MODELINPUT minpt;
  minpt.mx = mx;
  minpt.my = my;

  /*compude q2y*/
  initMatrix(&py);
  if(varg.vtype == LOO){
    LeaveOneOut(&minpt, _LDA_, py, NULL, nthreads, NULL, 0);
  }
  else{
    BootstrapRandomGroupsCV(&minpt, 3, 100, _LDA_, py, NULL, 4, NULL, 0);
  }

  /*Computer the model in recalculation */
  NewLDAModel(&tmpmod);
  LDA(mx, my, tmpmod);
  initMatrix(&yrec);

  initMatrix(&pfeatures);
  initMatrix(&probability);
  initMatrix(&mnpdf);

  LDAPrediction(mx, tmpmod, pfeatures, probability, mnpdf, yrec);

  DelMatrix(&pfeatures);
  DelMatrix(&probability);
  DelMatrix(&mnpdf);

  initDVector(&roc_auc);
  initDVector(&pr_auc);
  LDAMulticlassStatistics(my,
                          yrec,
                          NULL,
                          roc_auc,
                          NULL,
                          pr_auc);

  for(i = 0; i < roc_auc->size; i++){
    auc_recalc_->data[i] = roc_auc->data[i];
  }

  DelDVector(&roc_auc);
  DelDVector(&pr_auc);
  initDVector(&roc_auc);
  initDVector(&pr_auc);
  LDAMulticlassStatistics(my,
                          py,
                          NULL,
                          roc_auc,
                          NULL,
                          pr_auc);
  for(i = 0; i < roc_auc->size; i++){
    auc_validation_->data[i] = roc_auc->data[i];
  }
  DelDVector(&roc_auc);
  DelDVector(&pr_auc);

  DelMatrix(&yrec);
  DelMatrix(&py);
  DelLDAModel(&tmpmod);
}

void PermutedObservetCorrelation(matrix *my_true,
                                 matrix *my_perm,
                                 dvector *ccoef,
                                 AlgorithmType algo)
{
  size_t i, j;
  dvector *yt;
  dvector *yp;
  double score;

  for(j = 0; j < my_true->col; j++){
    initDVector(&yt);
    initDVector(&yp);

    for(i = 0; i < my_true->row; i++){
      DVectorAppend(yt, my_true->data[i][j]);
      DVectorAppend(yp, my_perm->data[i][j]);
    }

    if(algo == _PLS_ ||
       algo == _MLR_ ||
       algo == _EPLS_){
       score = R2(yt, yp);
    }
    else{
      matrix *roc;
      initMatrix(&roc);
      ROC(yt, yp, roc, &score);
      DelMatrix(&roc);
    }
    ccoef->data[j] = score;
    DelDVector(&yt);
    DelDVector(&yp);
  }
}

void YScrambling(MODELINPUT *input,
                 AlgorithmType algo,
                 ValidationArg varg,
                 size_t iterations,
                 matrix *ccoeff_yscrambling,
                 size_t nthreads,
                 ssignal *s)
{
  size_t i, j, k, it, outcols;
  double ytmp;

  matrix *mx, *my;
  matrix *randomY;
  dvector *corrpermobs, *coeff_int, *coeff_valid;
  mx = input->mx;
  my = input->my;

  srand_(mx->row+mx->col+my->col+iterations);

  initMatrix(&randomY);
  MatrixCopy(my, &randomY);

  /*Create a the r2q2scrambling matrix */
  /*
   * ccoeff_yscrambling columns:
   *  0: correlation between permuted and observed y
   *  1: r2 for the model
   *  2: q2 for the model
   */


  if(algo == _PLS_ || algo == _PLS_DA_ || algo == _MLR_){
    NewDVector(&corrpermobs, my->col);
    NewDVector(&coeff_int, my->col);
    NewDVector(&coeff_valid, my->col);
    outcols = my->col;
  }
  else{
    int n_classes = getNClasses(input->my);
    if(n_classes == 2){
      n_classes = 1;
    }
    NewDVector(&corrpermobs, n_classes);
    NewDVector(&coeff_int, n_classes);
    NewDVector(&coeff_valid, n_classes);
    outcols = n_classes;
  }


  ResizeMatrix(ccoeff_yscrambling, 1+iterations, 3*outcols);

  PermutedObservetCorrelation(my, my, corrpermobs, algo);
  /*The first row is the model not scrambled...*/
  if(algo == _PLS_){
    PLSRegressionYScramblingPipeline(mx,
                                     my,
                                     input->xautoscaling,
                                     input->yautoscaling,
                                     input->nlv,
                                     varg,
                                     nthreads,
                                     coeff_int,
                                     coeff_valid);
  }
  else if(algo == _PLS_DA_){
    PLSDiscriminantAnalysisYScramblingPipeline(mx,
                                               my,
                                               input->xautoscaling,
                                               input->yautoscaling,
                                               input->nlv,
                                               varg,
                                               nthreads,
                                               coeff_int,
                                               coeff_valid);
  }
  else if(algo == _MLR_){
    MLRYScramblingPipeline(mx,
                           my,
                           varg,
                           nthreads,
                           coeff_int,
                           coeff_valid);
  }
  else if(algo == _LDA_){
    LDAYScramblingPipeline(mx,
                           my,
                           varg,
                           nthreads,
                           coeff_int,
                           coeff_valid);
  }

  for(j = 0, k = 0; j < outcols; j++){
    ccoeff_yscrambling->data[0][k] = corrpermobs->data[j];
    k++;
    ccoeff_yscrambling->data[0][k] = coeff_int->data[j];
    k++;
    ccoeff_yscrambling->data[0][k] = coeff_valid->data[j];
    k++;
  }

  for(it = 0; it < iterations; it++){
    /* Sattolo's algorithm to shuffle the y*/
    i = randomY->row;
    while(i > 1){
      i = i - 1;
      j = randInt(0, i);
      for(k = 0; k < randomY->col; k++){
        ytmp = randomY->data[j][k];
        randomY->data[j][k] = randomY->data[i][k];
        randomY->data[i][k] = ytmp;
      }
    }

    PermutedObservetCorrelation(my, randomY, corrpermobs, algo);
    /*The first row is the model not scrambled...*/
    if(algo == _PLS_){
      PLSRegressionYScramblingPipeline(mx,
                                       randomY,
                                       input->xautoscaling,
                                       input->yautoscaling,
                                       input->nlv,
                                       varg,
                                       nthreads,
                                       coeff_int,
                                       coeff_valid);
    }
    else if(algo == _PLS_DA_){
      PLSDiscriminantAnalysisYScramblingPipeline(mx,
                                                 randomY,
                                                 input->xautoscaling,
                                                 input->yautoscaling,
                                                 input->nlv,
                                                 varg,
                                                 nthreads,
                                                 coeff_int,
                                                 coeff_valid);
    }
    else if(algo == _MLR_){
      MLRYScramblingPipeline(mx,
                             randomY,
                             varg,
                             nthreads,
                             coeff_int,
                             coeff_valid);
    }
    else if(algo == _LDA_){
      LDAYScramblingPipeline(mx,
                             randomY,
                             varg,
                             nthreads,
                             coeff_int,
                             coeff_valid);
    }


    for(j = 0, k = 0; j < outcols; j++){
      ccoeff_yscrambling->data[it+1][k] = corrpermobs->data[j];
      k++;
      ccoeff_yscrambling->data[it+1][k] = coeff_int->data[j];
      k++;
      ccoeff_yscrambling->data[it+1][k] = coeff_valid->data[j];
      k++;
    }

  }

  PrintDVector(corrpermobs);
  DelDVector(&coeff_valid);
  DelDVector(&coeff_int);
  DelDVector(&corrpermobs);
  DelMatrix(&randomY);
}
