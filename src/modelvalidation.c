#include "modelvalidation.h"
#include "memwrapper.h"
#include "mlr.h"
#include "pls.h"
#include "pca.h" /*Using: MatrixAutoScaling(); and calcVarExpressed(); */
#include "epls.h"
#include "numeric.h" /* Using:  if(FLOAT_EQ(NumOne, NumTwo));*/
#include "metricspace.h"
#include <math.h>
#include <pthread.h>
#include <stdarg.h>

void random_kfold_group_generator(matrix **gid, size_t ngroups, size_t nobj, unsigned int *srand_init)
{
  ResizeMatrix(gid, ngroups, (size_t)ceil(nobj/(double)ngroups));
  MatrixSet((*gid), -1);
  size_t i, j, k = 0, n;
  for(i = 0; i <  (*gid)->row; i++){
    for(j = 0; j <  (*gid)->col; j++){
      do{
        n = (size_t)myrand_r(srand_init) % nobj;
      } while(ValInMatrix((*gid), n) == 1 && k < nobj);
      if(k < nobj){
        (*gid)->data[i][j] = n;
        k++;
      }
      else
        continue;
    }
  }
}

void kfold_group_train_test_split(matrix *x, matrix *y, matrix *gid, size_t group_id, matrix **x_train, matrix **y_train, matrix **x_test, matrix **y_test)
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
        size_t a =  (size_t)gid->data[i][j]; /* get the row index */
        if(a != -1){
          for(n = 0; n < x->col; n++){
            (*x_train)->data[k][n] = x->data[a][n];
          }
          for(n = 0; n < y->col; n++){
            (*y_train)->data[k][n] = y->data[a][n];
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
        size_t a = (size_t)gid->data[i][j];
        if(a != -1){
          for(n = 0; n < x->col; n++){
            (*x_test)->data[l][n] = x->data[a][n];
          }
          for(n = 0; n < y->col; n++){
            (*y_test)->data[l][n] = y->data[a][n];
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

void train_test_split(matrix *x, matrix *y, double testsize, matrix **x_train, matrix **y_train, matrix **x_test, matrix **y_test, uivector **testids, unsigned int *srand_init)
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
  for(i = 0; i < testsize_; i++){
    do{
      n = (size_t)myrand_r(srand_init) % x->row;
    } while(UIVectorHasValue((*testids), n) == 0);
    UIVectorAppend(testids, n);

    for(j = 0; j < x->col; j++){
      (*x_test)->data[i][j] = x->data[n][j];
    }
    for(j = 0; j < y->col; j++){
      (*y_test)->data[i][j] = y->data[n][j];
    }
  }
  /* Fill training set */
  for(i = 0, n = 0; i < x->row; i++){
    if(UIVectorHasValue((*testids), i) == 1){
      for(j = 0; j < x->col; j++){
        (*x_train)->data[n][j] = x->data[i][j];
      }
      for(j = 0; j < y->col; j++){
        (*y_train)->data[n][j] = y->data[i][j];
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
  random_kfold_group_generator(&gid, arg->group, arg->mx->row, &arg->srand_init);

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

    kfold_group_train_test_split(arg->mx, arg->my, gid, g, &x_train,&y_train,&x_test, &y_test);

    NewPLSModel(&subm);

    PLS(x_train, y_train, arg->nlv, arg->xautoscaling, arg->yautoscaling, subm, NULL);
    initMatrix(&y_test_predicted);
    PLSYPredictorAllLV(x_test, subm, NULL, &y_test_predicted);

    for(j = 0, k = 0; j < gid->col; j++){
      size_t a = (size_t)gid->data[g][j]; /*object id*/
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
  random_kfold_group_generator(&gid, arg->group, arg->mx->row, &arg->srand_init);

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

    kfold_group_train_test_split(arg->mx, arg->my, gid, g, &x_train,&y_train,&x_test, &y_test);

    NewMLRModel(&subm);

    MLR(x_train, y_train, subm, NULL);

    initMatrix(&y_test_predicted);
    MLRPredictY(x_test, NULL, subm, &y_test_predicted, NULL, NULL, NULL);

    for(j = 0, k = 0; j < gid->col; j++){
      size_t a = (size_t)gid->data[g][j]; /*object id*/
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
  random_kfold_group_generator(&gid, arg->group, arg->mx->row, &arg->srand_init);
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
    kfold_group_train_test_split(arg->mx, arg->my, gid, g, &x_train,&y_train,&x_test, &y_test);

    NewEPLSModel(&subm);

    EPLS(x_train, y_train, arg->nlv, arg->xautoscaling, arg->yautoscaling, subm, arg->eparm, NULL);

    initMatrix(&y_test_predicted);
    EPLSYPRedictorAllLV(x_test, subm, arg->crule, NULL, &y_test_predicted);

    for(j = 0, k = 0; j < gid->col; j++){
      size_t a = (size_t)gid->data[g][j]; /*object id*/
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

void BootstrapRandomGroupsCV(MODELINPUT *input, size_t group, size_t iterations, AlgorithmType algo, matrix **predicted_y, matrix **pred_residuals, size_t nthreads, ssignal *s, int arg, ...)
{
  matrix *mx = (*input->mx);
  matrix *my = (*input->my);
  size_t nlv;
  size_t xautoscaling = input->xautoscaling;
  size_t yautoscaling = input->yautoscaling;
  ELearningParameters eparm;
  CombinationRule crule = Averaging;

  if(input->nlv > mx->col){
    nlv = mx->col;
  }
  else{
    nlv = input->nlv;
  }

  if(arg > 0){
    va_list valist;
    va_start(valist, arg);
    for (size_t i = 0; i < arg; i++){
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
    size_t th, iterations_, i, j;

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
          arg[th].srand_init = (unsigned int) group + mx->row + my->col + iterations + th + iterations_;
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
          (*predicted_y)->data[i][j] = sum_ypredictions->data[i][j];

        if(pred_residuals != NULL)
          (*pred_residuals)->data[i][j] = sum_ypredictions->data[i][j] - my->data[i][(size_t)floor(j/nlv)];
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

  MLRPredictY(arg->x_test, NULL, subm, &arg->y_test_predicted, NULL, NULL, NULL);

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
  PLSYPredictorAllLV(arg->x_test, subm, NULL, &arg->y_test_predicted);
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

void LeaveOneOut(MODELINPUT *input, AlgorithmType algo, matrix** predicted_y, matrix **pred_residuals, size_t nthreads, ssignal *s, int arg_, ...)
{
  size_t i, j, k, l, th, model;
  va_list valist;
  pthread_t *threads;
  loocv_th_arg *loo_arg;
  size_t scol;
  matrix *loopredictedy;
  matrix *mx = (*input->mx);
  matrix *my = (*input->my);
  size_t nlv = input->nlv;
  size_t xautoscaling = input->xautoscaling;
  size_t yautoscaling = input->yautoscaling;
  ELearningParameters eparm;
  CombinationRule crule = Averaging;

  if(input->nlv > mx->col){
    nlv = mx->col;
  }
  else{
    nlv = input->nlv;
  }

  if(arg_ > 0){
    va_start(valist, arg_);
    for(i = 0; i < arg_; i++){
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
      MatrixCopy(loopredictedy, predicted_y);
    }

    if(pred_residuals != NULL){
      ResizeMatrix(pred_residuals, my->row, my->col*nlv); /* each component have my->col ypsilon */
      for(i = 0; i < loopredictedy->row; i++){
        for(j = 0; j < loopredictedy->col; j++){
          (*pred_residuals)->data[i][j] = loopredictedy->data[i][j] - my->data[i][(size_t)floor(j/nlv)];
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

/* Get the best r2 and q2 for PLS Model */
void PLSRegressionYScramblingPipeline(matrix *mx, matrix *my, size_t xautoscaling, size_t yautoscaling, size_t nlv, ValidationArg varg, size_t nthreads, dvector **r2, dvector **q2)
{
  size_t i;
  PLSMODEL *tmpmod;
  matrix *py;
  matrix *pres;
  matrix *yrec;
  matrix *tmpq2;

  MODELINPUT minpt;
  minpt.mx = &mx;
  minpt.my = &my;
  minpt.nlv = nlv;
  minpt.xautoscaling = xautoscaling;
  minpt.yautoscaling = yautoscaling;

  /*compude q2y*/
  initMatrix(&py);
  initMatrix(&pres);
  if(varg.vtype == LOO){
    LeaveOneOut(&minpt, _PLS_, &py, &pres, nthreads, NULL, 0);
  }
  else{
    BootstrapRandomGroupsCV(&minpt, 3, 100, _PLS_, &py, &pres, 4, NULL, 0);
  }

  initMatrix(&tmpq2);
  PLSRegressionStatistics(my, py, &tmpq2, NULL, NULL);

  /*Computer r2y*/
  NewPLSModel(&tmpmod);
  PLS(mx, my, nlv, xautoscaling, yautoscaling, tmpmod, NULL);

  initMatrix(&yrec);
  PLSYPredictorAllLV(mx, tmpmod, NULL, &yrec);
  PLSRegressionStatistics(my, yrec, &(tmpmod->r2y_model), NULL, NULL);

  /* Calculate y real vs yscrambled and add other r2 q2 */
  size_t r2cutoff = GetLVCCutoff(tmpmod->r2y_model);
  size_t q2cutoff = GetLVCCutoff(tmpq2);

  for(i = 0; i < my->col; i++){
    (*r2)->data[i] = tmpmod->r2y_model->data[r2cutoff][i];
    (*q2)->data[i] = tmpq2->data[q2cutoff][i];
  }

  DelMatrix(&yrec);
  DelPLSModel(&tmpmod);
  DelMatrix(&tmpq2);
  DelMatrix(&py);
  DelMatrix(&pres);
}

void PLSDiscriminantAnalysisYScramblingPipeline(matrix *mx, matrix *my, size_t xautoscaling, size_t yautoscaling, size_t nlv, ValidationArg varg, size_t nthreads, dvector **auc_recalc_, dvector **auc_validation_)
{
  size_t i;
  PLSMODEL *tmpmod;
  matrix *py;
  matrix *pres;
  matrix *yrec;
  matrix *auc_recalc;
  matrix *auc_validation;


  MODELINPUT minpt;
  minpt.mx = &mx;
  minpt.my = &my;
  minpt.nlv = nlv;
  minpt.xautoscaling = xautoscaling;
  minpt.yautoscaling = yautoscaling;

  /*compude q2y*/
  initMatrix(&py);
  initMatrix(&pres);
  if(varg.vtype == LOO){
    LeaveOneOut(&minpt, _PLS_, &py, &pres, nthreads, NULL, 0);
  }
  else{
    BootstrapRandomGroupsCV(&minpt, 3, 100, _PLS_, &py, &pres, 4, NULL, 0);
  }


  /*Computer the model in recalculation */
  NewPLSModel(&tmpmod);
  PLS(mx, my, nlv, xautoscaling, yautoscaling, tmpmod, NULL);
  initMatrix(&yrec);
  PLSYPredictorAllLV(mx, tmpmod, NULL, &yrec);
  PLSDiscriminantAnalysisStatistics(my, yrec, NULL, &auc_recalc, NULL, NULL);

  /*compute the model invalidation */
  initMatrix(&auc_validation);
  PLSDiscriminantAnalysisStatistics(my, py, NULL, &auc_validation, NULL, NULL);

  /* Calculate y real vs yscrambled and add other r2 q2 */
  size_t auc_cutoff = GetLVCCutoff(auc_validation);

  for(i = 0; i < my->col; i++){
    (*auc_recalc_)->data[i] = auc_recalc->data[auc_cutoff][i];
    (*auc_validation_)->data[i] = auc_validation->data[auc_cutoff][i];
  }

  DelMatrix(&yrec);
  DelPLSModel(&tmpmod);
  DelMatrix(&auc_validation);
  DelMatrix(&auc_recalc);
  DelMatrix(&py);
  DelMatrix(&pres);
}

void PermutedObservetCorrelation(matrix *my_true, matrix *my_perm, dvector *yaverage, dvector **ccoef)
{
  size_t i, j;
  for(j = 0; j < my_true->col; j++){
    double rss = 0.f, tss = 0.f;
    for(i = 0; i < my_true->row; i++){
      rss += square(my_true->data[i][j] - my_perm->data[i][j]);
      tss += square(my_true->data[i][j] - yaverage->data[j]);
    }
    (*ccoef)->data[j] = 1 - (rss/tss);
  }
}

void YScrambling(MODELINPUT *input, AlgorithmType algo, ValidationArg varg, size_t iterations,
                  matrix **ccoeff_yscrambling, size_t nthreads, ssignal *s)
{
  size_t i, j, k, it;
  double ytmp;

  matrix *mx, *my;
  matrix *randomY;
  dvector *yaverage;
  dvector *corrpermobs, *coeff_int, *coeff_valid;
  mx = (*input->mx);
  my = (*input->my);
  unsigned int srand_init = 1;

  initDVector(&yaverage);
  MatrixColAverage(my, &yaverage);

  srand(mx->row+mx->col+my->col+iterations);

  initMatrix(&randomY);
  MatrixCopy(my, &randomY);

  /*Create a the r2q2scrambling matrix */
  /*
   * ccoeff_yscrambling columns:
   *  0: correlation between permuted and observed y
   *  1: r2 for the model
   *  2: q2 for the model
   */
  ResizeMatrix(ccoeff_yscrambling, 1+iterations, 3);
  NewDVector(&corrpermobs, my->col);
  NewDVector(&coeff_int, my->col);
  NewDVector(&coeff_valid, my->col);

  PermutedObservetCorrelation(my, my, yaverage, &corrpermobs);
  /*The first row is the model not scrambled...*/
  if(algo == _PLS_){
    PLSRegressionYScramblingPipeline(mx, my, input->xautoscaling, input->yautoscaling, input->nlv,  varg, nthreads, &coeff_int, &coeff_valid);
  }
  else if(algo == _PLS_DA_){
    PLSDiscriminantAnalysisYScramblingPipeline(mx, my, input->xautoscaling, input->yautoscaling, input->nlv, varg, nthreads, &coeff_int, &coeff_valid);
  }
  /*else if(algo_ == _MLR_)
  else if(algo_ == _LDA_)*/
  for(j = 0; j < my->col; j++){
    size_t mult = j*my->col;
    (*ccoeff_yscrambling)->data[0][mult] = corrpermobs->data[j];
    (*ccoeff_yscrambling)->data[0][mult+1] = coeff_int->data[j];
    (*ccoeff_yscrambling)->data[0][mult+2] = coeff_valid->data[j];
  }

  for(it = 0; it < iterations; it++){
    /* Sattolo's algorithm to shuffle the y*/
    i = randomY->row;
    while(i > 1){
      i = i - 1;
      //j = randInt(0, i);
      j = myrand_r(&srand_init) % i;
      for(k = 0; k < randomY->col; k++){
        ytmp = randomY->data[j][k];
        randomY->data[j][k] = randomY->data[i][k];
        randomY->data[i][k] = ytmp;
      }
    }

    PermutedObservetCorrelation(my, randomY, yaverage, &corrpermobs);
    /*The first row is the model not scrambled...*/
    if(algo == _PLS_){
      PLSRegressionYScramblingPipeline(mx, randomY, input->xautoscaling, input->yautoscaling, input->nlv, varg, nthreads, &coeff_int, &coeff_valid);
    }
    else if(algo == _PLS_DA_){
      PLSDiscriminantAnalysisYScramblingPipeline(mx, my, input->xautoscaling, input->yautoscaling, input->nlv, varg, nthreads, &coeff_int, &coeff_valid);
    }
    /*else if(algo_ == _MLR_)
    else if(algo_ == _LDA_)*/
    for(j = 0; j < my->col; j++){
      size_t mult = j*my->col;
      (*ccoeff_yscrambling)->data[it+1][mult] = corrpermobs->data[j];
      (*ccoeff_yscrambling)->data[it+1][mult+1] = coeff_int->data[j];
      (*ccoeff_yscrambling)->data[it+1][mult+2] = coeff_valid->data[j];
    }
  }

  DelDVector(&coeff_valid);
  DelDVector(&coeff_int);
  DelDVector(&corrpermobs);
  DelDVector(&yaverage);
  DelMatrix(&randomY);
}

void YScrambling_OLDMETHOD(MODELINPUT *input, AlgorithmType algo, ValidationArg varg, size_t blocks,
                  matrix **ccoeff_yscrambling, size_t nthreads, ssignal *s)
{
  size_t scrambiterations, iterations_, i, j, k, n, y_, blocksize;
  int id;
  double temp;

  matrix *mx, *my;
  matrix *randomY, *sorted_y_id, *sorty, *gid;
  dvector *yaverage;
  dvector *corrpermobs, *r2mod, *q2mod;
  mx = (*input->mx);
  my = (*input->my);

  initDVector(&yaverage);
  MatrixColAverage(my, &yaverage);

  srand(mx->row*mx->col*my->col*blocks);
  NewMatrix(&randomY, my->row, my->col);
  NewMatrix(&sorted_y_id, my->row, my->col);

  NewMatrix(&sorty, my->row, 2);
  for(j = 0; j < my->col; j++){
    for(i = 0; i < my->row; i++){
      sorty->data[i][0] = my->data[i][j];
      sorty->data[i][1] = i;
    }
    MatrixSort(sorty, 0);

    for(i = 0; i < my->row; i++){
      sorted_y_id->data[i][j] = sorty->data[i][1];
    }
  }
  DelMatrix(&sorty);


  /*calcualte the block size for the rotate matrix*/
  blocksize = (size_t)ceil(mx->row/(double)blocks);
  blocksize += (size_t)ceil((float)((blocksize*blocks) - mx->row)/  blocks);

  NewMatrix(&gid, blocks, blocksize);
  MatrixSet(gid, -2);
  /* Crate the boxes to fill -2 means no value to fill, -1 means value to fill*/
  for(i = 0, j = 0, k = 0; i < mx->row; i++){
    if(j < blocks){
      gid->data[j][k] = -1;
      j++;
    }
    else{
      j = 0;
      k++;
      gid->data[j][k] = -1;
      j++;
    }
  }

  /*get number of iterations*/
  scrambiterations = 0;
  for(i = 0; i < gid->row; i++){
    iterations_ = 0;
    for(j = 0; j < gid->col; j++){
      if((int)gid->data[i][j] == -1){
        iterations_++;
      }
      else{
        continue;
      }
    }

    if(iterations_ > scrambiterations){
      scrambiterations = iterations_;
    }
    else{
      continue;
    }
  }

  /*Create a the r2q2scrambling matrix */
  /*
   * ccoeff_yscrambling columns:
   *  0: correlation between permuted and observed y
   *  1: r2 for the model
   *  2: q2 for the model
   */
  ResizeMatrix(ccoeff_yscrambling, 1+scrambiterations, 3);
  NewDVector(&corrpermobs, my->col);
  NewDVector(&r2mod, my->col);
  NewDVector(&q2mod, my->col);

  PermutedObservetCorrelation(my, my, yaverage, &corrpermobs);
  /*The first row is the model not scrambled...*/
  if(algo == _PLS_){
    PLSRegressionYScramblingPipeline(mx, my, input->xautoscaling, input->yautoscaling, input->nlv,  varg, nthreads, &r2mod, &q2mod);
  }
  /*else if(algo_ == _MLR_)
  else if(algo_ == _LDA_)*/
  for(j = 0; j < my->col; j++){
    size_t mult = j*my->col;
    (*ccoeff_yscrambling)->data[0][mult] = corrpermobs->data[j];
    (*ccoeff_yscrambling)->data[0][mult+1] = r2mod->data[j];
    (*ccoeff_yscrambling)->data[0][mult+2] = q2mod->data[j];
  }


  for(y_ = 0; y_ < sorted_y_id->col; y_++){
    /* START WITH THE ORDERED Y_*/
    k = 0;
    for(i = 0; i < gid->row; i++){
      for(j = 0; j < gid->col; j++){
        if(gid->data[i][j] >= -1){
          gid->data[i][j] = sorted_y_id->data[k][y_];
          k++;
        }
        else{
          continue;
        }
      }
    }
    /*
    puts("GID Y_");
    PrintMatrix(gid);
    */

    iterations_ = 0;
    while(iterations_ <  scrambiterations){
      if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
        break;
      }
      else{
        /* Shuffle the Y bloks */
        for(i = 0; i < gid->row; i++){
          for(j = (gid->col-1); j > 0; j--){
            if(gid->data[i][j] > -1){
              /*then shift from this id to the first value..*/
              temp = gid->data[i][j];
              for(k = j; k > 0; k--){
                gid->data[i][k] = gid->data[i][k-1];
              }
              gid->data[i][0] = temp;
              break;
            }
            else{
              continue;
            }
          }
        }

        /*
        printf("Shuffled Y ID %d\n", (int)iterations_);
        PrintMatrix(gid);
        */

        /*Fill the shifted y*/
        n = 0;
        for(i = 0; i < gid->row; i++){
          for(j = 0; j < gid->col; j++){
            id = gid->data[i][j];
            if(id > -1){
              for(k = 0; k < my->col; k++){
                randomY->data[n][k] = my->data[id][k];
              }
              n++;
            }
            else{
              continue;
            }
          }
        }

        PermutedObservetCorrelation(my, randomY, yaverage, &corrpermobs);
        /*The first row is the model not scrambled...*/
        if(algo == _PLS_){
          PLSRegressionYScramblingPipeline(mx, randomY, input->xautoscaling, input->yautoscaling, input->nlv, varg, nthreads, &r2mod, &q2mod);
        }
        /*else if(algo_ == _MLR_)
        else if(algo_ == _LDA_)*/
        for(j = 0; j < my->col; j++){
          size_t mult = j*my->col;
          (*ccoeff_yscrambling)->data[iterations_+1][mult] = corrpermobs->data[j];
          (*ccoeff_yscrambling)->data[iterations_+1][mult+1] = r2mod->data[j];
          (*ccoeff_yscrambling)->data[iterations_+1][mult+2] = q2mod->data[j];
        }
        iterations_++;
      }
    }
  }

  DelDVector(&q2mod);
  DelDVector(&r2mod);
  DelDVector(&corrpermobs);
  DelDVector(&yaverage);
  DelMatrix(&sorted_y_id);
  DelMatrix(&randomY);
  DelMatrix(&gid);
}
