#include "modelvalidation.h"
#include "memwrapper.h"
#include "mlr.h"
#include "pls.h"
#include "pca.h" /*Using: MatrixAutoScaling(); and calcVarExpressed(); */
#include "numeric.h" /* Using:  if(FLOAT_EQ(NumOne, NumTwo));*/
#include "metricspace.h"
#include <math.h>
#include <pthread.h>

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
  matrix *mx, *my; /*INPUT*/
  matrix *predicted_y;  /*OUPUT*/
  uivector *predictioncounter; /*OUPUT*/
  size_t xautoscaling, yautoscaling, nlv, group; /*INPUT*/
  unsigned int srand_init;
} rgcv_th_arg;

void *PLSRandomGroupCVModel(void *arg_)
{
  size_t j, k, n, g, lv;
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

    /* Predict Y for each latent variable */
    //initMatrix(&predicty);
    initMatrix(&y_test_predicted);

    for(lv = 1; lv <= arg->nlv; lv++){
      matrix *predicted_test_y;
      matrix *predicted_test_scores;

      initMatrix(&predicted_test_y);
      initMatrix(&predicted_test_scores);

      PLSScorePredictor(x_test, subm, lv, &predicted_test_scores);
      PLSYPredictor(predicted_test_scores, subm, lv, &predicted_test_y);

      for(j = 0; j < predicted_test_y->col; j++){
        dvector *tmp = getMatrixColumn(predicted_test_y, j);
        MatrixAppendCol(&y_test_predicted, tmp);
        DelDVector(&tmp);
      }

      DelMatrix(&predicted_test_y);
      DelMatrix(&predicted_test_scores);
    }

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

//void PLSRandomGroupsCV(matrix *mx, matrix *my, size_t xautoscaling, size_t yautoscaling, size_t nlv, size_t group, size_t iterations, matrix **q2y, matrix **sdep, matrix **bias, matrix **predicted_y, matrix **pred_residuals, size_t nthreads, ssignal *s)
void BootstrapRandomGroupsCV(MODELINPUT *input, size_t group, size_t iterations, AlgorithmType algo, matrix **predicted_y, matrix **pred_residuals, size_t nthreads, ssignal *s)
{
  matrix *mx = (*input->mx);
  matrix *my = (*input->my);
  size_t nlv = input->nlv;
  size_t xautoscaling = input->xautoscaling;
  size_t yautoscaling = input->yautoscaling;

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
          arg[th].mx = mx;
          arg[th].my = my;
          arg[th].group = group;
          arg[th].nlv = nlv;
          arg[th].xautoscaling = xautoscaling;
          arg[th].yautoscaling = yautoscaling;
          arg[th].srand_init = (unsigned int) group + mx->row + my->col + iterations + th + iterations_;
          NewMatrix(&arg[th].predicted_y, my->row, scol);
          NewUIVector(&arg[th].predictioncounter, my->row);
          if(algo == _PLS_){
            pthread_create(&threads[th], NULL, PLSRandomGroupCVModel, (void*) &arg[th]);
          }
          if(algo == _MLR_){
            pthread_create(&threads[th], NULL, MLRRandomGroupCVModel, (void*) &arg[th]);
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
  size_t lv, j;
  loocv_th_arg *arg;
  arg = (loocv_th_arg*) arg_;

  PLSMODEL *subm;
  NewPLSModel(&subm);

  PLS(arg->x_train, arg->y_train, arg->nlv, arg->xautoscaling, arg->yautoscaling, subm, NULL);

  /* Predict Y for each latent variable */
  for(lv = 1; lv <= arg->nlv; lv++){
    matrix *predicted_test_y;
    matrix *predicted_test_scores;

    initMatrix(&predicted_test_y);
    initMatrix(&predicted_test_scores);

    PLSScorePredictor(arg->x_test, subm, lv, &predicted_test_scores);
    PLSYPredictor(predicted_test_scores, subm, lv, &predicted_test_y);

    for(j = 0; j < predicted_test_y->col; j++){
      dvector *tmp = getMatrixColumn(predicted_test_y, j);
      MatrixAppendCol(&arg->y_test_predicted, tmp);
      DelDVector(&tmp);
    }

    DelMatrix(&predicted_test_y);
    DelMatrix(&predicted_test_scores);
  }

  DelPLSModel(&subm);
  return 0;
}

void LeaveOneOut(MODELINPUT *input, AlgorithmType algo, matrix** predicted_y, matrix **pred_residuals, size_t nthreads, ssignal *s)
{
  matrix *mx = (*input->mx);
  matrix *my = (*input->my);
  size_t nlv = input->nlv;
  size_t xautoscaling = input->xautoscaling;
  size_t yautoscaling = input->yautoscaling;

  if(mx->row == my->row){
    size_t i, j, k, l, th, model;
    pthread_t *threads;
    loocv_th_arg *arg;

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


    matrix *loopredictedy;
    NewMatrix(&loopredictedy, my->row, scol);

    threads = xmalloc(sizeof(pthread_t)*nthreads);
    arg = xmalloc(sizeof(loocv_th_arg)*nthreads);

    /* initialize threads arguments.. */
    for(th = 0; th < nthreads; th++){
      arg[th].nlv = nlv;
      arg[th].xautoscaling = xautoscaling;
      arg[th].yautoscaling = yautoscaling;
      NewMatrix(&arg[th].x_train, mx->row-1, mx->col);
      NewMatrix(&arg[th].y_train, my->row-1, my->col);
      NewMatrix(&arg[th].x_test, 1, mx->col);
      NewMatrix(&arg[th].y_test, 1, my->col);
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
                  arg[th].x_train->data[l][k] = mx->data[j][k];
                }
                for(k = 0; k < my->col; k++){
                  /*setMatrixValue(arg[th].submy, l, k, getMatrixValue(my, j, k));*/
                  arg[th].y_train->data[l][k] = my->data[j][k];
                }
                l++;
              }
              else{
                for(k = 0; k < mx->col; k++){
                  arg[th].x_test->data[0][k] = mx->data[j][k];
                }
                for(k = 0; k < my->col; k++){
                  arg[th].y_test->data[0][k] = my->data[j][k];
                }
              }
            }

            initMatrix(&arg[th].y_test_predicted);
            if(algo == _PLS_)
              pthread_create(&threads[th], NULL, PLSLOOModel_, (void*) &arg[th]);
            else if(algo == _MLR_)
              pthread_create(&threads[th], NULL, MLRLOOModel_, (void*) &arg[th]);
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
            for(j = 0; j < arg[th].y_test_predicted->col; j++){
              loopredictedy->data[model+th][j] = arg[th].y_test_predicted->data[0][j];
            }
            DelMatrix(&arg[th].y_test_predicted);
          }
        }
      }
    }

    /*Delete thread arguments*/

    for(th = 0; th < nthreads; th++){
      DelMatrix(&arg[th].x_train);
      DelMatrix(&arg[th].y_train);
      DelMatrix(&arg[th].x_test);
      DelMatrix(&arg[th].y_test);
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
    xfree(arg);
    xfree(threads);
  }
  else{
    fprintf(stderr, "Error!! Unable to compute PLS Leave One Out Validation!!\n");
  }
}
