
/* Cross Validation
 *
 * 1) Divide the dataset in "g" random group
 * 2) for each randomization run:
 *    for each group thake this out and run a PLS model with all the remains groups
 * 3) Predict the response value for the out group and compute PRESS (y_pred - y_r)^2 and the Sum of square of (y_r - y_mid)^2
 * 4) for each randomization divide all the sum by the number of group and sum this untill all the computation are done.
 */


typedef struct{

  matrix *mx, *my; /*INPUT*/

  matrix *predicted_y;  /*OUPUT*/
  uivector *predictioncounter; /*OUPUT*/

  size_t xautoscaling, yautoscaling, nlv, group; /*INPUT*/
  unsigned int srand_init;
} rgcv_th_arg;


void *RandomGroupCVModel(void *arg_)
{
  size_t i, j, k, n, g, lv;
  rgcv_th_arg *arg;
  matrix *gid; /* randomization and storing id for each random group into a matrix */

  /* Matrix for compute the PLS models for groups */
  matrix *subX;
  matrix *subY;
  PLSMODEL *subm;

  /*matrix to predict*/
  matrix *predictX;
  matrix *realY;
  matrix *predicty;


  arg = (rgcv_th_arg*) arg_;

  NewMatrix(&gid, arg->group, (size_t)ceil(arg->mx->row/(double)arg->group));

  /* Divide in group  all the Dataset */
  MatrixSet(gid, -1);

  /* step 1 generate the random groups */
  k = 0;
  for(i = 0; i <  gid->row; i++){
    for(j = 0; j <  gid->col; j++){
      do{
        /*n = randInt(0, arg->mx->row);*/
        n = (size_t)myrand_r(&arg->srand_init) % (arg->mx->row);
      } while(ValInMatrix(gid, n) == 1 && k < (arg->mx->row));
      if(k < arg->mx->row){
        gid->data[i][j] = n;
        k++;
      }
      else
        continue;
    }
  }

  /*
  puts("Gid Matrix");
  PrintMatrix(gid);
  */
      /*
    printf("Excuded the group number %u\n", (unsigned int)g);
    puts("Sub Model\nX:");
    PrintArray(subX);
    puts("Y:");
    PrintArray(subY);

    puts("\n\nPredict Group\nX:");
    PrintArray(predictX);
    puts("RealY:");
    PrintArray(realY);
    */
  /*step 2*/
  for(g = 0; g < gid->row; g++){ /*For aeach group */
    /* Estimate how many objects are inside the sub model without the group "g" */
    n = 0;
    for(i = 0; i < gid->row; i++){
      if(i != g){
        for(j = 0; j < gid->col; j++){
          if((int)gid->data[i][j] != -1)
            n++;
          else
            continue;
        }
      }
      else
        continue;
    }

    /*Allocate the submodel*/
    NewMatrix(&subX, n, arg->mx->col);
    NewMatrix(&subY, n, arg->my->col);

    /* Estimate how many objects are inside the group "g" to predict*/
    n = 0;
    for(j = 0; j < gid->col; j++){
      if((int)gid->data[g][j] != -1)
        n++;
      else
        continue;
    }


    /*Allocate the */
    NewMatrix(&predictX, n, arg->mx->col);
    NewMatrix(&realY, n, arg->my->col);

    /* copy the submodel values */

    for(i = 0, k = 0; i < gid->row; i++){
      if(i != g){
        for(j = 0; j < gid->col; j++){
          size_t a =  (size_t)gid->data[i][j]; /* get the row index */
          if(a != -1){
            for(n = 0; n < arg->mx->col; n++){
              subX->data[k][n] = arg->mx->data[a][n];
            }
            for(n = 0; n < arg->my->col; n++){
              subY->data[k][n] = arg->my->data[a][n];
            }
            k++;
          }
          else{
            continue;
          }
        }
      }
      else{
        continue;
      }
    }

    /* copy the objects to predict into predictmx*/
    for(j = 0, k = 0; j < gid->col; j++){
      size_t a = (size_t)gid->data[g][j];
      if(a != -1){
        for(n = 0; n < arg->mx->col; n++){
          predictX->data[k][n] = arg->mx->data[a][n];
        }
        for(n = 0; n < arg->my->col; n++){
          realY->data[k][n] = arg->my->data[a][n];
        }
        k++;
      }
      else{
        continue;
      }
    }

    NewPLSModel(&subm);

    PLS(subX, subY, arg->nlv, arg->xautoscaling, arg->yautoscaling, subm, NULL);

    /* Predict Y for each latent variable */
    initMatrix(&predicty);

    for(lv = 1; lv <= arg->nlv; lv++){
      matrix *recalcy;
      matrix *recalcscores;

      initMatrix(&recalcy);
      initMatrix(&recalcscores);

      PLSScorePredictor(predictX, subm, lv, &recalcscores);
      PLSYPredictor(recalcscores, subm, lv, &recalcy);

      for(j = 0; j < recalcy->col; j++){
        dvector *tmp = getMatrixColumn(recalcy, j);
        MatrixAppendCol(&predicty, tmp);
        DelDVector(&tmp);
      }

      DelMatrix(&recalcy);
      DelMatrix(&recalcscores);
    }

    for(j = 0, k = 0; j < gid->col; j++){
      size_t a = (size_t)gid->data[g][j]; /*riga dell'oggetto....*/
      if(a != -1){
        arg->predictioncounter->data[a] += 1; /* this object was visited */
        /* updating y */
        for(n = 0; n < predicty->col; n++){
          arg->predicted_y->data[a][n] +=  predicty->data[k][n];
        }

        k++;
      }
      else{
        continue;
      }
    }

    DelMatrix(&predicty);
    DelPLSModel(&subm);
    DelMatrix(&subX);
    DelMatrix(&subY);
    DelMatrix(&predictX);
    DelMatrix(&realY);
  }
  DelMatrix(&gid);
  return 0;
}

void PLSRandomGroupsCV(matrix *mx, matrix *my, size_t xautoscaling, size_t yautoscaling, size_t nlv, size_t group, size_t iterations, matrix **q2y, matrix **sdep, matrix **bias, matrix **predicted_y, matrix **pred_residuals, size_t nthreads, ssignal *s)
{
  if(nlv > 0 && mx->row == my->row && group > 0 && iterations > 0){
    size_t th, iterations_, i, j;

    pthread_t *threads;

    dvector *ymean;
    uivector *predictcounter;
    matrix *sum_ypredictions;

    if(nlv > mx->col){
      nlv = mx->col;
    }

    NewMatrix(&sum_ypredictions, my->row, my->col*nlv); /* each component have my->col ypsilon */
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
          NewMatrix(&arg[th].predicted_y, my->row, my->col*nlv); /* each component have my->col ypsilon */
          NewUIVector(&arg[th].predictioncounter, my->row);
          pthread_create(&threads[th], NULL, RandomGroupCVModel, (void*) &arg[th]);
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

    /*Calculate the Q2 and SDEP and Bias */
    if(q2y != NULL)
      ResizeMatrix(q2y, nlv, my->col);

    if(sdep != NULL)
      ResizeMatrix(sdep, nlv, my->col);

    if(bias != NULL)
      ResizeMatrix(bias, nlv, my->col);

    if(q2y != NULL || sdep != NULL){
      initDVector(&ymean);
      MatrixColAverage(my, &ymean);

      for(size_t lv = 0; lv < nlv; lv++){
        for(j = 0; j < my->col; j++){
          double ssreg = 0.f;
          double sstot = 0.f;

          /*y = m x + k: bias is the m angular coefficient */
          /*double ypredaverage = 0.f; used to calculate the k */
          for(i = 0; i < my->row; i++){
            ssreg += square(my->data[i][j] - sum_ypredictions->data[i][my->col*lv+j]);
            sstot += square(my->data[i][j] - ymean->data[j]);
            /*ypredaverage = sum_ypredictions->data[i][j+lv]; */
          }

          if(bias != NULL){
            double sum_xi = 0.f, sum_yi = 0.f;
            /*ypredaverage /= (double)my->row;*/
            for(i = 0; i < my->row; i++){
              sum_yi+=(sum_ypredictions->data[i][j+lv]*(my->data[i][j]-ymean->data[j]));
              sum_xi+=(my->data[i][j]*(my->data[i][j]-ymean->data[j]));
            }

            (*bias)->data[lv][j] = fabs(1 - sum_yi/sum_xi);
            /*k = Y-(X*b);*/
          }

          if(q2y != NULL)
            (*q2y)->data[lv][j] = 1.f - (ssreg/sstot);

          if(sdep != NULL)
            (*sdep)->data[lv][j] = sqrt(ssreg/(double)my->row);


        }
      }

      DelDVector(&ymean);
    }

    DelMatrix(&sum_ypredictions);
    DelUIVector(&predictcounter);
    xfree(threads);
    xfree(arg);
  }
  else{
    fprintf(stderr, "Error!! Unable to compute PLS Random Group Cross Validation!!\n");
  }
}


/* Leave One Ouut Validation
 *
 * 1) remove one object
 * 2) calculate the model
 * 3) predict the removed object and so the r2 and q2
 */

typedef struct{
  matrix *submx, *submy, *predmx, *predmy, *pred_y;
  size_t nlv, xautoscaling, yautoscaling;
} loocv_th_arg;

void *PLSLOOModel(void *arg_)
{
  size_t lv, j;
  loocv_th_arg *arg;
  arg = (loocv_th_arg*) arg_;

  PLSMODEL *subm;
  NewPLSModel(&subm);

  PLS(arg->submx, arg->submy, arg->nlv, arg->xautoscaling, arg->yautoscaling, subm, NULL);

  /* Predict Y for each latent variable */
  for(lv = 1; lv <= arg->nlv; lv++){
    matrix *recalcy;
    matrix *recalcscores;

    initMatrix(&recalcy);
    initMatrix(&recalcscores);

    PLSScorePredictor(arg->predmx, subm, lv, &recalcscores);
    PLSYPredictor(recalcscores, subm, lv, &recalcy);

    for(j = 0; j < recalcy->col; j++){
      dvector *tmp = getMatrixColumn(recalcy, j);
      MatrixAppendCol(&arg->pred_y, tmp);
      DelDVector(&tmp);
    }

    DelMatrix(&recalcy);
    DelMatrix(&recalcscores);
  }

  DelPLSModel(&subm);
  return 0;
}

void PLSLOOCV(matrix* mx, matrix* my, size_t xautoscaling, size_t yautoscaling, size_t nlv, matrix** q2y, matrix** sdep, matrix **bias, matrix** predicted_y, matrix **pred_residuals, size_t nthreads, ssignal *s)
{
 if(nlv > 0 && mx->row == my->row){
    size_t i, j, k, l, th, model;
    pthread_t *threads;
    loocv_th_arg *arg;

    if(nlv > mx->col){
      nlv = mx->col;
    }

    dvector *ymean;
    matrix *loopredictedy;
    NewMatrix(&loopredictedy, my->row, my->col*nlv);

    threads = xmalloc(sizeof(pthread_t)*nthreads);
    arg = xmalloc(sizeof(loocv_th_arg)*nthreads);

    /* initialize threads arguments.. */
    for(th = 0; th < nthreads; th++){
      arg[th].nlv = nlv;
      arg[th].xautoscaling = xautoscaling;
      arg[th].yautoscaling = yautoscaling;
      NewMatrix(&arg[th].submx, mx->row-1, mx->col);
      NewMatrix(&arg[th].submy, my->row-1, my->col);
      NewMatrix(&arg[th].predmx, 1, mx->col);
      NewMatrix(&arg[th].predmy, 1, my->col);
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
                  arg[th].submx->data[l][k] = mx->data[j][k];
                }
                for(k = 0; k < my->col; k++){
                  /*setMatrixValue(arg[th].submy, l, k, getMatrixValue(my, j, k));*/
                  arg[th].submy->data[l][k] = my->data[j][k];
                }
                l++;
              }
              else{
                for(k = 0; k < mx->col; k++){
                  arg[th].predmx->data[0][k] = mx->data[j][k];
                }
                for(k = 0; k < my->col; k++){
                  arg[th].predmy->data[0][k] = my->data[j][k];
                }
              }
            }

            initMatrix(&arg[th].pred_y);
            pthread_create(&threads[th], NULL, PLSLOOModel, (void*) &arg[th]);
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
            for(j = 0; j < arg[th].pred_y->col; j++){
              loopredictedy->data[model+th][j] = arg[th].pred_y->data[0][j];
            }
            DelMatrix(&arg[th].pred_y);
          }
        }
      }
    }

    /*Delete thread arguments*/

    for(th = 0; th < nthreads; th++){
      DelMatrix(&arg[th].submx);
      DelMatrix(&arg[th].submy);
      DelMatrix(&arg[th].predmx);
      DelMatrix(&arg[th].predmy);
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

    /*Calculate the Q2 and SDEP */
    if(q2y != NULL)
      ResizeMatrix(q2y, nlv, my->col);

    if(sdep != NULL)
      ResizeMatrix(sdep, nlv, my->col);

    if(bias != NULL)
      ResizeMatrix(bias, nlv, my->col);

    if(sdep != NULL || q2y != NULL){
      initDVector(&ymean);
      MatrixColAverage(my, &ymean);

      for(size_t lv = 0; lv < nlv; lv++){
        for(j = 0; j < my->col; j++){
          double ssreg = 0.f;
          double sstot = 0.f;
          for(i = 0; i < my->row; i++){
            ssreg += square(loopredictedy->data[i][my->col*lv+j] - my->data[i][j]);
            sstot += square(my->data[i][j] - ymean->data[j]);
          }

          if(bias != NULL){
            double sum_yi = 0.f, sum_xi = 0.f;
            /*ypredaverage /= (double)my->row;*/
            for(i = 0; i < my->row; i++){
              sum_yi+=(loopredictedy->data[i][my->col*lv+j]*(my->data[i][j]-ymean->data[j]));
              sum_xi+=(my->data[i][j]*(my->data[i][j]-ymean->data[j]));
            }
            /*sum_yi/sum_xi = m */
            (*bias)->data[lv][j] = fabs(1 - sum_yi/sum_xi);
            /*k = Y-(X*b);*/
          }

          if(q2y != NULL)
            (*q2y)->data[lv][j] = 1.f - (ssreg/sstot);

          if(sdep != NULL)
            (*sdep)->data[lv][j] = sqrt(ssreg/(double)my->row);
        }

      }
      DelDVector(&ymean);
    }

    DelMatrix(&loopredictedy);
    xfree(arg);
    xfree(threads);
  }
  else{
    fprintf(stderr, "Error!! Unable to compute PLS Leave One Out Validation!!\n");
  }
}


void PLSStaticSampleValidator(matrix *mx, matrix *my, uivector *obj_class,
                        size_t xautoscaling, size_t yautoscaling,
                        size_t nlv, size_t sample_size, size_t niters,
                        size_t rgcv_group, size_t rgcv_iterations, size_t nthreads,
                        matrix **q2_distr, matrix **sdep_distr, uivector **bestid, ssignal *s)
{
  if(obj_class->size == mx->row){
    size_t iter, i, j, nclass = 0, srand_init;
    matrix *xsample, *ysample, *tmpq2, *tmpsdep;
    uivector *tmpid;
    double bestq2 = -9999;
    for(i = 0; i < obj_class->size; i++){
      if(obj_class->data[i] > nclass){
        nclass = obj_class->data[i];
      }
      else{
        continue;
      }
    }

    initMatrix(&tmpq2);
    initMatrix(&tmpsdep);

    srand_init = mx->row+mx->col+my->col+sample_size+niters+nclass;

    NewMatrix(&xsample, sample_size, mx->col);
    NewMatrix(&ysample, sample_size, my->col);
    if(bestid != NULL)
      UIVectorResize(bestid, sample_size);
    NewUIVector(&tmpid, sample_size);
    for(iter = 0; iter < niters; iter++){
      if(nclass == 0){ /* simple bootstrap on all objects */
        for(i = 0; i < sample_size; i++){
          do{
            size_t n = (size_t)myrand_r((unsigned *)&srand_init) % (mx->row);
            tmpid->data[i] = n;
            if(n > obj_class->size){
              continue;
            }
            else{
              for(j = 0; j < mx->col; j++){
                xsample->data[i][j] = mx->data[n][j];
              }

              for(j = 0; j < my->col; j++){
                ysample->data[i][j] = my->data[n][j];
              }
              break;
            }
          }while(1);
        }
      }
      else{ /* select a random object from each class */
        size_t a_class = 0;
        for(i = 0; i < sample_size; i++){
          do{
            size_t n = (size_t)myrand_r((unsigned *)&srand_init) % (mx->row);
            tmpid->data[i] = n;
            if(n > obj_class->size){
              continue;
            }
            else{
              if(obj_class->data[n] == a_class){
                if(a_class == nclass){
                  a_class =  (size_t)myrand_r((unsigned *)&srand_init) % (nclass); /*all class are already completed then restart random...*/
                }
                else{
                  a_class++;
                }

                for(j = 0; j < mx->col; j++){
                  xsample->data[i][j] = mx->data[n][j];
                }

                for(j = 0; j < my->col; j++){
                  ysample->data[i][j] = my->data[n][j];
                }
                break;
              }
              else{
                continue;
              }
            }
          }while(1);
        }
      }

      /*Now Compute validation...*/
      matrix *q2y;
      matrix *sdep;
      initMatrix(&q2y);
      initMatrix(&sdep);

      if(rgcv_group == 0 || rgcv_iterations == 0){
        PLSLOOCV(xsample, ysample, xautoscaling, yautoscaling, nlv,
                  &q2y,
                  &sdep,
                  NULL,
                  NULL,
                  NULL, nthreads, s);
      }
      else{
        PLSRandomGroupsCV(xsample, ysample, xautoscaling, yautoscaling, nlv, rgcv_group, rgcv_iterations,
                      &q2y,
                      &sdep,
                      NULL,
                      NULL,
                      NULL, nthreads, s);
      }

      if(bestid != NULL){
        size_t pc = GetLVCCutoff(q2y);
        double cumq2 = 0.f;
        for(i = 0; i < q2y->col; i++)
          cumq2 += q2y->data[pc][i];

        if(cumq2 > bestq2){
          bestq2 = cumq2;
          for(i = 0; i < tmpid->size; i++)
            (*bestid)->data[i] = tmpid->data[i];
        }
      }

      //Search for the best Q2 Model and Select ID to return
      for(i = 0; i < q2y->col; i++){
        dvector *q2row = getMatrixColumn(q2y, i);
        MatrixAppendRow(&tmpq2, q2row);
        DelDVector(&q2row);

        dvector *sdep_row = getMatrixColumn(sdep, i);
        MatrixAppendRow(&tmpsdep, sdep_row);
        DelDVector(&sdep_row);
      }

      DelMatrix(&q2y);
      DelMatrix(&sdep);
    }
    DelMatrix(&xsample);
    DelMatrix(&ysample);

    /*Finalize output */
    ResizeMatrix(q2_distr, tmpq2->col, tmpq2->row);
    ResizeMatrix(sdep_distr, tmpsdep->col, tmpsdep->row);

    MatrixTranspose(tmpq2, (*q2_distr));
    MatrixTranspose(tmpsdep, (*sdep_distr));
    DelMatrix(&tmpsdep);
    DelMatrix(&tmpq2);
    DelUIVector(&tmpid);
  }
  else{
    fprintf(stderr, "Error!! The number of objects differ from the number of class objects\n");
  }
}

void PLSDynamicSampleValidator(matrix *mx, matrix *my,
                        size_t xautoscaling, size_t yautoscaling,
                        size_t nlv, size_t niters,
                        uivector *obj_class, size_t deltaobj, size_t maxobj,
                        size_t rgcv_group, size_t rgcv_iterations, size_t nthreads,
                        matrix **q2_surface, matrix **sdep_surface, uivector **bestid, ssignal *s)
{
  if(obj_class->size == mx->row && mx->row > 10){
    size_t iter, i, j, k, l, nclass = 0, srand_init, incobj;
    matrix *xsample, *ysample;
    uivector *tmpflag, *tmpbestflag;
    double bestq2 = -9999;
    for(i = 0; i < obj_class->size; i++){
      if(obj_class->data[i] > nclass){
        nclass = obj_class->data[i];
      }
      else{
        continue;
      }
    }

    if(maxobj > mx->row)
      maxobj = mx->row;

    NewMatrix(q2_surface, ceil(maxobj/deltaobj)*niters*nlv, 2+my->col);
    NewMatrix(sdep_surface, ceil(maxobj/deltaobj)*niters*nlv, 2+my->col);

    /*first column: nobj
     second column: nlv
     from third to the end each serie of my->col is a q2 for each y*/

    srand_init = mx->row+mx->col+my->col+niters+nclass+maxobj;

    if(bestid != NULL){
      NewUIVector(&tmpflag, mx->row);
      NewUIVector(&tmpbestflag, mx->row);
    }

    l = 0;
    for(iter = 0; iter < niters; iter++){
      incobj = 10;
      for(k = 0; k < ceil(maxobj/deltaobj); k++){
        NewMatrix(&xsample, incobj, mx->col);
        NewMatrix(&ysample, incobj, my->col);

        if(nclass == 0){ /* simple bootstrap on all objects */
          for(i = 0; i < incobj; i++){
            do{
              size_t n = (size_t)myrand_r((unsigned *)&srand_init) % (mx->row);
              if(n > mx->row){
                continue;
              }
              else{
                if(bestid != NULL)
                  tmpflag->data[n] = 1;

                for(j = 0; j < mx->col; j++){
                  xsample->data[i][j] = mx->data[n][j];
                }

                for(j = 0; j < my->col; j++){
                  ysample->data[i][j] = my->data[n][j];
                }
                break;
              }
            }while(1);
          }
        }
        else{ /* select a random object from each class */
          size_t a_class = 0;
          for(i = 0; i < incobj; i++){
            do{
              size_t n = (size_t)myrand_r((unsigned *)&srand_init) % (mx->row);
              if(n > mx->row){
                continue;
              }
              else{
                if(obj_class->data[n] == a_class){
                  if(a_class == nclass){
                    a_class =  (size_t)myrand_r((unsigned *)&srand_init) % (nclass); /*all class are already completed then restart random...*/
                  }
                  else{
                    a_class++;
                  }

                  if(bestid != NULL)
                    tmpflag->data[n] = 1;

                  for(j = 0; j < mx->col; j++){
                    xsample->data[i][j] = mx->data[n][j];
                  }

                  for(j = 0; j < my->col; j++){
                    ysample->data[i][j] = my->data[n][j];
                  }
                  break;
                }
                else{
                  continue;
                }
              }
            }while(1);
          }
        }

        /*Now Compute validation...*/
        matrix *q2y;
        matrix *sdep;
        initMatrix(&q2y);
        initMatrix(&sdep);

        if(rgcv_group == 0 || rgcv_iterations == 0){
          PLSLOOCV(xsample, ysample, xautoscaling, yautoscaling, nlv,
                    &q2y,
                    &sdep,
                    NULL,
                    NULL,
                    NULL, nthreads, s);
        }
        else{
          PLSRandomGroupsCV(xsample, ysample, xautoscaling, yautoscaling, nlv, rgcv_group, rgcv_iterations,
                        &q2y,
                        &sdep,
                        NULL,
                        NULL,
                        NULL, nthreads, s);
        }

        if(bestid != NULL){
          size_t pc = GetLVCCutoff(q2y);
          double cumq2 = 0.f;
          for(i = 0; i < q2y->col; i++)
            cumq2 += q2y->data[pc][i];

          if(cumq2 > bestq2){
            bestq2 = cumq2;
            for(i = 0; i < tmpflag->size; i++){
              tmpbestflag->data[i] = tmpflag->data[i];
              tmpflag->data[i] = 0;
            }
          }
        }

        for(i = 0; i < q2y->row; i++){
          (*q2_surface)->data[l][0] = incobj;
          (*q2_surface)->data[l][1] = i+1;
          (*sdep_surface)->data[l][0] = incobj;
          (*sdep_surface)->data[l][1] = i+1;
          for(j = 0; j < q2y->col; j++){
            if(q2y->data[i][j] < -1){
              (*q2_surface)->data[l][j+2] = 0;
            }
            else
              (*q2_surface)->data[l][j+2] = q2y->data[i][j];
            (*sdep_surface)->data[l][j+2] = sdep->data[i][j];
          }
          l++;
        }

        if(incobj + deltaobj < mx->row)
          incobj += deltaobj;
        else
          incobj = mx->row;

        DelMatrix(&q2y);
        DelMatrix(&sdep);
        DelMatrix(&xsample);
        DelMatrix(&ysample);
      }
    }

    if(bestid != NULL){

      for(i = 0, j = 0; i < tmpbestflag->size; i++){
        if(tmpbestflag->data[i] == 1)
          j++;
        else
          continue;
      }

      UIVectorResize(bestid, j);
      for(i = 0, j = 0; i < tmpbestflag->size; i++){
        if(tmpbestflag->data[i] == 1){
          (*bestid)->data[j] = i;
          j++;
        }
        else
          continue;
      }

      DelUIVector(&tmpflag);
      DelUIVector(&tmpbestflag);
    }
  }
  else{
    if(mx->row > 10)
      fprintf(stderr, "Error!! Insufficient number of objects\n");
    else
      fprintf(stderr, "Error!! The number of objects differ from the number of class objects\n");
  }
}
