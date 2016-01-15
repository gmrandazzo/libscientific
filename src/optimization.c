#include "optimization.h"


void SIMPLEX(matrix* mx, matrix* my, size_t npc, size_t xautoscaling, size_t yautoscaling, size_t vtype, size_t ngroup, size_t niter, uivector* var_on, uivector *modobjects, double error, double resolution, size_t nloop, size_t operator_, matrix** scoeff, ssignal* s)
{
  size_t i, j, loop;
  matrix *mx_;
  /* dvector *coefficients; */
  matrix *coefficients;

  NewMatrix(&mx_, mx->row, mx->col);

  loop = 0;
  while(loop < nloop){
    if(loop == 0){
      initMatrix(&coefficients);
      for(i = 0; i < mx->row; i++){
        for(j = 0; j < mx->col; j++){
          setMatrixValue(mx_, i, j, getMatrixValue(mx, i, j));
        }
      }

      if(modobjects == NULL || modobjects->size == 0){
        /*LOO OVERFITTING*/
        SINGLESIMPLS(mx_, my, npc, xautoscaling, yautoscaling,
                              0, ngroup, niter,
                              var_on, error, resolution, operator_, coefficients, NULL);
      }
      else{
        /*LOO OVERFITTING*/
        SINGLESIMPLS_OBJMOD(mx_, my, npc, xautoscaling, yautoscaling,
                              0, ngroup, niter,
                              var_on, modobjects, resolution, operator_, coefficients, NULL);
      }


      /*RGCV OVERFITTING
      SIMPLS(x, y, 5, 1, 0,
                            1, 5, 20,
                            var_on, 80, coefficients, NULL);
                            */

      MatrixCopy(coefficients, scoeff);
    }
    else{
      for(i = 0; i < mx_->row; i++){
        for(j = 0; j < mx_->col; j++){
          setMatrixValue(mx_, i, j, getMatrixValue(mx_, i, j)+getMatrixValue(coefficients, i, j));
          setMatrixValue((*scoeff), i, j, getMatrixValue((*scoeff), i ,j) + getMatrixValue(coefficients, i, j));
        }
      }

      DelMatrix(&coefficients);
      initMatrix(&coefficients);
      if(modobjects == NULL){
        SINGLESIMPLS(mx_, my, npc, xautoscaling, yautoscaling,
                            0, ngroup, niter,
                            var_on, error, resolution, operator_, coefficients, NULL);
      }
      else{
          SINGLESIMPLS_OBJMOD(mx_, my, npc, xautoscaling, yautoscaling,
                            0, ngroup, niter,
                            var_on, modobjects, resolution, operator_, coefficients, NULL);
      }
    }
    loop++;
  }

  for(i = 0; i < mx->row; i++){
    for(j = 0; j < mx->col; j++){
      setMatrixValue((*scoeff), i, j, getMatrixValue((*scoeff), i ,j) + getMatrixValue(coefficients, i, j));
    }
  }

  DelMatrix(&coefficients);
  DelMatrix(&mx_);
}

void SINGLESIMPLS(matrix *mx, matrix *my, size_t npc, size_t xautoscaling, size_t yautoscaling,
                           size_t vtype, size_t ngroup, size_t niter,
                           uivector *var_on, double error, double resolution, size_t operator_, matrix *coefficients_, ssignal *s)
{
  size_t obj, i, j, k, a, loop, mod_cutoff;
  double n, d, yerror, new_yerror, step, coeff, res, worst, second_worst, best, old_simplex_worst, old_simplex_second_worst, old_simplex_best;
  dvector *new_response, *colaverage, *sdeprowaverage, *tmp;
  matrix *coefficients, *submx, *submy, *predmx, *predscore, *predmy, *predictedy, *response, *mod_sdep, *mod_pred_y;

  PLSMODEL *subm;

  /*var_on->siye == mx->col*/
  NewMatrix(&response, var_on->size+1, var_on->size+1); /* +1 because the last is the response value */
  NewMatrix(&coefficients, mx->row, var_on->size+1); /*row vector to store the best response*/

  NewDVector(&new_response, var_on->size); /*row vector to store the best response*/

  /*
  initDVector(&colaverage);
  MatrixColAverage(mx, &colaverage);
  */
  NewMatrix(&predmx, 1, mx->col); /* 1 because we have only an object to make a good prediction */
  NewMatrix(&predmy, 1, my->col);

  /* Search outlier and do simplex for each outlier... */
  initMatrix(&mod_sdep);
  initMatrix(&mod_pred_y);

  if(vtype == 0){
    PLSLOOCV(mx, my, xautoscaling, yautoscaling, npc, NULL, NULL, &mod_sdep, &mod_pred_y, NULL, s);
  }
  else{
    PLSRandomGroupsCV(mx, my, xautoscaling, yautoscaling, npc, ngroup, niter, NULL, NULL, &mod_sdep, &mod_pred_y, NULL, s);
  }

  /*get the PC with a low error: mod_cutoff

  if only one y
    MatrixGetMinValue(sdep, &mod_cutoff, NULL);
  else
    get the medium error value
  */

  initDVector(&sdeprowaverage);
  MatrixRowAverage(mod_sdep, &sdeprowaverage); /*medium of sdep if more than one y*/
  mod_cutoff = 0;
  for(i = 1; i < sdeprowaverage->size; i++){
    if(getDVectorValue(sdeprowaverage, i) < getDVectorValue(sdeprowaverage, mod_cutoff)){
      mod_cutoff = i;
    }
    else{
      continue;
    }
  }

  /*
  puts("Sdep Row Average");
  PrintDVector(sdeprowaverage);

  printf("Best SDEP %f\n", getDVectorValue(sdeprowaverage, mod_cutoff));
  */
  DelDVector(&sdeprowaverage);

  /*printf("best model component... %lu\n", mod_cutoff+1);*/

  initMatrix(&submx);
  initMatrix(&submy);
  for(obj = 0; obj < my->row; obj++){
  /*
  Prediction for example are organized as this example
  y = 2
  number_of_pc = 4

      pc1    pc2     pc3       pc4
    |  |   |   |   |   |    |    |
    y1 y2  y1  y2  y1  y2   y1   y2
ID  0 1    2   3   4   5    6    7

so:  ID = a * n_y
    a = PC
    n_y = id of y

    Example get y in PC3:
    a = 2
    n_y1 = 0
    n_y2 = 1
    pc3 y1 = 4 + 0 = 4
    pc3 y2 = 4 + 1 = 5
  */
    a = mod_cutoff * my->col;
    /*
    printf("Get the predicted y in matrix with id %lu \n", a);
    */

    /*sum the error of all y and make a media of error*/
    yerror = 0.f;
    n = d = 0.f;
    for(j = 0; j < my->col; j++){
      n += sqrt(square(getMatrixValue(my, obj, j) - getMatrixValue(mod_pred_y, obj, a+j)));
      d += getMatrixValue(my, obj, j);
    }
    n /= (double)my->col;
    d /= (double)my->col;
    yerror = (n/d) * 100;

    /*printf("The ERROR for the %lu object is %lf > %lf?\n", obj, yerror, error);*/

    if(yerror > error){
      continue;
    }
    else{
      tmp = getMatrixRow(mx, obj);
      MatrixAppendRow(&submx, tmp);
      DelDVector(&tmp);
      tmp = getMatrixRow(my, obj);
      MatrixAppendRow(&submy, tmp);
      DelDVector(&tmp);
    }
  }

  /*
  puts("Matrix Training model");
  PrintMatrix(submx);
  PrintMatrix(submy);
 */
  if(submx->row > 0){
    /*build  submodel in order to do the training  for the other outliers... */
    initDVector(&colaverage);
    MatrixColAverage(submx, &colaverage);

    NewPLSModel(&subm);
    PLS(submx, submy, npc, xautoscaling, yautoscaling, subm, s);

    for(obj = 0; obj < my->row; obj++){
      /*puts("\n #### SEARCH FOR NEW OBJECT #### \n");*/
      a = mod_cutoff * my->col;

      /*printf("sum for get the y in prediction matrix %lu \n", a);*/

      /*sum the error of all y and make a media of error*/
      yerror = 0.f;
      n = d = 0.f;
      for(j = 0; j < my->col; j++){
        n += sqrt(square(getMatrixValue(my, obj, j) - getMatrixValue(mod_pred_y, obj, a+j)));
        d += getMatrixValue(my, obj, j);
      }
      n /= (double)my->col;
      d /= (double)my->col;
      yerror = (n/d) * 100;

      if(yerror > error){
        /*
        printf("Object  out of error %lu\n", obj);
        printf("Error %lf\n", yerror);
        */

        /*Perform k+1 experiments on the vertice of a simplex
        * Build the Simplex matrix experiment with the response value as last column of response
        *
        *  coeff = 0.f if we sum in new_mx
        *  coeff = 1.f if we moltiplicate in new_mx
        */

        if(operator_ == 0){ /*sum*/
          coeff = 0.f;
        }
        else{
          coeff = 1.f;
        }

        for(k = 0; k < response->row; k++){
          if(k == 0){
            for(j = 0; j < response->col-1; j++){/*-1 because the last column is the score*/
              if(1 == getUIVectorValue(var_on, j)){
                setMatrixValue(response, k, j, coeff);
              }
              else{
                continue;
              }
            }
          }
          else{
            for(j = 0; j < response->col-1; j++){
              if(1 == getUIVectorValue(var_on, j)){
                step = (getDVectorValue(colaverage, j)*resolution) / 100.f;
                if(j < k-1){
                  setMatrixValue(response, k, j, (step)/2.);
                }
                else if(j == k-1){
                  setMatrixValue(response, k, j, step);
                }
                else{
                  continue;
                }
              }
              else{
                setMatrixValue(response, k, j, coeff);
              }
            }
          }
        }

        /*
        puts("Simplex response matrix coefficients");
        PrintMatrix(response);
        */

        /* CALCULATE THE FIRST CICLE OF RESPONSE AND GET THE LOW ERROR IN PREDICTION*/
        for(k = 0; k < response->row; k++){

          /* build an prediction matrix in predmx, predmy with the simplex optimization... */
          for(j = 0; j < mx->col; j++){
            if(operator_ == 0){ /*sum*/
              setMatrixValue(predmx, 0, j, getMatrixValue(mx, obj, j) + getMatrixValue(response, k, j));
            }
            else{
              setMatrixValue(predmx, 0, j, getMatrixValue(mx, obj, j) * getMatrixValue(response, k, j));
            }
          }

          for(j = 0; j < my->col; j++){
            setMatrixValue(predmy, 0, j, getMatrixValue(my, obj, j));
          }

          /*
          puts("New Matrix");
          PrintMatrix(predmx);
          PrintMatrix(predmy);
          */
          initMatrix(&predictedy);
          initMatrix(&predscore);
          PLSScorePredictor(predmx, subm, npc, &predscore);
          PLSYPredictor(predscore, subm, mod_cutoff, &predictedy);


          new_yerror = 0.f;
          n = d = 0.f;
          for(j = 0; j < my->col; j++){
            n += sqrt(square(getMatrixValue(predmy, 0, j) - getMatrixValue(predictedy, 0, j)));
            d += getMatrixValue(predmy, 0, j);
          }
          n /= (double)my->col;
          d /= (double)my->col;
          new_yerror = n / d;

          setMatrixValue(response, k, response->col-1, new_yerror); /* yerror is the new response */
          DelMatrix(&predscore);
          DelMatrix(&predictedy);
        }

        MatrixSort(response, response->col-1); /*the worst response is at row col-1 and the position were wil be builded the ner  response*/
        old_simplex_worst = getMatrixValue(response, response->row-1, response->col-1);
        old_simplex_second_worst = getMatrixValue(response, 1, response->col-1);
        old_simplex_best = getMatrixValue(response, 0, response->col-1);

        /*old_simplex_worst = getMatrixValue(response, 0, response->col-1);
        old_simplex_second_worst = getMatrixValue(response, 1, response->col-1);
        old_simplex_best = getMatrixValue(response, response->row-1, response->col-1);*/

        /*set the best response...*/
        for(j = 0; j < response->col; j++){
          setMatrixValue(coefficients, obj, j, getMatrixValue(response, 0, j));
        }

        /* Generating the new response... as x_new = c + c + x_worst */
        for(k = 1; k < response->row; k++){
          for(j = 0; j < new_response->size; j++){
            setDVectorValue(new_response, j, getDVectorValue(new_response, j) + getMatrixValue(response, k, j));
          }
        }

        /*New response*/
        for(j = 0; j < new_response->size-1; j++){
          setDVectorValue(new_response, j, getDVectorValue(new_response, j)/(double)(response->row-1));
          /* xnew = c + c - xlow */
          setDVectorValue(new_response, j, getDVectorValue(new_response, j)+getDVectorValue(new_response, j)-getMatrixValue(response, 0, j));
        }

        /* puts("New Response"); */
        for(j = 0; j < response->col-1; j++){
          setMatrixValue(response, 0, j, getDVectorValue(new_response, j));
          /* printf("%f\t", getMatrixValue(response, low_q2_id, j)); */
        }

        /*puts("All Response");
        PrintMatrix(response);

        puts("Best Response");
        PrintDVector(coefficients);
        */

        loop = 0;
        while(loop < 3){
          i = 0;

          if(FLOAT_EQ(getMatrixValue(coefficients, obj, coefficients->col-1), getMatrixValue(response, 0, response->col-1), PLSCONVERGENCE)){
            loop++;
          }
          else{
            loop = 0;
            for(j = 0; j < response->col; j++){
              setMatrixValue(coefficients, obj, j, getMatrixValue(response, 0, j));
            }
          }

          for(j = 0; j < mx->col; j++){
            if(operator_ == 0){ /*sum*/
              setMatrixValue(predmx, 0, j, getMatrixValue(mx, obj, j) + getMatrixValue(response, response->row-1, j));
            }
            else{
              setMatrixValue(predmx, 0, j, getMatrixValue(mx, obj, j) * getMatrixValue(response, response->row-1, j));
            }
          }

          for(j = 0; j < my->col; j++){
            setMatrixValue(predmy, 0, j, getMatrixValue(my, obj, j));
          }
          /*
          puts("New Matrix");
          PrintMatrix(predmx);
          PrintMatrix(predmy);
          */
          initMatrix(&predictedy);
          initMatrix(&predscore);
          PLSScorePredictor(predmx, subm, npc, &predscore);
          PLSYPredictor(predscore, subm, mod_cutoff, &predictedy);

          new_yerror = 0.f;
          n = d = 0.f;
          for(j = 0; j < my->col; j++){
            n += sqrt(square(getMatrixValue(predmy, 0, j) - getMatrixValue(predictedy, 0, j)));
            d += getMatrixValue(predmy, 0, j);
          }
          n /= (double)my->col;
          d /= (double)my->col;
          new_yerror = n / d;
          res = new_yerror;

          setMatrixValue(response, response->row-1, response->col-1, new_yerror); /* new_yerror is the new response */
          DelMatrix(&predscore);
          DelMatrix(&predictedy);
          /*
          puts("Response Calculated");
          PrintMatrix(response);
          */

          /* sort matrix by response column
          * the best response is the latest row
          */
          MatrixSort(response, response->col-1);

          worst = getMatrixValue(response, response->row-1, response->col-1);
          second_worst = getMatrixValue(response, 1, response->col-1);
          best = getMatrixValue(response, 0, response->col-1);

        /* 0 is the worst model... so delete and replace with the new model with new coordinate as simplex modified:
          *
          * Now calculating the centroids
          */

          for(j = 0; j < new_response->size; j++){
            setDVectorValue(new_response, j, coeff);
          }

          for(k = 1; k < response->row; k++){
            for(j = 0; j < new_response->size; j++){
              setDVectorValue(new_response, j, getDVectorValue(new_response, j) + getMatrixValue(response, k, j));
            }
          }

          /*Create the new response*/

          if(res < old_simplex_best){
            /*
            * x_new = c + a(c - x_worst)
            *
            * a = 2
            */

            for(j = 0; j < new_response->size; j++){
              setDVectorValue(new_response, j, getDVectorValue(new_response, j)/(double)(response->row-1));
              /* xnew = c + c - x_low for the standard simplex*/
              setDVectorValue(new_response, j, getDVectorValue(new_response, j) + 2*(getDVectorValue(new_response, j)-getMatrixValue(response, response->row-1, j)));
            }
          }
          else if(res < old_simplex_worst && res > old_simplex_second_worst){
            /*
            * x_new = c + b(c - x_worst)
            *
            * b = 0.5
            */

            for(j = 0; j < new_response->size; j++){
              setDVectorValue(new_response, j, getDVectorValue(new_response, j)/(double)(response->row-1));
              /* xnew = c + c - x_low for the standard simplex*/
              setDVectorValue(new_response, j, getDVectorValue(new_response, j) + 0.5*(getDVectorValue(new_response, j)-getMatrixValue(response, response->row-1, j)));
            }
          }
          else if(res > old_simplex_worst){
            /*
            * x_new = c - b(c - x_worst)
            *
            * b = 0.5
            */
            for(j = 0; j < new_response->size; j++){
              setDVectorValue(new_response, j, getDVectorValue(new_response, j)/(double)(response->row-1));
              /* xnew = c + c - x_low for the standard simplex*/
              setDVectorValue(new_response, j, getDVectorValue(new_response, j) - 0.5*(getDVectorValue(new_response, j)-getMatrixValue(response, response->row-1, j)));
            }
          }
          else{
            /*
            * x_new = c + c - x_worst
            */

            for(j = 0; j < new_response->size; j++){
              setDVectorValue(new_response, j, getDVectorValue(new_response, j)/(double)(response->row-1));
              /* xnew = c + c - x_low for the standard simplex*/
              setDVectorValue(new_response, j, getDVectorValue(new_response, j)+getDVectorValue(new_response, j)-getMatrixValue(response, response->row-1, j));
            }
          }

          old_simplex_best = best;
          old_simplex_worst = worst;
          old_simplex_second_worst = second_worst;

          /*puts("New Response");*/
          for(j = 0; j < response->col-1; j++){
            setMatrixValue(response, response->row-1, j, getDVectorValue(new_response, j));
            /* printf("%f\t", getMatrixValue(response, 0, j)); */
          }

          /*
          puts("\nAll response with the new..");
          PrintMatrix(response);
          puts("----------------\n");

          printf("%f %f %f\n", yerror, new_yerror, error);
          */
        }
      }
      else{ /* Skip this object because the error in y prediction is < 0.5*/
        continue;
      }
    }

    DelPLSModel(&subm);
    DelDVector(&colaverage);
  }


  ResizeMatrix(&coefficients_, mx->row, mx->col);

  for(i = 0; i < mx->row; i++){
    for(j = 0; j < mx->col; j++){
      setMatrixValue(coefficients_, i, j, getMatrixValue(coefficients, i, j));
    }
  }

  DelMatrix(&coefficients);
  DelMatrix(&mod_sdep);
  DelMatrix(&mod_pred_y);
  DelMatrix(&submx);
  DelMatrix(&submy);
  DelMatrix(&predmx);
  DelMatrix(&predmy);
  DelDVector(&new_response);
  DelMatrix(&response);
}

void SINGLESIMPLS_OBJMOD(matrix *mx, matrix *my, size_t npc, size_t xautoscaling, size_t yautoscaling,
                           size_t vtype, size_t ngroup, size_t niter,
                           uivector *var_on, uivector *modobjects, double resolution, size_t operator_, matrix *sumcoeff, ssignal *s)
{
  size_t obj, i, j, k, loop, mod_cutoff;
  double n, d, new_yerror, step, coeff, res, worst, second_worst, best, old_simplex_worst, old_simplex_second_worst, old_simplex_best;
  dvector *new_response, *colaverage, *sdeprowaverage, *tmp;
  matrix *coefficients, *submx, *submy, *subpredmx, *subpredmy, *predmx, *predscore, *predmy, *predictedy, *response;

  PLSMODEL *subm;

  /*var_on->siye == mx->col*/
  NewMatrix(&response, var_on->size+1, var_on->size+1); /* +1 because the last is the response value */

  NewDVector(&new_response, var_on->size); /*row vector to store the best response*/

  NewMatrix(&predmx, 1, mx->col); /* 1 because we have only an object to make a good prediction */
  NewMatrix(&predmy, 1, my->col);

  /* Search outlier and do simplex for each outlier... */
  initMatrix(&submx);
  initMatrix(&submy);
  initMatrix(&subpredmx);
  initMatrix(&subpredmy);

  for(obj = 0; obj < modobjects->size; obj++){
    if(getUIVectorValue(modobjects, obj) == 1){ /*testset*/
      tmp = getMatrixRow(mx, obj);
      MatrixAppendRow(&subpredmx, tmp);
      DelDVector(&tmp);
      tmp = getMatrixRow(my, obj);
      MatrixAppendRow(&subpredmy, tmp);
      DelDVector(&tmp);
    }
    else{ /*dataset*/
      tmp = getMatrixRow(mx, obj);
      MatrixAppendRow(&submx, tmp);
      DelDVector(&tmp);
      tmp = getMatrixRow(my, obj);
      MatrixAppendRow(&submy, tmp);
      DelDVector(&tmp);
    }
  }

  NewMatrix(&coefficients, subpredmx->row, var_on->size+1); /*row vector to store the best response*/

  /*
  puts("Matrix Training model");
  PrintMatrix(submx);
  PrintMatrix(submy);

  PrintMatrix(subpredmx);
  PrintMatrix(subpredmy);

  */
  if(submx->row > 0){
    /*build  submodel in order to do the training  for the other outliers... */
    initDVector(&colaverage);
    MatrixColAverage(submx, &colaverage);

    NewPLSModel(&subm);
    PLS(submx, submy, npc, xautoscaling, yautoscaling, subm, s);

    if(vtype == 0){
      PLSLOOCV(submx, submy, xautoscaling, yautoscaling, npc, NULL, &subm->q2y, &subm->sdep, &subm->predicted_y, NULL, s);
    }
    else{
      PLSRandomGroupsCV(submx, submy, xautoscaling, yautoscaling, npc, ngroup, niter, NULL, &subm->q2y, &subm->sdep, &subm->predicted_y, NULL, s);
    }

    /*get the PC with a low error: mod_cutoff

    if only one y
      MatrixGetMinValue(sdep, &mod_cutoff, NULL);
    else
      get the medium error value
    */

    initDVector(&sdeprowaverage);
    MatrixRowAverage(subm->sdep, &sdeprowaverage); /*medium of sdep if more than one y*/
    mod_cutoff = 0;
    for(i = 1; i < sdeprowaverage->size; i++){
      if(getDVectorValue(sdeprowaverage, i) < getDVectorValue(sdeprowaverage, mod_cutoff)){
        mod_cutoff = i;
      }
      else{
        continue;
      }
    }
    DelDVector(&sdeprowaverage);

    for(obj = 0; obj < subpredmx->row; obj++){
      /* printf("#### OPTIMIZE OBJECT No.: %lu ####\n", obj);*/
      /*Perform k+1 experiments on the vertice of a simplex
      * Build the Simplex matrix experiment with the response value as last column of response
      *
      *  coeff = 0.f if we sum in new_mx
      *  coeff = 1.f if we moltiplicate in new_mx
      */

      if(operator_ == 0){ /*sum*/
        coeff = 0.f;
      }
      else{
        coeff = 1.f;
      }

      for(k = 0; k < response->row; k++){
        if(k == 0){
          for(j = 0; j < response->col-1; j++){/*-1 because the last column is the score*/
            if(1 == getUIVectorValue(var_on, j)){
              setMatrixValue(response, k, j, coeff);
            }
            else{
              continue;
            }
          }
        }
        else{
          for(j = 0; j < response->col-1; j++){
            if(1 == getUIVectorValue(var_on, j)){
              step = (getDVectorValue(colaverage, j)*resolution) / 100.f;
              if(j < k-1){
                setMatrixValue(response, k, j, (step)/2.);
              }
              else if(j == k-1){
                setMatrixValue(response, k, j, step);
              }
              else{
                continue;
              }
            }
            else{
              setMatrixValue(response, k, j, coeff);
            }
          }
        }
      }

      /* CALCULATE THE FIRST CICLE OF RESPONSE AND GET THE LOW ERROR IN PREDICTION*/
      for(k = 0; k < response->row; k++){
        /* build an prediction matrix in predmx, predmy with the simplex optimization... */
        for(j = 0; j < subpredmx->col; j++){
          if(operator_ == 0){ /*sum*/
            setMatrixValue(predmx, 0, j, getMatrixValue(subpredmx, obj, j) + getMatrixValue(response, k, j));
          }
          else{
            setMatrixValue(predmx, 0, j, getMatrixValue(subpredmx, obj, j) * getMatrixValue(response, k, j));
          }
        }

        for(j = 0; j < subpredmy->col; j++){
          setMatrixValue(predmy, 0, j, getMatrixValue(subpredmy, obj, j));
        }


        /*
        puts("New Matrix");
        PrintMatrix(predmx);
        PrintMatrix(predmy);
        */
        initMatrix(&predictedy);
        initMatrix(&predscore);
        PLSScorePredictor(predmx, subm, npc, &predscore);
        PLSYPredictor(predscore, subm, mod_cutoff, &predictedy);

        new_yerror = 0.f;
        n = d = 0.f;
        for(j = 0; j < subpredmy->col; j++){
          n += sqrt(square(getMatrixValue(predmy, 0, j) - getMatrixValue(predictedy, 0, j)));
          d += getMatrixValue(predmy, 0, j);
        }
        n /= (double)subpredmy->col;
        d /= (double)subpredmy->col;
        new_yerror = n / d;

        setMatrixValue(response, k, response->col-1, new_yerror); /* yerror is the new response */
        DelMatrix(&predscore);
        DelMatrix(&predictedy);
      }


      /*
      puts("Simplex response matrix coefficients");
      PrintMatrix(response);
      usleep(100000);
      */
      MatrixSort(response, response->col-1); /*the worst response is at row col-1 and the position were wil be builded the ner  response*/
      old_simplex_worst = getMatrixValue(response, response->row-1, response->col-1);
      old_simplex_second_worst = getMatrixValue(response, 1, response->col-1);
      old_simplex_best = getMatrixValue(response, 0, response->col-1);

      /*old_simplex_worst = getMatrixValue(response, 0, response->col-1);
      old_simplex_second_worst = getMatrixValue(response, 1, response->col-1);
      old_simplex_best = getMatrixValue(response, response->row-1, response->col-1);*/

      /*set the best response...*/
      for(j = 0; j < response->col; j++){
        setMatrixValue(coefficients, obj, j, getMatrixValue(response, 0, j));
      }

      /* Generating the new response... as x_new = c + c + x_worst */
      for(k = 1; k < response->row; k++){
        for(j = 0; j < new_response->size; j++){
          setDVectorValue(new_response, j, getDVectorValue(new_response, j) + getMatrixValue(response, k, j));
        }
      }

      /*New response*/
      for(j = 0; j < new_response->size-1; j++){
        setDVectorValue(new_response, j, getDVectorValue(new_response, j)/(double)(response->row-1));
        /* xnew = c + c - xlow */
        setDVectorValue(new_response, j, getDVectorValue(new_response, j)+getDVectorValue(new_response, j)-getMatrixValue(response, 0, j));
      }

      /* puts("New Response"); */
      for(j = 0; j < response->col-1; j++){
        setMatrixValue(response, 0, j, getDVectorValue(new_response, j));
        /* printf("%f\t", getMatrixValue(response, low_q2_id, j)); */
      }

      loop = 0;
      while(loop < 3){
        i = 0;
        if(FLOAT_EQ(getMatrixValue(coefficients, obj, coefficients->col-1), getMatrixValue(response, 0, response->col-1), PLSCONVERGENCE)){
          loop++;
        }
        else{
          loop = 0;
          for(j = 0; j < response->col; j++){
            setMatrixValue(coefficients, obj, j, getMatrixValue(response, 0, j));
          }
        }

        for(j = 0; j < subpredmx->col; j++){
          if(operator_ == 0){ /*sum*/
            setMatrixValue(predmx, 0, j, getMatrixValue(subpredmx, obj, j) + getMatrixValue(response, response->row-1, j));
          }
          else{
            setMatrixValue(predmx, 0, j, getMatrixValue(subpredmx, obj, j) * getMatrixValue(response, response->row-1, j));
          }
        }

        for(j = 0; j < subpredmy->col; j++){
          setMatrixValue(predmy, 0, j, getMatrixValue(subpredmy, obj, j));
        }
        /*
        puts("New Matrix");
        PrintMatrix(predmx);
        PrintMatrix(predmy);
        */
        initMatrix(&predictedy);
        initMatrix(&predscore);
        PLSScorePredictor(predmx, subm, npc, &predscore);
        PLSYPredictor(predscore, subm, mod_cutoff, &predictedy);

        new_yerror = 0.f;
        n = d = 0.f;
        for(j = 0; j < subpredmy->col; j++){
          n += sqrt(square(getMatrixValue(predmy, 0, j) - getMatrixValue(predictedy, 0, j)));
          d += getMatrixValue(predmy, 0, j);
        }
        n /= (double)subpredmy->col;
        d /= (double)subpredmy->col;
        new_yerror = n / d;
        res = new_yerror;

        setMatrixValue(response, response->row-1, response->col-1, new_yerror); /* new_yerror is the new response */
        DelMatrix(&predscore);
        DelMatrix(&predictedy);
        /*
        puts("Response Calculated");
        PrintMatrix(response);
        */

        /* sort matrix by response column
        * the best response is the latest row
        */
        MatrixSort(response, response->col-1);

        worst = getMatrixValue(response, response->row-1, response->col-1);
        second_worst = getMatrixValue(response, 1, response->col-1);
        best = getMatrixValue(response, 0, response->col-1);

      /* 0 is the worst model... so delete and replace with the new model with new coordinate as simplex modified:
        *
        * Now calculating the centroids
        */

        for(j = 0; j < new_response->size; j++){
          setDVectorValue(new_response, j, coeff);
        }

        for(k = 1; k < response->row; k++){
          for(j = 0; j < new_response->size; j++){
            setDVectorValue(new_response, j, getDVectorValue(new_response, j) + getMatrixValue(response, k, j));
          }
        }

        /*Create the new response*/
        if(res < old_simplex_best){
          /*
          * x_new = c + a(c - x_worst)
          *
          * a = 2
          */

          for(j = 0; j < new_response->size; j++){
            setDVectorValue(new_response, j, getDVectorValue(new_response, j)/(double)(response->row-1));
            /* xnew = c + c - x_low for the standard simplex*/
            setDVectorValue(new_response, j, getDVectorValue(new_response, j) + 2*(getDVectorValue(new_response, j)-getMatrixValue(response, response->row-1, j)));
          }
        }
        else if(res < old_simplex_worst && res > old_simplex_second_worst){
          /*
          * x_new = c + b(c - x_worst)
          *
          * b = 0.5
          */

          for(j = 0; j < new_response->size; j++){
            setDVectorValue(new_response, j, getDVectorValue(new_response, j)/(double)(response->row-1));
            /* xnew = c + c - x_low for the standard simplex*/
            setDVectorValue(new_response, j, getDVectorValue(new_response, j) + 0.5*(getDVectorValue(new_response, j)-getMatrixValue(response, response->row-1, j)));
          }
        }
        else if(res > old_simplex_worst){
          /*
          * x_new = c - b(c - x_worst)
          *
          * b = 0.5
          */
          for(j = 0; j < new_response->size; j++){
            setDVectorValue(new_response, j, getDVectorValue(new_response, j)/(double)(response->row-1));
            /* xnew = c + c - x_low for the standard simplex*/
            setDVectorValue(new_response, j, getDVectorValue(new_response, j) - 0.5*(getDVectorValue(new_response, j)-getMatrixValue(response, response->row-1, j)));
          }
        }
        else{
          /*
          * x_new = c + c - x_worst
          */

          for(j = 0; j < new_response->size; j++){
            setDVectorValue(new_response, j, getDVectorValue(new_response, j)/(double)(response->row-1));
            /* xnew = c + c - x_low for the standard simplex*/
            setDVectorValue(new_response, j, getDVectorValue(new_response, j)+getDVectorValue(new_response, j)-getMatrixValue(response, response->row-1, j));
          }
        }

        old_simplex_best = best;
        old_simplex_worst = worst;
        old_simplex_second_worst = second_worst;

        /*puts("New Response");*/
        for(j = 0; j < response->col-1; j++){
          setMatrixValue(response, response->row-1, j, getDVectorValue(new_response, j));
          /* printf("%f\t", getMatrixValue(response, 0, j)); */
        }

        /*
        puts("\nAll response with the new..");
        PrintMatrix(response);
        puts("----------------\n");
        usleep(300000);

        printf("%f %f %f\n", yerror, new_yerror, error);
        */
      }

    }

    DelPLSModel(&subm);
    DelDVector(&colaverage);
  }

  ResizeMatrix(&sumcoeff, mx->row, mx->col);

  k = 0;
  for(i = 0; i < mx->row; i++){
    if(getUIVectorValue(modobjects, i) == 1){
      for(j = 0; j < mx->col; j++){
        setMatrixValue(sumcoeff, i, j, getMatrixValue(coefficients, k, j));
      }
      k++;
    }
    else{
      continue;
    }
  }

  DelMatrix(&subpredmx);
  DelMatrix(&subpredmy);
  DelMatrix(&coefficients);
  DelMatrix(&submx);
  DelMatrix(&submy);
  DelMatrix(&predmx);
  DelMatrix(&predmy);
  DelDVector(&new_response);
  DelMatrix(&response);
}
