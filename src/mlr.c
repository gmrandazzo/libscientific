/* mlr.c
*
* Copyright (C) <2016>  Giuseppe Marco Randazzo
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "mlr.h"
#include "memwrapper.h"
#include "matrix.h"

#include "algebra.h"
#include "numeric.h"

#include <math.h>


void NewMLRModel(MLRMODEL** m)
{
  (*m) = xmalloc(sizeof(MLRMODEL));
  initMatrix(&(*m)->b);
  initMatrix(&(*m)->recalculated_y);
  initMatrix(&(*m)->predicted_y);
  initMatrix(&(*m)->recalc_residuals);
  initMatrix(&(*m)->pred_residuals);
  initDVector(&(*m)->ymean);
  initDVector(&(*m)->r2y_model);
  initDVector(&(*m)->q2y);
  initDVector(&(*m)->sdep);/* Standard Deviation over Prediction */
  initDVector(&(*m)->sdec); /* Standard Deviation over Recalculating */
  initDVector(&(*m)->bias);
  initMatrix(&(*m)->r2q2scrambling);
}

void DelMLRModel(MLRMODEL** m)
{
  DelMatrix(&(*m)->b);
  DelMatrix(&(*m)->recalculated_y);
  DelMatrix(&(*m)->predicted_y);
  DelMatrix(&(*m)->recalc_residuals);
  DelMatrix(&(*m)->pred_residuals);
  DelDVector(&(*m)->ymean);
  DelDVector(&(*m)->r2y_model);
  DelDVector(&(*m)->q2y);
  DelDVector(&(*m)->sdep);/* Standard Deviation over Prediction */
  DelDVector(&(*m)->sdec); /* Standard Deviation over Recalculating */
  DelDVector(&(*m)->bias);
  DelMatrix(&(*m)->r2q2scrambling);
  xfree((*m));
}

void MLR(matrix* mx, matrix* my, MLRMODEL* model, ssignal *s)
{
  size_t i, j;
  matrix *mx_;
  dvector *y_;
  dvector *b_;
  NewMatrix(&mx_, mx->row, mx->col+1);
  NewDVector(&y_, mx->row);

  /* Preparate the matrix to correlate with y trough OLS */
  MatrixCheck(mx);

  /*Adding constant coefcient for polinomial equation */
  for(j = 1; j < mx_->col; j++){
    for(i = 0; i < mx_->row; i++){
      setMatrixValue(mx_, i, j, getMatrixValue(mx, i, j-1));
    }
  }

  for(i = 0; i < mx_->row; i++){
    setMatrixValue(mx_, i, 0, 1);
  }

  /*build models with more than one y*/
  for(j = 0; j < my->col; j++){
    for(i = 0; i < my->row; i++){
      setDVectorValue(y_, i, getMatrixValue(my, i, j));
    }

    initDVector(&b_);

    OrdinaryLeastSquares(mx_, y_, b_);

    MatrixAppendCol(&model->b, b_);
    DelDVector(&b_);
  }

  /* Recalculate y from model */
  MatrixColAverage(my, &model->ymean);
  MLRPredictY(mx, my, model, &model->recalculated_y, &model->recalc_residuals, &model->r2y_model, &model->sdec);

  DelDVector(&y_);
  DelMatrix(&mx_);
}

/* Predict y using the equation:
 *
 * y = c + ax1 + bx2 + cx3 + dx4 + ...
 *
 */
void MLRPredictY(matrix* mx, matrix *my, MLRMODEL* model, matrix** predicted_y, matrix** predicted_residuals, dvector** r2y, dvector** sdep)
{
  if(mx->col == (model->b->row-1)){
    size_t i, j, k;
    double ypred, rss, tss;

    if((*predicted_y)->row != mx->row && (*predicted_y)->col != model->b->col){
      ResizeMatrix(predicted_y, mx->row, model->b->col);
    }

    for(k = 0; k < model->b->col; k++){ /*FOR EACH Y*/
      for(i = 0; i < mx->row; i++){ /* FOR EACH OBJECT */
        ypred = getMatrixValue(model->b, 0, k); /* CONSTANT c + ...*/
        for(j = 0; j < mx->col; j++){
          ypred += getMatrixValue(mx, i, j)*getMatrixValue(model->b, j+1, k);
        }
        setMatrixValue((*predicted_y), i, k, ypred);
      }
    }

    if(my != 0){
      if(predicted_residuals != 0){
        ResizeMatrix(predicted_residuals, my->row, my->col);
      }
      /*Calculate Residual Error
      *      ei = y^i - yi.
      * where:
      *   yi = real y
      *   y^i = predicted y
      *
      * Calculate R^2 and SDEC
      */

      for(j = 0; j < my->col; j++){
        rss = tss = 0.f;
        for(i = 0; i < my->row; i++){
          if(predicted_residuals != 0){
            setMatrixValue((*predicted_residuals), i, j, getMatrixValue((*predicted_y), i, j) - getMatrixValue(my, i, j));
          }
          rss += square(getMatrixValue(my, i, j) - getMatrixValue((*predicted_y), i, j));
          tss += square(getMatrixValue(my, i, j) - getDVectorValue(model->ymean, j));
        }

        if(r2y != 0){
          DVectorAppend(r2y, 1-(rss/tss));
        }

        if(sdep != 0){
          DVectorAppend(sdep, sqrt(rss/my->row));
        }
      }
    }
  }
  else{
    fprintf(stderr, "Error!! Unable to compute MLR Prediction.\n The number of variable differ with the variables in model.\n");
  }
}


void MLRRegressionStatistics(matrix *my_true, matrix *my_pred, dvector** ccoeff, dvector **stdev, dvector **bias)
{
  size_t i, j;
  dvector *ymean;

  /*Calculate the Q2 and SDEP */
  if(ccoeff != NULL)
    DVectorResize(ccoeff,my_true->col);

  if(stdev != NULL)
    DVectorResize(stdev, my_true->col);

  if(bias != NULL)
    DVectorResize(bias, my_true->col);

  initDVector(&ymean);
  MatrixColAverage(my_true, &ymean);

  for(j = 0; j < my_true->col; j++){
    double ssreg = 0.f;
    double sstot = 0.f;
    for(i = 0; i < my_true->row; i++){
      ssreg += square(my_pred->data[i][j] - my_true->data[i][j]);
      sstot += square(my_true->data[i][j] - ymean->data[j]);
    }

    if(bias != NULL){
      double sum_yi = 0.f, sum_xi = 0.f;
      /*ypredaverage /= (double)my->row;*/
      for(i = 0; i < my_true->row; i++){
        sum_yi+=(my_pred->data[i][j]*(my_true->data[i][j]-ymean->data[j]));
        sum_xi+=(my_true->data[i][j]*(my_true->data[i][j]-ymean->data[j]));
      }
      /*sum_yi/sum_xi = m */
      (*bias)->data[j] = fabs(1 - sum_yi/sum_xi);
      /*k = Y-(X*b);*/
    }

    if(ccoeff != NULL)
      (*ccoeff)->data[j] = 1.f - (ssreg/sstot);

    if(stdev != NULL)
      (*stdev)->data[j] = sqrt(ssreg/(double)my_true->row);
  }

  DelDVector(&ymean);
}

void MLRYScrambling(matrix *mx, matrix *my,
                        size_t block, size_t valtype, size_t rgcv_group, size_t rgcv_iterations,
                        matrix **r2q2scrambling, ssignal *s)
{
  size_t scrambiterations, iterations_, i, j, k, n, y_, blocksize;
  int id;
  double temp;
  matrix *randomY, *sorted_y_id, *sorty, *gid;
  dvector *tmpq2, *yaverage;
  MLRMODEL *tmpmod;

  srand(mx->row*mx->col*my->col*block);
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
  blocksize = (size_t)ceil(mx->row/(double)block);
  blocksize += (size_t)ceil((float)((blocksize*block) - mx->row)/  block);

  NewMatrix(&gid, block, blocksize);
  MatrixSet(gid, -2);
  /* Crate the boxes to fill -2 means no value to fill, -1 means value to fill*/
  for(i = 0, j = 0, k = 0; i < mx->row; i++){
    if(j < block){
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


  /*Create a the r2q2scrambling matrix */
  ResizeMatrix(r2q2scrambling, 1+scrambiterations, my->col*3); /* each row is a model. first my->col columns are the r2 scrambled/real, second are the r2 scrabled/scrambled and third the q2 scrabled/scrabmled */

  /*First row is the model not scrambled...*/
  initDVector(&tmpq2);

  if(valtype == 0){
    MLRLOOCV(mx, my, &tmpq2, NULL, NULL, NULL, NULL, s);
  }
  else{
    MLRRandomGroupsCV(mx, my, rgcv_group, rgcv_iterations, &tmpq2, NULL, NULL, NULL, NULL, s);
  }

  NewMLRModel(&tmpmod);
  MLR(mx, my, tmpmod, s);

  /* Calculate y real vs yscrambled and add other r2 q2 */
  initDVector(&yaverage);
  MatrixColAverage(my, &yaverage);
  for(j = 0; j < my->col; j++){
    double rss = 0.f, tss = 0.f;
    for(i = 0; i < my->row; i++){
      rss += square(my->data[i][j] - my->data[i][j]);
      tss += square(my->data[i][j] - yaverage->data[j]);
    }
    (*r2q2scrambling)->data[0][j] = 1 - (rss/tss);
    (*r2q2scrambling)->data[0][j+my->col] = tmpmod->r2y_model->data[j];
    (*r2q2scrambling)->data[0][j+my->col+my->col] = tmpq2->data[j];
  }

  DelMLRModel(&tmpmod);
  DelDVector(&tmpq2);

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

        /*
        puts("MX");
        PrintMatrix(mx);
        puts("MY");
        PrintMatrix(my);
        puts("RANDOMY");
        PrintMatrix(randomY);
        */
        /* Calculate calculate Q2 for y predicted...*/

        initDVector(&tmpq2);

        if(valtype == 0){
          MLRLOOCV(mx, randomY, &tmpq2, NULL, NULL, NULL, NULL, s);
        }
        else{
          MLRRandomGroupsCV(mx, randomY, rgcv_group, rgcv_iterations, &tmpq2, NULL, NULL, NULL, NULL, s);
        }

        NewMLRModel(&tmpmod);
        MLR(mx, randomY, tmpmod, s);

        /* Calculate y real vs yscrambled and add other r2 q2 */
        for(j = 0; j < my->col; j++){
          double rss = 0.f, tss = 0.f;
          for(i = 0; i < my->row; i++){
            rss += square(my->data[i][j] - randomY->data[i][j]);
            tss += square(my->data[i][j] - yaverage->data[j]);
          }
          (*r2q2scrambling)->data[iterations_+1][j] = 1 - (rss/tss);
          (*r2q2scrambling)->data[iterations_+1][j+my->col] = tmpmod->r2y_model->data[j];
          (*r2q2scrambling)->data[iterations_+1][j+my->col+my->col] = tmpq2->data[j];
        }

        DelMLRModel(&tmpmod);
        DelDVector(&tmpq2);
        iterations_++;
      }
    }
  }

  DelDVector(&yaverage);
  DelMatrix(&sorted_y_id);
  DelMatrix(&randomY);
  DelMatrix(&gid);
}

/*cv_ are the cross validated coefficients used to plot predicted vs experimental. if is null is not calculated.*/

void MLRRandomGroupsCV(matrix *mx, matrix *my,
                        size_t group, size_t iterations,
                        dvector **q2y, dvector **sdep, dvector **bias, matrix **predicted_y, matrix** pred_residuals, ssignal *s)
{
  if(mx->row == my->row && group > 0 && iterations > 0){
    size_t iterations_, i, j, k, n, g;
    matrix *gid; /* randomization and storing id for each random group into a matrix */

    /* Matrix for compute the PLS models for groups */
    matrix *subX;
    matrix *subY;
    MLRMODEL *subm;


    /* matrix for the randomg group to predict */
    matrix *predictX;
    matrix *realY;
    matrix *predictY;
    dvector *mean_col_y;

    dvector *rss;
    dvector *tss;
    uivector *predictcounter; /*count how many times the y are predicted out of the model*/

    matrix *_predicted_y_;


    NewMatrix(&gid, group, (size_t)ceil(mx->row/(double)group));

    NewDVector(&rss, my->col);
    NewDVector(&tss, my->col);

    if(predicted_y != NULL){
      ResizeMatrix(predicted_y, my->row, my->col);
      _predicted_y_ = (*predicted_y);
    }
    else{
      NewMatrix(&_predicted_y_, my->row, my->col);
    }

    if(pred_residuals != NULL)
      ResizeMatrix(pred_residuals, my->row, my->col);

    NewUIVector(&predictcounter, my->row);

    srand(group*mx->row*iterations);

    iterations_ = 0;
    while(iterations_ <  iterations){
      if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
        break;
      }
      else{
        /* Divide in group  all the Dataset */
        MatrixSet(gid, -1);

        /* step 1 generate the random groups */
        k = 0;
        for(i = 0; i <  gid->row; i++){
          for(j = 0; j <  gid->col; j++){
            do{
              n = (size_t)rand() % (mx->row);
            } while(ValInMatrix(gid, n) == 1 && k < (mx->row));
            if(k < mx->row){
              setMatrixValue(gid, i, j, n);
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

        /*step 2*/
        for(g = 0; g < gid->row; g++){ /*For aeach group */
          /* Estimate how many objects are inside the sub model without the group "g" */
          n = 0;
          for(i = 0; i < gid->row; i++){
            if(i != g){
              for(j = 0; j < gid->col; j++){
                if((int)getMatrixValue(gid, i, j) != -1)
                  n++;
                else
                  continue;
              }
            }
            else
              continue;
          }

          /*Allocate the submodel*/
          NewMatrix(&subX, n, mx->col);
          NewMatrix(&subY, n, my->col);

          /* Estimate how many objects are inside the group "g" to predict*/
          n = 0;
          for(j = 0; j < gid->col; j++){
            if((int)getMatrixValue(gid, g, j) != -1)
              n++;
            else
              continue;
          }


          /*Allocate the */
          NewMatrix(&predictX, n, mx->col);
          NewMatrix(&realY, n, my->col);


          /* copy the submodel values */

          for(i = 0, k = 0; i < gid->row; i++){
            if(i != g){
              for(j = 0; j < gid->col; j++){
                size_t a =  (size_t)getMatrixValue(gid, i, j); /* get the row index */
                if(a != -1){
                  for(n = 0; n < mx->col; n++){
                    setMatrixValue(subX, k, n, getMatrixValue(mx, a, n));
                  }
                  for(n = 0; n < my->col; n++){
                    setMatrixValue(subY, k, n, getMatrixValue(my, a, n));
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
            size_t a = (size_t)getMatrixValue(gid, g, j);
            if(a != -1){
              for(n = 0; n < mx->col; n++){
                setMatrixValue(predictX, k, n, getMatrixValue(mx, a, n));
              }
              for(n = 0; n < my->col; n++){
                setMatrixValue(realY, k, n, getMatrixValue(my, a, n));
              }
              k++;
            }
            else{
              continue;
            }
          }

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

          initDVector(&mean_col_y);
          MatrixColAverage(subY, &mean_col_y);

          NewMLRModel(&subm);

          MLR(subX, subY, subm, s);

          initMatrix(&predictY);
          MLRPredictY(predictX, NULL, subm, &predictY, NULL, NULL, NULL);

          for(j = 0; j < realY->col; j++){
            for(i = 0; i < predictX->row; i++){
              setDVectorValue(rss, j, getDVectorValue(rss, j) + square(getMatrixValue(realY, i, j) - getMatrixValue(predictY, i, j)));
              setDVectorValue(tss, j, getDVectorValue(tss, j) + square(getMatrixValue(realY, i, j) - getDVectorValue(mean_col_y, j)));
            }
          }


          for(j = 0, k = 0; j < gid->col; j++){
            size_t a = (size_t)getMatrixValue(gid, g, j); /*riga dell'oggetto....*/
            if(a != -1){
              setUIVectorValue(predictcounter, a, getUIVectorValue(predictcounter, a)+1);
              for(n = 0; n < predictY->col; n++){
                _predicted_y_->data[a][n] += predictY->data[k][n];
                /*setMatrixValue((*predicted_y), a, n, (getMatrixValue((*predicted_y), a, n) + getMatrixValue(predictY, k, n)));*/
              }
              k++;
            }
            else{
              continue;
            }
          }

          DelDVector(&mean_col_y);
          DelMLRModel(&subm);
          DelMatrix(&predictY);
          DelMatrix(&realY);
          DelMatrix(&predictX);
          DelMatrix(&subY);
          DelMatrix(&subX);
        }
        iterations_++;
      }
    }

    /*Finalize the output by dividing for the number of iterations*/
    for(j = 0; j < my->col; j++){
      DVectorAppend(q2y, 1 - (getDVectorValue(rss, j)/getDVectorValue(tss, j)));
    }

    if(sdep != NULL)
      for(j = 0; j < my->col; j++)
        DVectorAppend(sdep, sqrt((getDVectorValue(rss, j)/(iterations))/my->row));


    for(i = 0; i < _predicted_y_->row; i++){
      for(j = 0; j < _predicted_y_->col; j++){
        _predicted_y_->data[i][j] /= predictcounter->data[i];
        if(pred_residuals != NULL)
          (*pred_residuals)->data[i][j] = _predicted_y_->data[i][j] - my->data[i][j];
      }
    }

    /*for(i = 0; i < (*predicted_y)->row; i++){
      for(j = 0; j < (*predicted_y)->col; j++){
        setMatrixValue((*predicted_y), i, j, getMatrixValue((*predicted_y), i, j)/getUIVectorValue(predictcounter, i));
        setMatrixValue((*pred_residuals), i, j,  getMatrixValue((*predicted_y), i, j) - getMatrixValue(my, i, j));
      }
    }*/


    if(bias != NULL){
      dvector *ymean;
      initDVector(&ymean);
      MatrixColAverage(my, &ymean);

      for(j = 0; j < my->col; j++){
        double sum_xi = 0.f, sum_yi = 0.f;
        /*ypredaverage /= (double)my->row;*/
        for(size_t i = 0; i < my->row; i++){
          sum_yi+=(_predicted_y_->data[i][j]*(my->data[i][j]-ymean->data[j]));
          sum_xi+=(my->data[i][j]*(my->data[i][j]-ymean->data[j]));
        }

        DVectorAppend(bias, fabs(1 - sum_yi/sum_xi));
        /*k = Y-(X*b);*/
      }
      DelDVector(&ymean);
    }

    if(predicted_y == NULL)
      DelMatrix(&_predicted_y_);
    DelUIVector(&predictcounter);

    DelMatrix(&gid);
    DelDVector(&rss);
    DelDVector(&tss);
  }
  else{
    fprintf(stderr, "Error!! Unable to compute MLR Cross Validation!!\n");
  }
}

/*loov_ are the loo validated coefficients used to plot predicted vs experimental. if is null is not calculated.*/
void MLRLOOCV(matrix *mx, matrix *my,
                        dvector **q2y, dvector **sdep, dvector **bias, matrix **predicted_y, matrix **pred_residuals, ssignal *s)
{
 if(mx->row == my->row){
    size_t j, k, l, model;

    /* Matrix for compute the PLS models for groups */
    matrix *subX;
    matrix *subY;
    MLRMODEL *subm;

    /* matrix for the randomg group to predict */
    matrix *predictX;
    matrix *realY;
    matrix *predictY;
    dvector *mean_col_y;

    matrix *predictedy;

    dvector *rss;
    dvector *tss;

    NewMatrix(&subX, mx->row-1, mx->col);
    NewMatrix(&subY, my->row-1, my->col);

    NewMatrix(&predictX, 1, mx->col);
    NewMatrix(&realY, 1, my->col);
    NewMatrix(&predictY, 1, my->col);

    NewDVector(&rss, my->col);
    NewDVector(&tss, my->col);

    if(predicted_y != NULL){
      ResizeMatrix(predicted_y, my->row, my->col);
    }
    else{
      NewMatrix(&predictedy, my->row, my->col);
    }

    if(pred_residuals!= NULL){
      ResizeMatrix(pred_residuals, my->row, my->col);
    }


    for(model = 0; model < mx->row; model++){ /* we compute mx->row models  */
      if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
        break;
      }
      else{
        /* copy data into subX, subY, predictX and realY */
        l = 0;
        for(j = 0; j < mx->row; j++){
          if(j != model){
            for(k = 0; k < mx->col; k++){
              setMatrixValue(subX, l, k, getMatrixValue(mx, j, k));
            }
            for(k = 0; k < my->col; k++){
              setMatrixValue(subY, l, k, getMatrixValue(my, j, k));
            }
            l++;
          }
          else{
            for(k = 0; k < mx->col; k++){
              setMatrixValue(predictX, 0, k, getMatrixValue(mx, j, k));
            }
            for(k = 0; k < my->col; k++){
              setMatrixValue(realY, 0, k, getMatrixValue(my, j, k));
            }
          }
        }


        /*Compute the MLRModel*/
        initDVector(&mean_col_y);
        MatrixColAverage(subY, &mean_col_y);

        NewMLRModel(&subm);
        MLR(subX, subY, subm, s);
        MLRPredictY(predictX, NULL, subm, &predictY, NULL, NULL, NULL);

        for(j = 0; j < realY->col; j++){
            setDVectorValue(rss, j, getDVectorValue(rss, j) + square(getMatrixValue(realY, 0, j) - getMatrixValue(predictY, 0, j)));
            setDVectorValue(tss, j, getDVectorValue(tss, j) + square(getMatrixValue(realY, 0, j) - getDVectorValue(mean_col_y, j)));
        }

        for(j = 0; j < predictY->col; j++){
          if(predicted_y != NULL)
            setMatrixValue((*predicted_y), model, j, getMatrixValue(predictY, 0, j));
          else
            predictedy->data[model][j] = predictY->data[0][j];

          if(pred_residuals != NULL){
            setMatrixValue((*pred_residuals), model, j, getMatrixValue((*predicted_y), model, j) - getMatrixValue(my, model, j));
          }
        }

        DelMLRModel(&subm);
        DelDVector(&mean_col_y);
      }
    }

    /*Finalize the output by dividing for the number of models*/
    for(j = 0; j < my->col; j++){
      DVectorAppend(q2y, 1 - (getDVectorValue(rss, j)/getDVectorValue(tss, j)));
      if(sdep != NULL)
        DVectorAppend(sdep, sqrt(getDVectorValue(rss, j)/mx->row));
    }


    if(bias != NULL){
      dvector *ymean;
      initDVector(&ymean);
      MatrixColAverage(my, &ymean);

      for(j = 0; j < my->col; j++){
        double sum_xi = 0.f, sum_yi = 0.f;
        /*ypredaverage /= (double)my->row;*/
        for(size_t i = 0; i < my->row; i++){
          if(predicted_y != NULL)
            sum_yi+=((*predicted_y)->data[i][j]*(my->data[i][j]-ymean->data[j]));
          else
            sum_yi+=(predictedy->data[i][j]*(my->data[i][j]-ymean->data[j]));

          sum_xi+=(my->data[i][j]*(my->data[i][j]-ymean->data[j]));
        }

        DVectorAppend(bias, fabs(1 - sum_yi/sum_xi));
        /*k = Y-(X*b);*/
      }
      DelDVector(&ymean);
    }

    if(predicted_y == NULL)
      DelMatrix(&predictedy);

    DelDVector(&rss);
    DelDVector(&tss);
    DelMatrix(&realY);
    DelMatrix(&predictX);
    DelMatrix(&predictY);
    DelMatrix(&subX);
    DelMatrix(&subY);
  }
  else{
    fprintf(stderr, "Error!! Unable to compute MLR Leave One Out Validation!!\n");
  }
}


void PrintMLR(MLRMODEL *m)
{
  puts("b Coeffcicient");
  PrintMatrix(m->b);

  puts("R^2");
  PrintDVector(m->r2y_model);

  puts("SDEC");
  PrintDVector(m->sdec);

  puts("Recalculated y");
  PrintMatrix(m->recalculated_y);

  puts("Recalculated Residuals");
  PrintMatrix(m->recalc_residuals);

  if(m->q2y->size > 0){
    puts("Q^2");
    PrintDVector(m->q2y);

    puts("SDEP");
    PrintDVector(m->sdep);

    puts("BIAS");
    PrintDVector(m->bias);

    puts("Predicted y");
    PrintMatrix(m->predicted_y);

    puts("Predicted Residuals");
    PrintMatrix(m->pred_residuals);
  }

  if(m->r2q2scrambling->row > 0){
    puts("R2 Q^2 Y Scrambling");
    PrintMatrix(m->r2q2scrambling);
  }
}
