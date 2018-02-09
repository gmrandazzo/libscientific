/* lda.c
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

#include <stdio.h>
#include <math.h>
#include "lda.h"
#include "memwrapper.h"
#include "numeric.h"
#include "tensor.h"
#include <pthread.h>

/* extern int rand_r(unsigned int *seedp);*/

void NewLDAModel(LDAMODEL** m)
{
  (*m) = xmalloc(sizeof(LDAMODEL));
  initUIVector(&((*m)->classid));
  initMatrix(&((*m)->evect));
  initTensor(&((*m)->mnpdf));
  initTensor(&((*m)->features));
  initMatrix(&((*m)->inv_cov));
  initMatrix(&((*m)->mu));
  initMatrix(&((*m)->fsdev));
  initMatrix(&((*m)->fmean));
  initDVector(&((*m)->eval));
  initDVector(&((*m)->pprob));
  initDVector(&((*m)->sens));
  initDVector(&((*m)->spec));
  initDVector(&((*m)->ppv)); /* o precision */
  initDVector(&((*m)->npv));
  initDVector(&((*m)->acc));
  (*m)->nclass = (*m)->class_start = 0;
}

void DelLDAModel(LDAMODEL** m)
{
  DelDVector(&((*m)->acc));
  DelDVector(&((*m)->sens));
  DelDVector(&((*m)->spec));
  DelDVector(&((*m)->ppv));
  DelDVector(&((*m)->npv));
  DelDVector(&((*m)->pprob));
  DelDVector(&((*m)->eval));
  DelMatrix(&((*m)->mu));
  DelMatrix(&((*m)->fsdev));
  DelMatrix(&((*m)->fmean));
  DelMatrix(&((*m)->inv_cov));
  DelMatrix(&((*m)->evect));
  DelUIVector(&((*m)->classid));
  DelTensor(&((*m)->features));
  DelTensor(&((*m)->mnpdf));
  xfree((*m));
}

void PrintLDAModel(LDAMODEL* m)
{
  puts("Eigenvalues"); PrintDVector(m->eval);
  puts("Eigenvectors"); PrintMatrix(m->evect);
  puts("Features");
  PrintTensor(m->features);
  puts("Multivariate Normal Distribution");
  PrintTensor(m->mnpdf);

  puts("Class Average. Each row represent a class");
  PrintMatrix(m->mu);

  puts("Validation...");
  puts("Senstivity"); PrintDVector(m->sens);
  puts("Specificity"); PrintDVector(m->spec);
  puts("Positive Predicted Value"); PrintDVector(m->ppv);
  puts("Negative Predicted Value"); PrintDVector(m->npv);
  puts("Accuracy"); PrintDVector(m->acc);
}


/* last column must be an integer column which define the classes */
void LDA(matrix *mx, uivector *y, LDAMODEL *lda)
{
  size_t i, j, l, k, cc, imin, imax;
  tensor *classes;
  tensor *S;
  matrix *X, *X_T, *Sb, *Sw, *InvSw_Sb;
  dvector *mutot;
  dvector *classmu;
  dvector *evect_, *ldfeature;
  matrix *covmx;


  lda->nclass = 0;


  imin = imax = y->data[0];

  for(i = 1; i < y->size; i++){
    if(y->data[i] > imax){
      imax = y->data[i];
    }

    if(y->data[i] < imin){
      imin = y->data[i];
    }
  }

  /* get the number of classes */
  if(imin == 0){
    lda->class_start = 0;
    lda->nclass = imax + 1;
  }
  else{
    lda->class_start = 1;
    lda->nclass = imax;
  }

  /*printf("nclass %d\n", (int)lda->nclass);*/

  /* Copy data */
  NewMatrix(&X, mx->row, mx->col);
  MatrixCopy(mx, &X);
  MatrixCheck(X);
  /*
  for(j = 0; j < mx->col-1; j++){
    for(i = 0; i < mx->row; i++){
      X->data[i][j] = mx->data[i][j];
    }
  }
  */

  /*create classes of objects */
  NewTensor(&classes, lda->nclass);
  UIVectorResize(&lda->classid, mx->row);
  j = 0;
  if(imin == 0){
    for(k = 0; k < lda->nclass; k++){
      cc = 0;
      for(i = 0; i < X->row; i++){
        if(y->data[i] == k)
          cc++;
        else
          continue;
      }
      NewTensorMatrix(&classes, k, cc, X->col);

      cc = 0;
      for(i = 0; i < X->row; i++){
        if(y->data[i] == k){
          for(l = 0; l < X->col; l++){
            classes->m[k]->data[cc][l] = X->data[i][l];
          }
          lda->classid->data[j] = i;
          j++;
          cc++;
        }
        else{
          continue;
        }
      }
    }
  }
  else{
    for(k = 0; k < lda->nclass; k++){
      cc = 0;
      for(i = 0; i < X->row; i++){
        if(y->data[i] == k+1)
          cc++;
        else
          continue;
      }
      NewTensorMatrix(&classes, k, cc, X->col);

      cc = 0;
      for(i = 0; i < X->row; i++){
        if(y->data[i] == k+1){
          for(l = 0; l < X->col; l++){
            classes->m[k]->data[cc][l] = X->data[i][l];
          }
          lda->classid->data[j] = i;
          j++;
          cc++;
        }
        else{
          continue;
        }
      }
    }
  }

  /*puts("Classes"); PrintTensor(classes);*/

  /* Compute the prior probability */
  for(k = 0; k < classes->order; k++){
    DVectorAppend(&lda->pprob, (classes->m[k]->row/(double)X->row));
  }

  /*
  puts("Prior Probability");
  PrintDVector(lda->pprob);
  */

  /*Compute the mean of each class*/
  for(k = 0; k < classes->order; k++){
    initDVector(&classmu);
    MatrixColAverage(classes->m[k], &classmu);

    MatrixAppendRow(&lda->mu, classmu);

    DelDVector(&classmu);
  }

  /*puts("Class Mu");
  FindNan(lda->mu);
  PrintMatrix(lda->mu);*/

  /*Calculate the total mean of samples..*/
  initDVector(&mutot);
  MatrixColAverage(X, &mutot);

  /*puts("Mu tot");
  PrintDVector(mutot);*/

  /*NewDVector(&mutot, mu->col);

  for(k = 0; k < mu->row; k++){
    for(i = 0; i < mu->col; i++){
      mutot->data[i] += mu->data[k][i];
    }
  }

  for(i = 0; i < mutot->size; i++){
    mutot->data[i] /= nclasses;
  }*/


  /*Centering data before computing the scatter matrix*/
  for(k = 0; k < classes->order; k++){
    for(i = 0; i < classes->m[k]->row; i++){
      for(j = 0; j < classes->m[k]->col; j++){
        classes->m[k]->data[i][j] -= mutot->data[j];
      }
    }
  }

  /*
  puts("Classes");
  for(i = 0; i < classes->order; i++){
    FindNan(classes->m[i]);
  }
   PrintTensor(classes);
  */
  /*Compute the scatter matrix
   * S = nobj - 1 * covmx
   */
  initTensor(&S);
  NewMatrix(&covmx, X->col, X->col);
  for(k = 0; k < classes->order; k++){
    matrix *m_T;
    NewMatrix(&m_T, classes->m[k]->col, classes->m[k]->row);
    MatrixTranspose(classes->m[k], m_T);

    MatrixDotProduct(m_T, classes->m[k], covmx);

    for(i = 0; i < covmx->row; i++){
      for(j = 0; j < covmx->col; j++){
        covmx->data[i][j] /= classes->m[k]->row;
      }
    }
    TensorAppendMatrix(&S, covmx);
    MatrixSet(covmx, 0.f);
    DelMatrix(&m_T);
  }
  /*
  puts("Scatter Matrix");
  for(i = 0; i < S->order; i++)
    FindNan(S->m[i]);

  PrintTensor(S);*/

  /* Compute the class scatter which represent the covariance matrix */
  NewMatrix(&Sw, X->col, X->col);

  for(k = 0; k < S->order; k++){
    for(i = 0; i  < S->m[k]->row; i++){
      for(j = 0; j  < S->m[k]->col; j++){
        Sw->data[i][j] += (double)(classes->m[k]->row/(double)X->row) *  S->m[k]->data[i][j];
      }
    }
  }
  /*
  puts("Class scatter matrix Sw");
  FindNan(Sw);
  PrintMatrix(Sw);
  */
  /*Compute the between class scatter matrix Sb*/
  NewMatrix(&Sb, X->col, X->col);
  for(k = 0; k < lda->mu->row; k++){ /*for each class of object*/
    cc = classes->m[k]->row;
    for(i = 0; i < Sb->row; i++){
      for(j = 0; j < Sb->col; j++){
        Sb->data[i][j] += cc * (lda->mu->data[k][i] - mutot->data[i]) * (lda->mu->data[k][j] - mutot->data[j]);
      }
    }
  }

  /*
  puts("Between class scatter matrix Sb");
  FindNan(Sb);
  PrintMatrix(Sb); */

  /* Computing the LDA projection */
  /*puts("Compute Matrix Inversion");*/
  MatrixInversion(Sw, &lda->inv_cov);

  double ss = 0.f;
  for(i = 0; i < lda->inv_cov->row; i++){
    for(j = 0; j < lda->inv_cov->col; j++){
      ss += square(lda->inv_cov->data[i][j]);
    }
    if(_isnan_(ss))
      break;
  }

  if(FLOAT_EQ(ss, 0.f, EPSILON) || _isnan_(ss)){
    /*Do SVD as pseudoinversion accordin to Liu et al. because matrix is nonsingular
     *
     * JUN LIU et al, Int. J. Patt. Recogn. Artif. Intell. 21, 1265 (2007). DOI: 10.1142/S0218001407005946
     * EFFICIENT PSEUDOINVERSE LINEAR DISCRIMINANT ANALYSIS AND ITS NONLINEAR FORM FOR FACE RECOGNITION
     *
     *
     * Sw`^-1 = Q * G^-1 * Q_T
     * Q G AND Q_T come from SVD
     */
    MatrixPseudoinversion(Sw, &lda->inv_cov);

    /*
    NewMatrix(&A_T, Sw->col, Sw->row);
    MatrixInversion(Sw, &A_T);
    NewMatrix(&A_T_Sw, A_T->row, Sw->col);
    MatrixDotProduct(A_T, Sw, A_T_Sw);

    initMatrix(&A_T_Sw_inv);
    MatrixInversion(A_T_Sw, &A_T_Sw_inv);

    MatrixDotProduct(A_T_Sw_inv, A_T, lda->inv_cov);
    DelMatrix(&A_T);
    DelMatrix(&A_T_Sw);
    DelMatrix(&A_T_Sw_inv);
    */
  }

  /*puts("Inverted Covariance Matrix from Sw");
   FindNan(lda->inv_cov);
   PrintMatrix(lda->inv_cov);
  */
  /*puts("Compute Matrix Dot Product");*/
  NewMatrix(&InvSw_Sb, lda->inv_cov->row, Sb->col);
  MatrixDotProduct(lda->inv_cov, Sb, InvSw_Sb);

  /*puts("InvSw_Sb"); PrintMatrix(InvSw_Sb);*/
  /*puts("Compute Eigenvectors");*/
  EVectEval(InvSw_Sb, &lda->eval, &lda->evect);
  /*EvectEval3(InvSw_Sb, InvSw_Sb->row, &lda->eval, &lda->evect);*/
  /*EVectEval(InvSw_Sb, &lda->eval, &lda->evect); */

  /* Calculate the new projection in the feature space
   *
   * and the multivariate normal distribution
   */
/*       Remove centering data   */
  for(k = 0; k < classes->order; k++){
    for(i = 0; i < classes->m[k]->row; i++){
      for(j = 0; j < classes->m[k]->col; j++){
        classes->m[k]->data[i][j] += mutot->data[j];
      }
    }
  }

  initMatrix(&X_T);

  for(k = 0; k < classes->order; k++){
    /*printf("row %d  col %d\n", (int)classes->m[k]->row, (int)classes->m[k]->col);*/
    AddTensorMatrix(&lda->features, classes->m[k]->row, classes->m[k]->col);
    AddTensorMatrix(&lda->mnpdf, classes->m[k]->row, classes->m[k]->col);
  }

  NewDVector(&evect_, lda->evect->row);
  initDVector(&ldfeature);

  ResizeMatrix(&lda->fmean, classes->order, lda->evect->col);
  ResizeMatrix(&lda->fsdev, classes->order, lda->evect->col);

  for(l = 0; l < lda->evect->col; l++){

    for(i = 0; i < lda->evect->row; i++){
      evect_->data[i] = lda->evect->data[i][l];
    }

    for(k = 0; k < classes->order; k++){

      ResizeMatrix(&X_T, classes->m[k]->col, classes->m[k]->row);
      MatrixTranspose(classes->m[k], X_T);
      DVectorResize(&ldfeature, classes->m[k]->row);

      DVectorMatrixDotProduct(X_T, evect_, ldfeature);

      for(i = 0; i < ldfeature->size; i++){
        lda->features->m[k]->data[i][l] = ldfeature->data[i];
      }

/*        Calculate the multivariate normal distribution  */
      double mean = 0.f, sdev = 0.f;
      DVectorMean(ldfeature, &mean);
      DVectorSDEV(ldfeature, &sdev);

      lda->fmean->data[k][l] = mean;
      lda->fsdev->data[k][l] = sdev;
      for(i = 0; i < ldfeature->size; i++){
        lda->mnpdf->m[k]->data[i][l] = 1./sqrt(2 * _pi_* sdev) * exp(-square((ldfeature->data[i] - mean)/sdev)/2.f);
      }
    }
  }

  DelDVector(&evect_);
  DelMatrix(&covmx);
  DelDVector(&ldfeature);
  DelDVector(&mutot);
  DelMatrix(&Sb);
  DelMatrix(&InvSw_Sb);
  DelTensor(&classes);
  DelTensor(&S);
  DelMatrix(&Sw);
  DelMatrix(&X_T);
  DelMatrix(&X);
}

void LDAPrediction(matrix *mx, LDAMODEL *lda, matrix **pfeatures, matrix **probability, matrix **mnpdf, uivector **classprediction)
{
  /* probability function:
   * f = mu_class * inv_cov * mx.T  - 1/2. * mu_class * inv_cov * mu_class.T + ln(prior_probability_class)
   */
  size_t i, j, argmax, id;
  matrix *mx_T;
  dvector *x, *mu, *C_x_T, *C_mu_T, *ldfeature, *evect_;

  NewMatrix(&mx_T, mx->col, mx->row);
  MatrixTranspose(mx, mx_T);

  NewDVector(&C_x_T, lda->inv_cov->row);
  NewDVector(&C_mu_T, lda->inv_cov->row);

  ResizeMatrix(probability, mx->row, lda->nclass);

  for(i = 0; i < mx_T->col; i++){
    x = getMatrixColumn(mx_T, i); /* object k*/
    for(j = 0; j < lda->nclass; j++){
      mu = getMatrixRow(lda->mu, j); /*each row correspond to a class of objects*/
      MatrixDVectorDotProduct(lda->inv_cov, x, C_x_T);
      MatrixDVectorDotProduct(lda->inv_cov, mu, C_mu_T);
      (*probability)->data[i][j] = DVectorDVectorDotProd(mu, C_x_T) - (0.5 * DVectorDVectorDotProd(mu, C_mu_T)) + log(lda->pprob->data[j]);
      DVectorSet(C_x_T, 0.f);
      DVectorSet(C_mu_T, 0.f);
      DelDVector(&mu);
    }
    DelDVector(&x);
  }

  for(i = 0; i < (*probability)->row; i++){
    argmax = 0;
    for(j = 1; j < (*probability)->col; j++){
      if((*probability)->data[i][j] > (*probability)->data[i][argmax]){
        argmax = j;
      }
      else{
        continue;
      }
    }

    if(lda->class_start  == 0){
      UIVectorAppend(classprediction, argmax);
    }
    else{
      UIVectorAppend(classprediction, argmax+1);
    }
  }

  /* Predict the the new projection in the feature space */
  ResizeMatrix(mnpdf, mx->row, lda->evect->col);
  NewDVector(&ldfeature, mx->row);
  for(i = 0; i < lda->evect->col; i++){
    evect_ = getMatrixColumn(lda->evect, i);
    DVectorMatrixDotProduct(mx_T, evect_, ldfeature);
    DelDVector(&evect_);
    MatrixAppendCol(pfeatures, ldfeature);

    for(j = 0; j < ldfeature->size; j++){
      id = (*classprediction)->data[j];
      if(lda->class_start == 1)
        id -= 1;

      (*mnpdf)->data[j][i] = 1./sqrt(2 * _pi_* lda->fsdev->data[id][i]) * exp(-square((ldfeature->data[j] - lda->fmean->data[id][i])/lda->fsdev->data[id][i])/2.f);
    }

    DVectorSet(ldfeature, 0.f);
  }

  DelDVector(&ldfeature);
  DelMatrix(&mx_T);
  DelDVector(&C_x_T);
  DelDVector(&C_mu_T);
}

void LDAError(matrix *mx, uivector *my, LDAMODEL *lda, dvector **sens, dvector **spec, dvector **ppv, dvector **npv, dvector **acc)
{
  size_t  i, k, pos;
  matrix *pfeatures;
  matrix *probability;
  matrix *mnpdf;
  uivector *classprediction;
  uivector* tp, *fp, *fn, *tn;
  double d, p, n;

  initMatrix(&pfeatures);
  initMatrix(&probability);
  initMatrix(&mnpdf);
  initUIVector(&classprediction);
  LDAPrediction(mx, lda, &pfeatures, &probability, &mnpdf, &classprediction);

  NewUIVector(&tp, lda->nclass);
  NewUIVector(&fp, lda->nclass);
  NewUIVector(&tn, lda->nclass);
  NewUIVector(&fn, lda->nclass);

  if(lda->class_start == 1){
    pos = -1;
  }
  else{
    pos = 0;
  }

  for(k = 0; k < lda->nclass; k++){ /*for each class*/
    for(i = 0; i < my->size; i++){
      /*printf("class %d  object class %d  predicted class %d\n", (int)k, (int)(my->data[i] - pos), (int)(classprediction->data[i] - pos));*/
      if(k == (my->data[i] - pos)){ /* oggetto della stessa classe */
        if(my->data[i] == classprediction->data[i]){
          tp->data[k]++;
          /*puts("TP");*/
        }
        else{
          fn->data[k]++;
          /*puts("FN");*/
        }
      }
      else{ /*oggetto di diversa classe dal quella cercata*/
        if(classprediction->data[i] - pos == k){ /*false positive!*/
          fp->data[k]++;
          /*puts("FP");*/
        }
        else{
          tn->data[k]++;
          /*puts("TN");*/
        }
      }
    }
  }

  if((*sens)->size != lda->nclass)
    DVectorResize(sens, lda->nclass);

  if((*spec)->size != lda->nclass)
    DVectorResize(spec, lda->nclass);

  if((*ppv)->size != lda->nclass)
    DVectorResize(ppv, lda->nclass);

  if((*npv)->size != lda->nclass)
    DVectorResize(npv, lda->nclass);

  if((*acc)->size != lda->nclass)
    DVectorResize(acc, lda->nclass);

  for(i = 0; i < tp->size; i++){ /* for each class*/
    /* sensitivity calculation */
    d = tp->data[i] + fn->data[i];
    if(FLOAT_EQ(d, 0, 1e-4)){
      (*sens)->data[i] = 0;
    }
    else{
      (*sens)->data[i] = tp->data[i] / d;
    }

    /* specificity calculation */
    d = tn->data[i] + fp->data[i];
    if(FLOAT_EQ(d, 0, 1e-4)){
      (*spec)->data[i] = 0;
    }
    else{
      (*spec)->data[i] = tn->data[i] / d;
    }

    /* ppv calculation */
    d = tp->data[i] + fp->data[i];
    if(FLOAT_EQ(d, 0, 1e-4)){
      (*ppv)->data[i] = 0;
    }
    else{
      (*ppv)->data[i] = tp->data[i] / d;
    }

    /* npv calculation */
    d = fn->data[i] + tn->data[i];
    if(FLOAT_EQ(d, 0, 1e-4)){
      (*npv)->data[i] = 0;
    }
    else{
      (*npv)->data[i] = tn->data[i] / d;
    }

    /* acc calculation */
    p = tp->data[i] + fn->data[i];
    n = tn->data[i] + fp->data[i];
    (*acc)->data[i] = (tp->data[i] + tn->data[i]) / (double)(p+n);
  }

  DelUIVector(&tp);
  DelUIVector(&fp);
  DelUIVector(&fn);
  DelUIVector(&tn);
  DelMatrix(&mnpdf);
  DelMatrix(&pfeatures);
  DelMatrix(&probability);
  DelUIVector(&classprediction);
}

void LDAStatistics(dvector *y_true, dvector *y_score, matrix **roc, double *roc_auc, matrix **precision_recal, double *pr_auc)
{
  ROC(y_true, y_score,  roc, roc_auc);
  PrecisionRecall(y_true, y_score,  precision_recal, pr_auc);
}


typedef struct{
  matrix *mx;
  uivector *my;
  dvector *sens, *spec, *ppv, *npv, *acc;  /*OUPUT*/
  size_t group; /*INPUT*/
  unsigned int srand_init;
} lda_rgcv_th_arg;


void *LDARandomGroupCVModel(void *arg_)
{
  size_t i, j, k, n, g;
  lda_rgcv_th_arg *arg;
  matrix *gid; /* randomization and storing id for each random group into a matrix */

  /* Matrix for compute the PLS models for groups */
  matrix *subX;
  uivector *subY;
  LDAMODEL *subm;

  /*matrix to predict*/
  matrix *predictX;
  uivector *realY;


  dvector *sens, *spec, *ppv, *npv, *acc;

  initDVector(&sens);
  initDVector(&spec);
  initDVector(&ppv);
  initDVector(&npv);
  initDVector(&acc);

  arg = (lda_rgcv_th_arg*) arg_;

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
    PrintTensor(subX);
    puts("Y:");
    PrintTensor(subY);

    puts("\n\nPredict Group\nX:");
    PrintTensor(predictX);
    puts("RealY:");
    PrintTensor(realY);
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
    NewUIVector(&subY, n);

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
    NewUIVector(&realY, n);

    /* copy the submodel values */

    for(i = 0, k = 0; i < gid->row; i++){
      if(i != g){
        for(j = 0; j < gid->col; j++){
          size_t a =  (size_t)gid->data[i][j]; /* get the row index */
          if(a != -1){
            for(n = 0; n < arg->mx->col; n++){
              /*setMatrixValue(subX, k, n, getMatrixValue(arg->mx, a, n));*/
              subX->data[k][n] = arg->mx->data[a][n];
            }
            subY->data[k] = arg->my->data[a];
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
        realY->data[k] = arg->my->data[a];
        k++;
      }
      else{
        continue;
      }
    }

    NewLDAModel(&subm);

    LDA(subX, subY, subm);

    LDAError(predictX, realY, subm, &sens, &spec, &ppv, &npv, &acc);

    if(arg->sens->size == 0){
      DVectorCopy(sens, &arg->sens);
    }
    else{
      for(i = 0; i < sens->size; i++){
        arg->sens->data[i] += sens->data[i];
      }
    }

    if(arg->spec->size == 0){
      DVectorCopy(spec, &arg->spec);
    }
    else{
      for(i = 0; i < spec->size; i++){
        arg->spec->data[i] += spec->data[i];
      }
    }

    if(arg->ppv->size == 0){
      DVectorCopy(ppv, &arg->ppv);
    }
    else{
      for(i = 0; i < ppv->size; i++){
        arg->ppv->data[i] += ppv->data[i];
      }
    }

    if(arg->npv->size == 0){
      DVectorCopy(npv, &arg->npv);
    }
    else{
      for(i = 0; i < npv->size; i++){
        arg->npv->data[i] += npv->data[i];
      }
    }

    if(arg->acc->size == 0){
      DVectorCopy(acc, &arg->acc);
    }
    else{
      for(i = 0; i < acc->size; i++){
        arg->acc->data[i] += acc->data[i];
      }
    }

    DelLDAModel(&subm);
    DelMatrix(&subX);
    DelUIVector(&subY);
    DelMatrix(&predictX);
    DelUIVector(&realY);
  }

  DelDVector(&acc);
  DelDVector(&sens);
  DelDVector(&spec);
  DelDVector(&ppv);
  DelDVector(&npv);
  DelMatrix(&gid);
  return 0;
}

void LDARandomGroupsCV(matrix *mx, uivector *my, size_t group, size_t iterations, dvector **sens, dvector **spec, dvector **ppv, dvector **npv, dvector **acc, size_t nthreads, ssignal *s)
{
  if(mx->row == my->size && group > 0 && iterations > 0){
    size_t th, iterations_, i;
    pthread_t *threads;

    /* each thread have its argument type */
    threads = xmalloc(sizeof(pthread_t)*nthreads);
    lda_rgcv_th_arg *arg;
    arg = xmalloc(sizeof(lda_rgcv_th_arg)*nthreads);


    for(th = 0; th < nthreads; th++){
      initDVector(&arg[th].sens);
      initDVector(&arg[th].spec);
      initDVector(&arg[th].ppv);
      initDVector(&arg[th].npv);
      initDVector(&arg[th].acc);
    }

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
          arg[th].srand_init = (unsigned int) group + mx->row + iterations + th + iterations_;

          /*Reset everithings*/
          for(i = 0; i < arg[th].sens->size; i++){
            arg[th].sens->data[i] = arg[th].spec->data[i] = arg[th].ppv->data[i] = arg[th].npv->data[i] = 0.f;
          }

          pthread_create(&threads[th], NULL, LDARandomGroupCVModel, (void*) &arg[th]);
        }

        /* Wait till threads are complete before main continues. Unless we  */
        /* wait we run the risk of executing an exit which will terminate   */
        /* the process and all threads before the threads have completed.   */
        for(th = 0; th < nthreads; th++){
          pthread_join(threads[th], NULL);
        }

        /* finalize thread outputs and free the memory.....*/
        for(th = 0; th < nthreads; th++){
          if((*sens)->size == 0){
            DVectorCopy(arg[th].sens, sens);
          }
          else{
            for(i = 0; i < arg[th].sens->size; i++){
              (*sens)->data[i] += arg[th].sens->data[i];
            }
          }

          if((*spec)->size == 0){
            DVectorCopy(arg[th].spec, spec);
          }
          else{
            for(i = 0; i < arg[th].spec->size; i++){
              (*spec)->data[i] += arg[th].spec->data[i];
            }
          }

          if((*ppv)->size == 0){
            DVectorCopy(arg[th].ppv, ppv);
          }
          else{
            for(i = 0; i < arg[th].ppv->size; i++){
              (*ppv)->data[i] += arg[th].ppv->data[i];
            }
          }

          if((*npv)->size == 0){
            DVectorCopy(arg[th].npv, npv);
          }
          else{
            for(i = 0; i < arg[th].npv->size; i++){
              (*npv)->data[i] += arg[th].npv->data[i];
            }
          }

          if((*acc)->size == 0){
            DVectorCopy(arg[th].acc, acc);
          }
          else{
            for(i = 0; i < arg[th].acc->size; i++){
              (*acc)->data[i] += arg[th].acc->data[i];
            }
          }
        }
      }
    }

    /*Finalize the output by dividing for the number of iterations*/

    double d = group*iterations;
    for(i = 0; i < (*sens)->size; i++){
      (*sens)->data[i] /= d;
      (*spec)->data[i] /= d;
      (*ppv)->data[i] /= d;
      (*npv)->data[i] /= d;
      (*acc)->data[i] /= d;
    }

    for(th = 0; th < nthreads; th++){
      DelDVector(&arg[th].sens);
      DelDVector(&arg[th].spec);
      DelDVector(&arg[th].ppv);
      DelDVector(&arg[th].npv);
      DelDVector(&arg[th].acc);
    }
    xfree(threads);
    xfree(arg);
  }
  else{
    fprintf(stderr, "Error!! Unable to compute LDA Cross Validation!!\n");
  }
}


/* Leave One Ouut Validation
 *
 * 1) remove one object
 * 2) calculate the model
 * 3) predict the removed object and so the r2 and q2
 */

typedef struct{
  matrix *submx, *predmx;
  uivector *submy, *predmy;
  dvector *sens, *spec, *ppv, *npv, *acc;
} lda_loocv_th_arg;

void *LDALOOModel(void *arg_)
{
  lda_loocv_th_arg *arg;
  arg = (lda_loocv_th_arg*) arg_;

  LDAMODEL *subm;
  NewLDAModel(&subm);

  LDA(arg->submx, arg->submy, subm);

  LDAError(arg->predmx, arg->predmy, subm, &arg->sens, &arg->spec, &arg->ppv, &arg->npv, &arg->acc);

  DelLDAModel(&subm);
  return 0;
}

void LDALOOCV(matrix* mx, uivector* my, dvector** sens, dvector** spec, dvector** ppv, dvector** npv, dvector **acc, size_t nthreads, ssignal *s)
{
 if(mx->row == my->size){
    size_t i, j, k, l, th, model;
    pthread_t *threads;
    lda_loocv_th_arg *arg;

    threads = xmalloc(sizeof(pthread_t)*nthreads);
    arg = xmalloc(sizeof(lda_loocv_th_arg)*nthreads);

    /* initialize threads arguments.. */
    for(th = 0; th < nthreads; th++){
      NewMatrix(&arg[th].submx, mx->row-1, mx->col);
      NewUIVector(&arg[th].submy, my->size-1);
      NewMatrix(&arg[th].predmx, 1, mx->col);
      NewUIVector(&arg[th].predmy, 1);
      initDVector(&arg[th].sens);
      initDVector(&arg[th].spec);
      initDVector(&arg[th].ppv);
      initDVector(&arg[th].npv);
      initDVector(&arg[th].acc);
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
                arg[th].submy->data[l] = my->data[j];
                l++;
              }
              else{
                for(k = 0; k < mx->col; k++){
                  arg[th].predmx->data[0][k] = mx->data[j][k];
                }
                arg[th].predmy->data[0] = my->data[j];
              }
            }

            pthread_create(&threads[th], NULL, LDALOOModel, (void*) &arg[th]);
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
          if((*sens)->size == 0){
            DVectorCopy(arg[th].sens, sens);
          }
          else{
            for(i = 0; i < arg[th].sens->size; i++){
              (*sens)->data[i] += arg[th].sens->data[i];
            }
          }

          if((*spec)->size == 0){
            DVectorCopy(arg[th].spec, spec);
          }
          else{
            for(i = 0; i < arg[th].spec->size; i++){
              (*spec)->data[i] += arg[th].spec->data[i];
            }
          }

          if((*ppv)->size == 0){
            DVectorCopy(arg[th].ppv, ppv);
          }
          else{
            for(i = 0; i < arg[th].ppv->size; i++){
              (*ppv)->data[i] += arg[th].ppv->data[i];
            }
          }

          if((*npv)->size == 0){
            DVectorCopy(arg[th].npv, npv);
          }
          else{
            for(i = 0; i < arg[th].npv->size; i++){
              (*npv)->data[i] += arg[th].npv->data[i];
            }
          }

          if((*acc)->size == 0){
            DVectorCopy(arg[th].acc, acc);
          }
          else{
            for(i = 0; i < arg[th].acc->size; i++){
              (*acc)->data[i] += arg[th].acc->data[i];
            }
          }
        }
      }
    }

    /*Delete thread arguments*/

    for(th = 0; th < nthreads; th++){
      DelMatrix(&arg[th].submx);
      DelUIVector(&arg[th].submy);
      DelMatrix(&arg[th].predmx);
      DelUIVector(&arg[th].predmy);
      DelDVector(&arg[th].sens);
      DelDVector(&arg[th].spec);
      DelDVector(&arg[th].ppv);
      DelDVector(&arg[th].npv);
      DelDVector(&arg[th].acc);
    }

    /*Finalize the output by dividing for the number of models*/
    for(i = 0; i < (*sens)->size; i++){
      (*sens)->data[i] /= mx->row;
      (*spec)->data[i] /= mx->row;
      (*ppv)->data[i] /= mx->row;
      (*npv)->data[i] /= mx->row;
      (*acc)->data[i] /= mx->row;
    }
    xfree(arg);
    xfree(threads);
  }
  else{
    fprintf(stderr, "Error!! Unable to compute LDA Leave One Out Validation!!\n");
  }
}
