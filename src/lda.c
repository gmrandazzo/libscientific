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
#include "metricspace.h"
#include "numeric.h"
#include "tensor.h"
#include "statistic.h"
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
  initMatrix(&((*m)->yscrambling));
  initMatrix(&((*m)->recalculated_y));
  initMatrix(&((*m)->recalculated_residuals));
  initMatrix(&((*m)->predicted_y));
  initMatrix(&((*m)->predicted_residuals));
  initDVector(&((*m)->eval));
  initDVector(&((*m)->pprob));

  initTensor(&((*m)->roc));
  initTensor(&((*m)->pr));
  initDVector(&((*m)->roc_aucs));
  initDVector(&((*m)->pr_aucs));
  (*m)->nclass = (*m)->class_start = 0;
}

void DelLDAModel(LDAMODEL** m)
{
  DelDVector(&((*m)->roc_aucs));
  DelDVector(&((*m)->pr_aucs));
  DelTensor(&((*m)->pr));
  DelTensor(&((*m)->roc));
  DelDVector(&((*m)->pprob));
  DelDVector(&((*m)->eval));
  DelMatrix(&((*m)->recalculated_y));
  DelMatrix(&((*m)->recalculated_residuals));
  DelMatrix(&((*m)->predicted_y));
  DelMatrix(&((*m)->predicted_residuals));
  DelMatrix(&((*m)->yscrambling));
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
  puts("ROC");
  PrintTensor(m->roc);
  PrintDVector(m->roc_aucs);
  puts("Precision-Recall");
  PrintTensor(m->pr);
  PrintDVector(m->pr_aucs);
}


/*
 * Compute LDA model
 * mx is the feature matrix
 * my is a nx1 vector that define the classes
 * lda is the output model
 */
void LDA(matrix *mx, matrix *my, LDAMODEL *lda)
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


  imin = imax = (int)my->data[0][0];

  for(i = 1; i < my->row; i++){
    if((int)my->data[i][0] > imax){
      imax = (int)my->data[i][0];
    }

    if((int)my->data[i][0] < imin){
      imin = (int)my->data[i][0];
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
        if((int)my->data[i][0] == k)
          cc++;
        else
          continue;
      }
      NewTensorMatrix(&classes, k, cc, X->col);

      cc = 0;
      for(i = 0; i < X->row; i++){
        if((int)my->data[i][0] == k){
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
        if((int)my->data[i][0] == k+1)
          cc++;
        else
          continue;
      }
      NewTensorMatrix(&classes, k, cc, X->col);

      cc = 0;
      for(i = 0; i < X->row; i++){
        if((int)my->data[i][0] == k+1){
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

void LDAPrediction(matrix *mx,
                   LDAMODEL *lda,
                   matrix **pfeatures,
                   matrix **probability,
                   matrix **mnpdf,
                   matrix **prediction)
{
  /* probability function:
   * f = mu_class * inv_cov * mx.T  - 1/2. * mu_class * inv_cov * mu_class.T + ln(prior_probability_class)
   */
  size_t i, j, argmax;
  int pos, id;
  matrix *mx_T;
  dvector *x, *mu, *C_x_T, *C_mu_T, *ldfeature, *evect_;

  NewMatrix(&mx_T, mx->col, mx->row);
  MatrixTranspose(mx, mx_T);

  NewDVector(&C_x_T, lda->inv_cov->row);
  NewDVector(&C_mu_T, lda->inv_cov->row);

  ResizeMatrix(probability, mx->row, lda->nclass);
  ResizeMatrix(prediction, mx->row, 1);

  if(lda->class_start == 1){
    pos = -1;
  }
  else{
    pos = 0;
  }

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

    (*prediction)->data[i][0] = (argmax+pos);
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
      id = (int)((*prediction)->data[j][0]+pos);
      (*mnpdf)->data[j][i] = 1./sqrt(2 * _pi_* lda->fsdev->data[id][i]) * exp(-square((ldfeature->data[j] - lda->fmean->data[id][i])/lda->fsdev->data[id][i])/2.f);
    }

    DVectorSet(ldfeature, 0.f);
  }

  DelDVector(&ldfeature);
  DelMatrix(&mx_T);
  DelDVector(&C_x_T);
  DelDVector(&C_mu_T);
}

void LDAError(matrix *mx,
              matrix *my,
              LDAMODEL *lda,
              dvector **sens,
              dvector **spec,
              dvector **ppv,
              dvector **npv,
              dvector **acc)
{
  size_t  i, k;
  int pos;
  matrix *pfeatures;
  matrix *probability;
  matrix *mnpdf;
  matrix *classprediction;
  uivector* tp, *fp, *fn, *tn;
  double d, p, n;

  initMatrix(&pfeatures);
  initMatrix(&probability);
  initMatrix(&mnpdf);
  initMatrix(&classprediction);
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
    for(i = 0; i < my->row; i++){
      /*printf("class %d  object class %d  predicted class %d\n", (int)k, (int)(my->data[i] + pos), (int)(classprediction->data[i] + pos));*/
      if(k == ((int)(my->data[i][0] + pos))){ /* oggetto della stessa classe */
        if((int)my->data[i][0] == (int)classprediction->data[i][0]){
          tp->data[k]++;
          /*puts("TP");*/
        }
        else{
          fn->data[k]++;
          /*puts("FN");*/
        }
      }
      else{ /*oggetto di diversa classe dal quella cercata*/
        if((int)(classprediction->data[i][0] + pos) == k){ /*false positive!*/
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
  DelMatrix(&classprediction);
}


void LDAStatistics(dvector *y_true,
                   dvector *y_pred,
                   matrix **roc,
                   double *roc_auc,
                   matrix **precision_recal,
                   double *pr_auc)
{

  ROC(y_true, y_pred,  roc, roc_auc);
  PrecisionRecall(y_true, y_pred, precision_recal, pr_auc);
}



void LDAMulticlassStatistics(matrix *y_true,
                             matrix *y_pred,
                             tensor **roc,
                             dvector **roc_aucs,
                             tensor **precision_recals,
                             dvector **pr_aucs)
{
  size_t i, j;
  size_t n_classes;
  dvector *ytrue;
  dvector *ypred;
  matrix *roc_;
  matrix *pr_;
  double roc_auc;
  double pr_auc;


  NewDVector(&ytrue, y_true->row);
  NewDVector(&ypred, y_true->row);

  n_classes = getNClasses(y_pred);
  if(n_classes == 2){
    n_classes = 1;
  }

  for(j = 0; j < n_classes; j++){
    for(i = 0; i < y_true->row; i++){
      if((int)y_true->data[i][0] == j){
        ytrue->data[i] = 1;
      }
      else{
        ytrue->data[i] = 0;
      }

      if((int)y_pred->data[i][0] == j){
        ytrue->data[i] = 1;
      }
      else{
        ytrue->data[i] = 0;
      }
    }

    initMatrix(&roc_);
    initMatrix(&pr_);

    ROC(ytrue, ypred,  &roc_, &roc_auc);
    PrecisionRecall(ytrue, ypred, &pr_, &pr_auc);

    DVectorAppend(pr_aucs, pr_auc);
    DVectorAppend(roc_aucs, roc_auc);

    if(roc != NULL){
      TensorAppendMatrix(roc, roc_);
    }

    if(precision_recals != NULL){
      TensorAppendMatrix(precision_recals, pr_);
    }

    DelMatrix(&roc_);
    DelMatrix(&pr_);
  }

  DelDVector(&ytrue);
  DelDVector(&ypred);

}

int getNClasses(matrix *my)
{
  size_t i;
  uivector *nclasses;
  int numclasses;
  initUIVector(&nclasses);
  for(i = 0; i < my->row; i++){
    if(UIVectorHasValue(nclasses, (int)my->data[i][0]) == 0){
      continue;
    }
    else{
        UIVectorAppend(&nclasses, (int)my->data[i][0]);
    }
  }
  numclasses = nclasses->size;
  DelUIVector(&nclasses);
  return numclasses;
}
