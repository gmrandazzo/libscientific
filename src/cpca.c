/* cpca.c
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

#include "cpca.h"
#include "memwrapper.h"
#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "preprocessing.h"
#include "pca.h"
#include "numeric.h"
#include "tensor.h"
#include "vector.h"


void NewCPCAModel(CPCAMODEL **m){
  (*m) = xmalloc(sizeof(CPCAMODEL));
  initTensor(&(*m)->block_scores);
  initTensor(&(*m)->block_loadings);
  initMatrix(&(*m)->super_scores);
  initMatrix(&(*m)->super_weights);
  initDVector(&(*m)->scaling_factor);
  initDVector(&(*m)->total_expvar);
  initDVectorList(&(*m)->block_expvar);
  initDVectorList(&(*m)->colaverage);
  initDVectorList(&(*m)->colscaling);
}


void DelCPCAModel(CPCAMODEL **m)
{
  DelTensor(&(*m)->block_scores);
  DelTensor(&(*m)->block_loadings);
  DelMatrix(&(*m)->super_scores);
  DelMatrix(&(*m)->super_weights);
  DelDVector(&(*m)->scaling_factor);
  DelDVector(&(*m)->total_expvar);
  DelDVectorList(&(*m)->block_expvar);
  DelDVectorList(&(*m)->colaverage);
  DelDVectorList(&(*m)->colscaling);
  xfree((*m));
}

static inline void CalcBlockLoadings(matrix *Xb, dvector *t, dvector *p){
  size_t i;
  double mod_t;
  matrix *Xb_T;

  NewMatrix(&Xb_T, Xb->col, Xb->row);
  MatrixTranspose(Xb, Xb_T);
  MT_MatrixDVectorDotProduct(Xb_T, t, p);
  DelMatrix(&Xb_T);
  mod_t = DVectorDVectorDotProd(t, t);
  for(i = 0; i < p->size; i++)
    p->data[i] /= mod_t;
}

/*
 * Consensus Principal Component Analysis
 *
 * ANALYSIS OF MULTIBLOCK AND HIERARCHICAL PCA AND PLS MODELS
 * JOHAN A. WESTERHUIS, THEODORA KOURTI* AND JOHN F. MACGREGOR
 * J. Chemometrics 12, 301–321 (1998)
 *
 * N.B.: The superscores of CPCA are the scores of a PCA!!
 */
void CPCA(tensor *x, int scaling, size_t npc, CPCAMODEL *model)
{
  size_t i;
  size_t j;
  size_t k;
  size_t pc;
  matrix *T_T;
  matrix *T;
  matrix *Eb_T;
  matrix *Eb_T_E;
  dvector *t_b;
  dvector *t;
  dvector *t_new;
  dvector *p_b;
  dvector *w_T;
  dvector *colvar;
  dvector *tr_orig;
  dvector *local_blockvexp;
  double ss;
  double mod_t;
  tensor *Eb; /* copy of the original data matrix to do autoscing and mean centring */
  size_t row = 0;
  size_t col = 0;

  NewTensor(&Eb, x->order);
  for(k = 0; k < Eb->order; k++){
    NewTensorMatrix(Eb, k, x->m[k]->row, x->m[k]->col);
    row = x->m[k]->row;
    col += x->m[k]->col;
  }

  #ifndef DEBUG
  /*Copy of the matrix to run the PCA on it*/
  matrix *X;
  NewMatrix(&X, row, col);
  #endif

  /* Tensor preprocess */
  TensorPreprocess(x, scaling, model->colaverage, model->colscaling, Eb);

  #ifndef DEBUG
 /* To maintain the same scaling in PCA as is used in CPCA,
  * we apply the following scaling to X for the PCA analysis
  * (step 1)
  */
  int c = 0; //Column
  for(k = 0; k < x->order; k++){
    double m = sqrt((double)x->m[k]->col);
    DVectorAppend(model->scaling_factor, m);
    for(j = 0; j < x->m[k]->col; j++){
      for(i = 0; i < x->m[k]->row; i++){
        X->data[i][c] = Eb->m[k]->data[i][j] / m; /* block scaling N.B.: this is applied to reproduce the pca score */
      }
      c += 1;
    }
  }
  #endif

 /* if the number of principal component selected is major of the permitted
  * set the value to the max value
  */
  for(k = 0; k < Eb->order; k++){
    if(npc > Eb->m[k]->col){
      npc = Eb->m[k]->col;
    }
    else{
      continue;
    }
  }


 /*
  * Calculate the sum of squares of each block
  */
  ss = 0.f;
  for(k = 0; k < x->order; k++){
    for(i = 0; i < Eb->m[k]->row; i++){
      for(j = 0; j < Eb->m[k]->col; j++){
        /* We apply the same scaling factor to mantain the same scaling of CPCA */
        ss += square(Eb->m[k]->data[i][j]/model->scaling_factor->data[k]);
      }
    }
    /*block_sum_of_squares->data[k] = block_ss;*/
  }

 /*
  * Prepare the model to store results
  */
  ResizeMatrix(model->super_scores, X->row, npc);
  ResizeMatrix(model->super_weights, Eb->order, npc);

  for(k = 0; k < Eb->order; k++){
    AddTensorMatrix(model->block_loadings, Eb->m[k]->col, npc);
  }

  #ifdef DEBUG
  PCAMODEL *m;
  NewPCAModel(&m);
  /* X is already centered and scaled PCA will recenter again.... */
  PCA(X, 0, npc, m, NULL);
  #endif

  NewDVector(&tr_orig, Eb->order);
  for(k = 0; k < Eb->order; k++){
    NewMatrix(&Eb_T, Eb->m[k]->col, Eb->m[k]->row);
    MatrixTranspose(Eb->m[k], Eb_T);
    NewMatrix(&Eb_T_E, Eb->m[k]->col, Eb->m[k]->col);
    MatrixDotProduct(Eb_T, Eb->m[k], Eb_T_E);
    tr_orig->data[k] = MatrixTrace(Eb_T_E);
    DelMatrix(&Eb_T_E);
    DelMatrix(&Eb_T);
  }

  NewDVector(&t_b, Eb->m[0]->row);
  NewDVector(&t_new, Eb->m[0]->row);
  NewMatrix(&T, Eb->m[0]->row, Eb->order);
  NewMatrix(&T_T, Eb->order, Eb->m[0]->row);
  NewDVector(&w_T, Eb->order);

  for(pc = 0; pc < npc; pc++){
   /*
    * The super score of CPCA equals the score of PCA, i.e. t_T : t
    */
    double best_colvar = 0.f;
    size_t best_colvar_id = 0;
    size_t best_block_id = 0;
    for(k = 0; k < Eb->order; k++){
      initDVector(&colvar);
      MatrixColVar(Eb->m[k], colvar);
      /* Step 1: select the column vector t with the largest column variance */
      j = 0;
      for(i = 1; i < Eb->m[k]->col; i++){
        if(colvar->data[i] > colvar->data[j])
          j = i;
        else
          continue;
      }

      if(colvar->data[j] > best_colvar){
        best_colvar = colvar->data[j];
        best_colvar_id = j;
        best_block_id = k;
      }

      DelDVector(&colvar);
    }

    /* copy the vector to the score mx.t_old for computing loadings */
    NewDVector(&t, Eb->m[best_block_id]->row);
    for(i = 0; i < Eb->m[best_block_id]->row; i++){
      t->data[i] = Eb->m[best_block_id]->data[i][best_colvar_id];
    }

    while(1){ /* loop until convergence of t */
      for(k = 0; k < Eb->order; k++){
        NewDVector(&p_b, Eb->m[k]->col);
       /*
        * pb = Xb_T x t / tT x t (step 2)
        */
        CalcBlockLoadings(Eb->m[k], t, p_b);

       /*
        * normalize pb to pb = 1 (step 3)
        */
        DVectNorm(p_b, p_b);

       /*
        * tb = Xb x p_b / cpca_scaling_factor_calculated_in_step_1
        */
        DVectorSet(t_b, 0.f);
        MT_MatrixDVectorDotProduct(Eb->m[k], p_b, t_b);

        for(i = 0; i < t_b->size; i++){
          t_b->data[i] /= model->scaling_factor->data[k];
          T_T->data[k][i] = t_b->data[i]; /*Combine all block scores in T */
        }
        DelDVector(&p_b);
      }

      /* Calculate the Super weights */
      DVectorSet(w_T, 0.f);
      MT_MatrixDVectorDotProduct(T_T, t, w_T);
      mod_t = DVectorDVectorDotProd(t, t);
      for(i = 0; i < w_T->size; i++)
        w_T->data[i] /= mod_t;
      DVectNorm(w_T, w_T);

      /*Calculate the super score*/
      DVectorSet(t_new, 0.f);
      MatrixTranspose(T_T, T);
      MT_MatrixDVectorDotProduct(T, w_T, t_new);
     
      /* check for convergence */
      if(calcConvergence(t_new, t) < CPCACONVERGENCE){
        #ifdef DEBUG
        printf("new score calculated\n");
        printf("pc: %zu\n", pc);

        dvector *tpca = getMatrixColumn(m->scores, pc);
        for(i = 0; i < tpca->size; i++)
          printf("%f %f\n", tpca->data[i], t_new->data[i]);
        DelDVector(&tpca);
        #endif
        /* store the block scores, the super scores and super weights */
        TensorAppendMatrix(model->block_scores, T);

        for(i = 0; i < model->super_scores->row; i++)
          model->super_scores->data[i][pc] = t_new->data[i];

        for(i = 0; i < model->super_weights->row; i++)
          model->super_weights->data[i][pc] = w_T->data[i];

        /* Deflation and calculation block variance explained */
        NewDVector(&local_blockvexp, Eb->order);

        for(k = 0; k < Eb->order; k++){
          /* Calculate the block loadings */
          NewDVector(&p_b, Eb->m[k]->col);
          CalcBlockLoadings(Eb->m[k], t_new, p_b);

          /* store the block of loadings */
          for(j = 0; j < p_b->size; j++)
            model->block_loadings->m[k]->data[j][pc] = p_b->data[j];

          /* Deflation */
          for(i = 0; i < Eb->m[k]->row; i++){
            for(j = 0; j < Eb->m[k]->col; j++){
              //Eb->m[k]->data[i][j] -= t_new->data[i]*p_b->data[j];
              Eb->m[k]->data[i][j] -= t_new->data[i]*model->block_loadings->m[k]->data[j][pc];
            }
          }
          DelDVector(&p_b);
        
       
        
         /*
          * Calculate cumulative percentage explained for each block
          */
          NewMatrix(&Eb_T, Eb->m[k]->col, Eb->m[k]->row);
          MatrixTranspose(Eb->m[k], Eb_T);
          NewMatrix(&Eb_T_E, Eb->m[k]->col, Eb->m[k]->col);
          MatrixDotProduct(Eb_T, Eb->m[k], Eb_T_E); /*SLOW ISNAN TEST +1SEC*/
          local_blockvexp->data[k] = (1.f-(MatrixTrace(Eb_T_E)/tr_orig->data[k]))*100.;
          DelMatrix(&Eb_T_E);
          DelMatrix(&Eb_T);


         /*
          * WARNING: EXPERIMENTAL PART!
          *
          * Calculate the explained variance
          *
          * The sum of the total variance explained
          * will give as results to total variance of a PCA model
          * over the merge of all the blocks.
          *

          bt = getMatrixColumn(T, k);
          double bt_mod = DVectorDVectorDotProd(bt, bt);
          local_blockvexp->data[k] = (bt_mod/block_sum_of_squares->data[k]) * 100.;
          DelDVector(&bt);
          */
        }

        DVectorAppend(model->total_expvar, (mod_t/ss) * 100.);
        DVectorListAppend(model->block_expvar, local_blockvexp);
        DelDVector(&local_blockvexp);
        break;
      }
      else{
        DVectorCopy(t_new, t);
      }
    }
    DelDVector(&t);
  }

  #ifdef DEBUG
  puts("pca EXPLAINED VARIANCE");
  for(pc = 0; pc < npc; pc++){
    printf("VarExp: %f\n", m->varexp->data[pc]);
  }

  puts("Block variance explained at every PC");
  for(k = 0; k < model->expvar->size; k++){
    printf("PC%zu\n", k);
    PrintDVector(model->expvar->d[k]);
  }
  #endif

  #ifdef DEBUG
  DelPCAModel(&m);
  #endif

  DelDVector(&tr_orig);
  DelDVector(&t_new);
  DelMatrix(&T);
  DelMatrix(&T_T);
  DelDVector(&w_T);
  DelMatrix(&X);
  DelTensor(&Eb);
  DelDVector(&t_b);
}

void CPCAScorePredictor(tensor *x,
                        CPCAMODEL *model,
                        size_t npc,
                        matrix *p_super_scores,
                        tensor *p_block_scores)
{
  size_t i ,j, k, pc;
  tensor *Eb;
  matrix *T;
  dvector *p_b;
  dvector *t_b;
  dvector *s_t;
  dvector *w_T;

  if(x->order != model->scaling_factor->size)
    abort();

  if(npc > model->super_scores->col)
    npc = model->super_scores->col;

  NewTensor(&Eb, x->order);
  for(k = 0; k < Eb->order; k++){
    NewTensorMatrix(Eb, k, x->m[k]->row, x->m[k]->col);
    /*Centering Eb*/
    for(j = 0; j < x->m[k]->col; j++){
      for(i = 0; i < x->m[k]->row; i++){
        Eb->m[k]->data[i][j] = x->m[k]->data[i][j] - model->colaverage->d[k]->data[j];
      }
    }

    /*Apply scaling to Eb */
    for(j = 0; j < Eb->m[k]->col; j++){
      if(FLOAT_EQ(model->colscaling->d[k]->data[j], 0.f, EPSILON)){
        for(i = 0; i< Eb->m[k]->row; i++){
          Eb->m[k]->data[i][j] = 0.f;
        }
      }
      else{
        for(i = 0; i < Eb->m[k]->row; i++){
          if(FLOAT_EQ(Eb->m[k]->data[i][j], MISSING, 1e-1)){
            continue;
          }
          else{
            Eb->m[k]->data[i][j] /= model->colscaling->d[k]->data[j];
          }
        }
      }
    }
  }

  NewDVector(&t_b, Eb->m[0]->row);
  NewDVector(&s_t, Eb->m[0]->row);
  NewDVector(&w_T, Eb->order);
  NewMatrix(&T, Eb->m[0]->row, Eb->order);

  /* Prepare the output to store */
  ResizeMatrix(p_super_scores, s_t->size, npc);

  for(pc = 0; pc < npc; pc++){
    for(k = 0; k < Eb->order; k++){
      /* Calculate the block loadings */
      p_b = getMatrixColumn(model->block_loadings->m[k], pc);
     /*
      * tb = Xb x p_b / cpca_scaling_factor_calculated_in_step_1
      */
      DVectNorm(p_b, p_b);

      DVectorSet(t_b, 0.f);
      MT_MatrixDVectorDotProduct(Eb->m[k], p_b, t_b);

      for(i = 0; i < t_b->size; i++){
        t_b->data[i] /= model->scaling_factor->data[k];
        T->data[i][k] = t_b->data[i]; /*Combine all block scores in T */
      }
      DelDVector(&p_b);
    }


    /*Calculate the super score*/
    for(k = 0; k < Eb->order; k++)
      w_T->data[k] = model->super_weights->data[k][pc];
    DVectorSet(s_t, 0.f);
    MT_MatrixDVectorDotProduct(T, w_T, s_t);

    /* Store the block_scores T and the super_score s_t */
    TensorAppendMatrix(p_block_scores, T);

    for(i = 0; i < s_t->size; i++)
      p_super_scores->data[i][pc] = s_t->data[i];

    for(k = 0; k < Eb->order; k++){
      /* Deflation */
      p_b = getMatrixColumn(model->block_loadings->m[k], pc);
      for(i = 0; i < Eb->m[k]->row; i++){
        for(j = 0; j < Eb->m[k]->col; j++){
          Eb->m[k]->data[i][j] -= s_t->data[i]*p_b->data[j];
        }
      }
      DelDVector(&p_b);
    }
  }

  DelTensor(&Eb);
  DelDVector(&w_T);
  DelDVector(&t_b);
  DelDVector(&s_t);
  DelMatrix(&T);
}

void PrintCPCA(CPCAMODEL *m)
{
  size_t k;
  puts("Super scores");
  PrintMatrix(m->super_scores);

  puts("Super weights");
  PrintMatrix(m->super_weights);

  puts("Block scores");
  for(k = 0; k < m->block_scores->order; k++){
    printf("Block %zu\n", k+1);
    PrintMatrix(m->block_scores->m[k]);
  }

  puts("Block loadings");
  for(k = 0; k < m->block_loadings->order; k++){
    printf("Block %zu\n", k+1);
    PrintMatrix(m->block_loadings->m[k]);
  }

  puts("Local Block variance explained at every PC");
  for(k = 0; k < m->block_expvar->size; k++){
    printf("PC%zu\n", k+1);
    PrintDVector(m->block_expvar->d[k]);
  }
  puts("Total variance explained at every PC");
  PrintDVector(m->total_expvar);
}
