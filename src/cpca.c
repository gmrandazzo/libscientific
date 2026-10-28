/* Implements Consensus Principal Component Analysis (CPCA).
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
  (*m) = NULL;
}

/*
* pb = Xb_T x t / tT x t (step 2)
*/
static inline void CalcBlockLoadings(const matrix *Xb, const dvector *t, dvector *p, double mod_t){
  if (mod_t <= EPSILON) {
    for (size_t j = 0; j < p->size; j++) p->data[j] = 0.0;
    return;
  }
  for (size_t j = 0; j < Xb->col; j++) p->data[j] = 0.0;
  for (size_t i = 0; i < Xb->row; i++) {
    double ti = t->data[i];
    const double *row_i = Xb->data[i];
    for (size_t j = 0; j < Xb->col; j++) {
      p->data[j] += row_i[j] * ti;
    }
  }
  for (size_t j = 0; j < Xb->col; j++) {
    p->data[j] /= mod_t;
  }
}

static void cpca_store_and_deflate(tensor *Eb,
                                   size_t pc,
                                   const dvector *t_new,
                                   const dvector *w_T,
                                   matrix *T_T,
                                   matrix *T,
                                   dvector **p_b_list,
                                   const dvector *tr_orig,
                                   double mod_t,
                                   double ss,
                                   CPCAMODEL *model)
{
  MatrixTranspose(T_T, T);
  TensorAppendMatrix(model->block_scores, T);

  for (size_t i = 0; i < model->super_scores->row; i++) {
    model->super_scores->data[i][pc] = t_new->data[i];
  }

  for (size_t i = 0; i < model->super_weights->row; i++) {
    model->super_weights->data[i][pc] = w_T->data[i];
  }

  dvector *local_blockvexp = NULL;
  NewDVector(&local_blockvexp, Eb->order);

  double mod_t_new = DVectorDVectorDotProd((dvector *)t_new, (dvector *)t_new);
  for (size_t k = 0; k < Eb->order; k++) {
    CalcBlockLoadings(Eb->m[k], t_new, p_b_list[k], mod_t_new);

    for (size_t j = 0; j < p_b_list[k]->size; j++) {
      model->block_loadings->m[k]->data[j][pc] = p_b_list[k]->data[j];
    }

    /* Deflation */
    for (size_t i = 0; i < Eb->m[k]->row; i++) {
      double ti = t_new->data[i];
      for (size_t j = 0; j < Eb->m[k]->col; j++) {
        Eb->m[k]->data[i][j] -= ti * p_b_list[k]->data[j];
      }
    }

    /* Calculate cumulative percentage explained for each block */
    double trace_val = 0.0;
    for (size_t i = 0; i < Eb->m[k]->row; i++) {
      for (size_t j = 0; j < Eb->m[k]->col; j++) {
        trace_val += square(Eb->m[k]->data[i][j]);
      }
    }
    local_blockvexp->data[k] = (1.0 - (trace_val / tr_orig->data[k])) * 100.0;
  }

  DVectorAppend(model->total_expvar, (mod_t / ss) * 100.0);
  DVectorListAppend(model->block_expvar, local_blockvexp);
  DelDVector(&local_blockvexp);
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
  dvector *t;
  dvector *t_new;
  dvector *w_T;
  dvector *colvar;
  dvector *tr_orig;
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
    double tr = 0.0;
    for (i = 0; i < Eb->m[k]->row; i++) {
      for (j = 0; j < Eb->m[k]->col; j++) {
        tr += square(Eb->m[k]->data[i][j]);
      }
    }
    tr_orig->data[k] = tr;
  }

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

    /* Pre-allocate vectors for the loop */
    dvector **p_b_list = xmalloc(sizeof(dvector*) * Eb->order);
    for(k=0; k<Eb->order; k++) NewDVector(&p_b_list[k], Eb->m[k]->col);

    while(1){ /* loop until convergence of t */
      mod_t = DVectorDVectorDotProd(t, t);
      for(k = 0; k < Eb->order; k++){
       /*
        * pb = Xb_T x t / tT x t (step 2)
        */
        CalcBlockLoadings(Eb->m[k], t, p_b_list[k], mod_t);

       /*
        * normalize pb to pb = 1 (step 3)
        */
        DVectNorm(p_b_list[k], p_b_list[k]);

       /*
        * tb = Xb x p_b / cpca_scaling_factor_calculated_in_step_1
        */
        double inv_scale = 1.0 / model->scaling_factor->data[k];
        for(i = 0; i < Eb->m[k]->row; i++){
          double sum = 0.0;
          const double *row_i = Eb->m[k]->data[i];
          for(j = 0; j < Eb->m[k]->col; j++){
            sum += row_i[j] * p_b_list[k]->data[j];
          }
          T_T->data[k][i] = sum * inv_scale; /* Combine all block scores in T */
        }
      }

      /* Calculate the Super weights */
      double inv_mod = (mod_t > EPSILON) ? (1.0 / mod_t) : 0.0;
      for(k = 0; k < Eb->order; k++){
        double sum = 0.0;
        const double *row_k = T_T->data[k];
        for(i = 0; i < T_T->col; i++){
          sum += row_k[i] * t->data[i];
        }
        w_T->data[k] = sum * inv_mod;
      }
      DVectNorm(w_T, w_T);

      /*Calculate the super score*/
      DVectorSet(t_new, 0.f);
      for(k = 0; k < Eb->order; k++) {
        double w = w_T->data[k];
        const double *row_k = T_T->data[k];
        for(i = 0; i < Eb->m[0]->row; i++) {
          t_new->data[i] += row_k[i] * w;
        }
      }

      /* check for convergence */
      if(calcConvergence(t_new, t) < CPCACONVERGENCE){
        #ifdef DEBUG
        printf("new score calculated\n");
        printf("pc: %zu\n", pc);
        #endif
        cpca_store_and_deflate(Eb, pc, t_new, w_T, T_T, T, p_b_list, tr_orig, mod_t, ss, model);
        break;
      }
      else{
        DVectorCopy(t_new, t);
      }
    }

    for(k=0; k<Eb->order; k++) DelDVector(&p_b_list[k]);
    xfree(p_b_list);

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

  NewDVector(&s_t, Eb->m[0]->row);
  NewDVector(&w_T, Eb->order);
  NewMatrix(&T, Eb->m[0]->row, Eb->order);

  /* Pre-allocate normalized loading vectors for each block */
  dvector **p_b_list = xmalloc(sizeof(dvector*) * Eb->order);
  for(k = 0; k < Eb->order; k++) {
    NewDVector(&p_b_list[k], model->block_loadings->m[k]->row);
  }

  /* Prepare the output to store */
  ResizeMatrix(p_super_scores, s_t->size, npc);

  for(pc = 0; pc < npc; pc++){
    for(k = 0; k < Eb->order; k++){
      /* Calculate the normalized block loadings */
      const matrix *loadings_k = model->block_loadings->m[k];
      for(j = 0; j < loadings_k->row; j++){
        p_b_list[k]->data[j] = loadings_k->data[j][pc];
      }
      DVectNorm(p_b_list[k], p_b_list[k]);

      /* tb = Xb x p_b / scaling_factor */
      double inv_scale = 1.0 / model->scaling_factor->data[k];
      for(i = 0; i < Eb->m[k]->row; i++){
        double sum = 0.0;
        const double *row_i = Eb->m[k]->data[i];
        for(j = 0; j < Eb->m[k]->col; j++){
          sum += row_i[j] * p_b_list[k]->data[j];
        }
        T->data[i][k] = sum * inv_scale;
      }
    }

    /* Calculate the super score: s_t = T * w_T */
    for(k = 0; k < Eb->order; k++) {
      w_T->data[k] = model->super_weights->data[k][pc];
    }
    for(i = 0; i < T->row; i++){
      double sum = 0.0;
      const double *row_i = T->data[i];
      for(k = 0; k < T->col; k++){
        sum += row_i[k] * w_T->data[k];
      }
      s_t->data[i] = sum;
      p_super_scores->data[i][pc] = sum;
    }

    /* Store the block_scores T */
    TensorAppendMatrix(p_block_scores, T);

    for(k = 0; k < Eb->order; k++){
      /* Deflation */
      const matrix *loadings_k = model->block_loadings->m[k];
      for(i = 0; i < Eb->m[k]->row; i++){
        double s = s_t->data[i];
        double *row_i = Eb->m[k]->data[i];
        for(j = 0; j < Eb->m[k]->col; j++){
          row_i[j] -= s * loadings_k->data[j][pc];
        }
      }
    }
  }

  for(k = 0; k < Eb->order; k++) {
    DelDVector(&p_b_list[k]);
  }
  xfree(p_b_list);

  DelTensor(&Eb);
  DelDVector(&w_T);
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
