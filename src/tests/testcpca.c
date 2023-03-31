#include <stdio.h>
#include <math.h>

#include "cpca.h"
#include "pca.h"
#include "tensor.h"
#include "numeric.h"

void PrepareMatrix4PCA(tensor *t, CPCAMODEL *model, matrix *X)
{
  size_t i, j, k, c, row, col;
  col = 0;
  row = 0;
  for(k = 0; k < t->order; k++){
    row = t->m[k]->row;
    col += t->m[k]->col;
  }

  ResizeMatrix(X, row, col);

  c = 0;
  for(k = 0; k < t->order; k++){
    double m = model->scaling_factor->data[k];
    for(j = 0; j < t->m[k]->col; j++){
      for(i = 0; i < t->m[k]->row; i++){
        X->data[i][c] = (t->m[k]->data[i][j] - model->colaverage->d[k]->data[j]);
        X->data[i][c] /= model->colscaling->d[k]->data[j];
        X->data[i][c] /= m; /* block scaling N.B.: this is applied to reproduce the pca score */
      }
      c += 1;
    }
  }
}

void test3()
{
  puts("Test 3: Projection test");
  size_t i, j, k;
  size_t rowsize=20;
  uivector *colsizes;
  tensor *t;
  matrix *p_super_scores;
  tensor *p_block_scores;

  initUIVector(&colsizes);
  UIVectorAppend(colsizes, 10);
  UIVectorAppend(colsizes, 7);
  UIVectorAppend(colsizes, 13);

  NewTensor(&t, colsizes->size);
  for(k = 0; k < colsizes->size; k++){
    NewTensorMatrix(t, k, rowsize, colsizes->data[k]);
    srand(rowsize+colsizes->data[k]);
    for(i = 0; i < rowsize; i++){
      for(j = 0; j < colsizes->data[k]; j++){
        t->m[k]->data[i][j] = rand() % 100;
      }
    }
  }

  CPCAMODEL *model;
  NewCPCAModel(&model);

  CPCA(t, 1, 5, model);
  //PrintCPCA(model);

  initMatrix(&p_super_scores);
  initTensor(&p_block_scores);
  CPCAScorePredictor(t,
                     model,
                     5,
                     p_super_scores,
                     p_block_scores);

  for(i = 0; i < model->super_scores->row; i++){
    for(j = 0; j < model->super_scores->col; j++){
      if(FLOAT_EQ(model->super_scores->data[i][j], p_super_scores->data[i][j], 1E-6)){
        continue;
      }
      else{
        abort();
      }
    }
  }

  for(k = 0; k < model->block_scores->order; k++){
    for(i = 0; i < model->block_scores->m[k]->row; i++){
      for(j = 0; j < model->block_scores->m[k]->col; j++){
        if(FLOAT_EQ(model->block_scores->m[k]->data[i][j], p_block_scores->m[k]->data[i][j], 1E-6)){
          continue;
        }
        else{
          abort();
        }
      }
    }
  }


  DelMatrix(&p_super_scores);
  DelTensor(&p_block_scores);
  DelCPCAModel(&model);
  DelTensor(&t);
  DelUIVector(&colsizes);
}

void test2()
{
  puts("Test 2: CPCA Variance test!");
 /*
  * we consider four blocks of five variables where one block
  * contains a strong direction that is not available in the other blocks.
  * X1 = [d1 d1 d1 d1 d1 ]
  * X2 = [d2 randn(4)]
  * X3 = [d2 randn(4)],
  * X4 = [d2 randn(4)]
  * All blocks have 50 observations
  * randn = random data
  *
  */
  size_t i, j, k;
  size_t rowsize=10;
  size_t colsize=5;
  size_t nblocks=4;
  tensor *t;

  /* d1 and d2 contains the cosine direction of the principal component */
  double d1 = -0.44;
  double d2 = +0.44;

  srand_(rowsize+colsize+nblocks);
  NewTensor(&t, nblocks);
  for(k = 0; k < nblocks; k++){
    NewTensorMatrix(t, k, rowsize, colsize);
    for(i = 0; i < rowsize; i++){
      for(j = 0; j < colsize; j++){
        t->m[k]->data[i][j] = randDouble(-1, 1);
      }
    }
  }

  for(i = 0; i < rowsize; i++){
    for(j = 0; j < colsize; j++){
      t->m[0]->data[i][j] = (t->m[0]->data[i][j]*d1*d1);
    }

    t->m[1]->data[i][0] = (t->m[1]->data[i][0]*d2*d2);
    t->m[2]->data[i][0] = (t->m[2]->data[i][0]*d2*d2);
    t->m[3]->data[i][0] = (t->m[3]->data[i][0]*d2*d2);
  }


  CPCAMODEL *model;
  NewCPCAModel(&model);

  CPCA(t, 1, 5, model);
  //PrintCPCA(model);

  DelCPCAModel(&model);
  DelTensor(&t);
}

void test1(){
  puts("Test 1: Algorithm Check. CPCA super scores = PCA scores");
  size_t i, j, k;
  size_t rowsize=20;
  uivector *colsizes;
  tensor *t;
  size_t npc = 5;
  initUIVector(&colsizes);
  UIVectorAppend(colsizes, 10);
  UIVectorAppend(colsizes, 7);
  UIVectorAppend(colsizes, 13);

  NewTensor(&t, colsizes->size);
  for(k = 0; k < colsizes->size; k++){
    NewTensorMatrix(t, k, rowsize, colsizes->data[k]);
    srand(rowsize+colsizes->data[k]);
    for(i = 0; i < rowsize; i++){
      for(j = 0; j < colsizes->data[k]; j++){
        t->m[k]->data[i][j] = rand() % 100;
      }
    }
  }

  CPCAMODEL *model;
  NewCPCAModel(&model);

  CPCA(t, 1, npc, model);
  //PrintCPCA(model);

  PCAMODEL *pca;
  matrix *X;
  initMatrix(&X);
  PrepareMatrix4PCA(t, model, X);

  NewPCAModel(&pca);
  PCA(X, -1, npc, pca, NULL);

  /* COMPARE THE SUPERSCORE WITH THE PCA RESULTS  */
  for(j = 0; j < model->super_scores->col; j++){
    for(i = 0; i < model->super_scores->row; i++){
      if(FLOAT_EQ(fabs(model->super_scores->data[i][j]), fabs(pca->scores->data[i][j]), 1E-1)){
        continue;
      }
      else{
        printf("%f =! %f\n", model->super_scores->data[i][j], pca->scores->data[i][j]);
        abort();
      }
    }
  }

  DelMatrix(&X);
  DelPCAModel(&pca);
  DelCPCAModel(&model);
  DelTensor(&t);
  DelUIVector(&colsizes);
}

int main(void){
  test1();
  test2();
  test3();
  return 0;
}
