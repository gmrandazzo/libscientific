/* testpca.c
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

#include <stdlib.h>
#include <stdlib.h>
#include "numeric.h"
#include "preprocessing.h"
#include "pca.h"
#include "pca.h"
#include "datasets.h"
#include "scientificinfo.h"

void test3()
{
  puts("Test3: Matrix Whitening");
  matrix *m;
  matrix *w;
  matrix *m_w;
  matrix *cov;
  size_t i;
  size_t j;

  NewMatrix(&m, 4, 3);
  for(i = 0 ; i < 4; i++){
    for(j = 0; j < 3; j++){
      m->data[i][j] = randDouble(-2, 2);
    }
  }
  initMatrix(&w);
  initMatrix(&m_w);
  MatrixWhitening(m, w, m_w);
  puts("whitening");
  PrintMatrix(w);

  /* test
   * the feature covariance matrix of the whitened signal
   * should equal the identity matrix
   */
  initMatrix(&cov);
  MatrixCovariance(m_w, cov);
  puts("covariance");
  PrintMatrix(cov);
  for(i = 0; i < cov->row; i++){
    for(j = 0; j < cov->col; j++){
      if(i == j){
        if(FLOAT_EQ(cov->data[i][j], 1.f, 1e-2))
          continue;
        else
          abort();
      }
      else{
        if(FLOAT_EQ(cov->data[i][j], 0.f, 1e-2))
          continue;
        else
          abort();
      }
    }
  }
  DelMatrix(&cov);
  DelMatrix(&m);
  DelMatrix(&m_w);
  DelMatrix(&w);
}
void test2()
{
    puts("Test1: MatrixPreprocess");
    size_t i;
    size_t j;
    size_t k;
    double v;
    matrix *m;
    matrix *p;
    dvector *cavg;
    dvector *cstdev;
    NewMatrix(&m, 4, 3);
    k = 1;
    for(i = 0 ; i < 4; i++){
      for(j = 0; j < 3; j++){
        m->data[i][j] = k;
        k++;
      }
    }

    NewMatrix(&p, 4, 3);
    initDVector(&cavg);
    initDVector(&cstdev);
    MatrixPreprocess(m, 3, cavg, cstdev, p);
    for(i = 0; i < cavg->size; i++){
      if(FLOAT_EQ(cavg->data[i], 5.5+i, 1e-2)){
        continue;
      }
      else{
        abort();
      }
    }

    for(i = 0; i < cstdev->size; i++){
      if(FLOAT_EQ(cstdev->data[i], 1.9680, 1e-4)){
        continue;
      }
      else{
        abort();
      }   
    }
    
    for(i = 0; i < 4; i++){
      if(i == 0){
        v = -2.287;
      }
      else if(i == 1){
        v = -0.762;
      }
      else if(i == 2){
        v = 0.762;
      }
      else{
        v = 2.287;
      }
      for(j = 0; j < 3; j++){
        if(FLOAT_EQ(p->data[i][j], v, 1e-3)){
          continue;
        }
        else{
          abort();
        }
      }
    }
    puts("MatrixPreprocess Pareto Scaling: OK");

    DelDVector(&cstdev);
    DelDVector(&cavg);
    DelMatrix(&p);
    DelMatrix(&m);
}

void test1()
{
    puts("Test1: MatrixPreprocess");
    size_t i;
    size_t j;
    size_t k;
    double v;
    matrix *m;
    matrix *p;
    dvector *cavg;
    dvector *cstdev;
    NewMatrix(&m, 4, 3);
    k = 1;
    for(i = 0 ; i < 4; i++){
      for(j = 0; j < 3; j++){
        m->data[i][j] = k;
        k++;
      }
    }

    NewMatrix(&p, 4, 3);
    initDVector(&cavg);
    initDVector(&cstdev);
    MatrixPreprocess(m, 1, cavg, cstdev, p);

    for(i = 0; i < cavg->size; i++){
      if(FLOAT_EQ(cavg->data[i], 5.5+i, 1e-2)){
        continue;
      }
      else{
        abort();
      }
    }

    for(i = 0; i < cstdev->size; i++){
      if(FLOAT_EQ(cstdev->data[i], 3.8730, 1e-4)){
        continue;
      }
      else{
        abort();
      }   
    }
    
    for(i = 0; i < 4; i++){
      if(i == 0){
        v = -1.162;
      }
      else if(i == 1){
        v = -0.387;
      }
      else if(i == 2){
        v = 0.387;
      }
      else{
        v = 1.162;
      }
      for(j = 0; j < 3; j++){
        if(FLOAT_EQ(p->data[i][j], v, 1e-3)){
          continue;
        }
        else{
          abort();
        }
      }
    }
    puts("MatrixPreprocess STEDV Scaling: OK");

    DelDVector(&cstdev);
    DelDVector(&cavg);
    DelMatrix(&p);
    DelMatrix(&m);
}

int main(void)
{
  test1();
  test2();
  /*test3();*/
  return 0;
}
