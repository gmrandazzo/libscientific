/* metricspace.c
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

#include "metricspace.h"
#include "numeric.h"
#include <math.h>


void EuclideanDistance(matrix* m1, matrix* m2, matrix **distances)
{
  if(m1->col == m2->col){
    size_t i, j, k;
    double dist;
    
    ResizeMatrix(distances, m2->row, m1->row); /* each column is a distance that correspond to m1->row */
    
    for(i = 0; i < m1->row; i++){
      for(k = 0; k < m2->row; k++){
        dist = 0.f;
        for(j = 0; j < m1->col; j++){ /* is the same of for(j = 0; j < m2->col; j++){ */
          dist += square(m1->data[i][j] - m2->data[k][j]);
        }
        (*distances)->data[k][i] = sqrt(dist);
      }
    }
  }
  else{
    fprintf(stderr, "Unable to compute Euclidean Distance. The number of variables differ\n");
    fflush(stderr);
    abort();
  }
}

void SquaredEuclideanDistance(matrix *m1, matrix *m2, matrix **distances)
{
  if(m1->col == m2->col){
    size_t i, j, k;
    double dist;
    
    ResizeMatrix(distances, m2->row, m1->row); /* each column is a distance that correspond to m1->row */
    
    for(i = 0; i < m1->row; i++){
      for(k = 0; k < m2->row; k++){
        dist = 0.f;
        for(j = 0; j < m1->col; j++){ /* is the same of for(j = 0; j < m2->col; j++){ */
          dist += square(m1->data[i][j] - m2->data[k][j]);
        }
        (*distances)->data[k][i] = dist;
      }
    }
  }
  else{
    fprintf(stderr, "Unable to compute Euclidean Distance. The number of variables differ\n");
    fflush(stderr);
    abort();
  }
}

void ManhattanDistance(matrix* m1, matrix* m2, matrix** distances)
{
  if(m1->col == m2->col){
    size_t i, j, k;
    double dist;
    
    ResizeMatrix(distances, m2->row, m1->row); /* each column is a distance that correspond to m1->row */
    
    for(i = 0; i < m1->row; i++){
      for(k = 0; k < m2->row; k++){
        dist = 0.f;
        for(j = 0; j < m1->col; j++){ /* is the same of for(j = 0; j < m2->col; j++){ */
          dist += fabs(m1->data[i][j] - m2->data[k][j]);
        }
        (*distances)->data[k][i] = dist;
      }
    }
  }
  else{
    fprintf(stderr, "Unable to compute Manhattan Distance. The number of variables differ\n");
    fflush(stderr);
    abort();
  } 
}

void CosineDistance(matrix* m1, matrix* m2, matrix** distances)
{
  if(m1->col == m2->col){
    size_t i, j, k;
    double n, d_a, d_b;
    
    ResizeMatrix(distances, m2->row, m1->row); /* each column is a distance that correspond to m1->row */
    
    for(i = 0; i < m1->row; i++){
      for(k = 0; k < m2->row; k++){
        n = 0.f; d_a = 0.f; d_b = 0.f;
        for(j = 0; j < m1->col; j++){ /* is the same of for(j = 0; j < m2->col; j++){ */
          n += m1->data[i][j] * m2->data[k][j];
          d_a += square(m1->data[i][j]);
          d_b += square(m2->data[k][j]);
        }
        (*distances)->data[k][i] = n/(sqrt(d_a)*sqrt(d_b));
      }
    }
  }
  else{
    fprintf(stderr, "Unable to compute Cosine Distance. The number of variables differ\n");
    fflush(stderr);
    abort();
  } 
}

double MahalanobisDistance(matrix* g1, matrix* g2)
{
  if(g1->col == g2->col){
    size_t i, j;
    double dist, a, b;
    dvector *mean1, *mean2, *xdiff, *x;
    matrix *cov1, *cov2, *invcov;
    
    initDVector(&mean1);
    MatrixColAverage(g1, &mean1);
    
    initDVector(&mean2);
    MatrixColAverage(g2, &mean2);
    
    initMatrix(&cov1);
    MatrixCovariance(g1, &cov1);
    
    initMatrix(&cov2);
    MatrixCovariance(g2, &cov2);
    
    /* pooled covariance matrix in cov1*/
    for(i = 0; i < cov1->row; i++){
      for(j = 0; j < cov1->col; j++){
        a = (double)g1->row/(double)(g1->row+g2->row);
        b = (double)g2->row/(double)(g1->row+g2->row);
        cov1->data[i][j] = (cov1->data[i][j] * a) + (cov2->data[i][j] * b);
      }
    }
    
    initMatrix(&invcov);
    
    MatrixInversion(cov1, &invcov);
    
    NewDVector(&xdiff, mean1->size);
    
    for(i = 0; i < mean1->size; i++){
      xdiff->data[i] = mean1->data[i] - mean2->data[i];
    }
    
    NewDVector(&x, mean1->size);
    DVectorMatrixDotProduct(invcov, xdiff, x);
    
    dist = sqrt(DVectorDVectorDotProd(xdiff, x));
    
    DelMatrix(&invcov);
    DelDVector(&mean2);
    DelDVector(&mean1);
    DelDVector(&xdiff);
    DelMatrix(&cov2);
    DelMatrix(&cov1);
    DelDVector(&x);
    
    return dist;
  }
  else{
    fprintf(stderr, "Unable to compute Euclidean Distance. The number of variables differ\n");
    fflush(stderr);
    abort();
    return 0.f;
  }
}
