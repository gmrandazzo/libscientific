/* Performance benchmarks for pca algorithms.
 * Copyright (C) 2020-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
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

#include <stdio.h>
#include <scientific.h>
#include <sys/time.h>

int main(int argc, char **argv)
{
    matrix *m; // Definition of the input matrix 
    PCAMODEL *model; // Definition of the PCA model
    int i, j;
    int nobj = atoi(argv[1]);
    int nvars = atoi(argv[2]);
    struct timeval  tv1, tv2;
    NewMatrix(&m, nobj, nvars);

    srand_(time(0));
    for(size_t i = 0; i < nobj; i++){
        for(size_t j = 0; j < nvars; j++){
            m->data[i][j] = randDouble(0,20);
        }
    }
    
    gettimeofday(&tv1, NULL);
    NewPCAModel(&model); // Allocation of the PCA model
    PCA(m, 1, 5, model, NULL); // Calculation of the PCA on matrix m using unit variance scaling (1) and the extraction of 5 principal components 
    gettimeofday(&tv2, NULL);
    fprintf( stdout, "Total time (sec): %f\n",
        (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
        (double) (tv2.tv_sec - tv1.tv_sec));
    DelPCAModel(&model);
    DelMatrix(&m);
}


