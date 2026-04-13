/* Performance benchmarks for cpca algorithms.
 * Copyright (C) 2023-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
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
    tensor *t; // Definition of the input matrix 
    CPCAMODEL *model; // Definition of the PCA model
    int i, j, k;
    int nobj = atoi(argv[1]);
    int nvars = atoi(argv[2]);
    int norder = atoi(argv[3]);
    struct timeval  tv1, tv2;
    NewTensor(&t, norder);
    

    for(k = 0; k < norder; k++){
        NewTensorMatrix(t, k, nobj, nvars);
        srand_(time(0));
        for(i = 0; i < nobj; i++){
            for(j = 0; j < nvars; j++){
                t->m[k]->data[i][j] = randDouble(0,100);
            }
        }
    }
    
    gettimeofday(&tv1, NULL);
    NewCPCAModel(&model);
    CPCA(t, 1, 5, model);
    gettimeofday(&tv2, NULL);
    fprintf( stdout, "Total time (sec): %f\n",
        (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
        (double) (tv2.tv_sec - tv1.tv_sec));
    DelCPCAModel(&model);
    DelTensor(&t);
}


