/* Performance benchmarks for mlr algorithms.
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
    matrix *x;
    matrix *y;
    MLRMODEL *model;
    int i, j;
    int nobj = atoi(argv[1]);
    int nvars = atoi(argv[2]);
    struct timeval  tv1, tv2;
    NewMatrix(&x, nobj, nvars);
    NewMatrix(&y, nobj, 1);
    srand_(time(0));
    
    for(size_t i = 0; i < nobj; i++){
        for(size_t j = 0; j < nvars; j++){
            x->data[i][j] = randDouble(0,20);
        }
        y->data[i][0] = 0.5*x->data[i][0] -3.4*x->data[i][1] +1.2*x->data[i][2]+2*x->data[i][3]+4.2;
     }
    
    gettimeofday(&tv1, NULL);
    NewMLRModel(&model);
    MLR(x, y, model, NULL);
    /*
    PrintPLSModel(model);
    dvector *betas;
    initDVector(&betas);
    PLSBetasCoeff(model, 5, betas);
    PrintDVector(betas);
    DelDVector(&betas);
    */
    gettimeofday(&tv2, NULL);
    fprintf( stdout, "Total time (sec): %f\n",
        (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
        (double) (tv2.tv_sec - tv1.tv_sec));
    DelMLRModel(&model);
    DelMatrix(&x);
    DelMatrix(&y);
}


