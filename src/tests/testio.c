/* testpls.c
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

#include "matrix.h"
#include "modelvalidation.h"
#include "io.h"
#include "pls.h"
#include "numeric.h"

#define ACCEPTABILITY 1e-8

void TestCPCA3(){
    puts("Test 3 - CPCA Model Saving/Reading");
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
    WriteCPCA("cpca.sqlite3", model);

    CPCAMODEL *mread;
    NewCPCAModel(&mread);
    ReadCPCA("cpca.sqlite3", mread);


    for(k = 0; k < model->block_scores->order; k++){
        for(j = 0; j < model->block_scores->m[k]->col; j++){
            for(i = 0; i < model->block_scores->m[k]->row; i++){
                if(FLOAT_EQ(model->block_scores->m[k]->data[i][j], mread->block_scores->m[k]->data[i][j], ACCEPTABILITY)){
                    continue;   
                }
                else{
                    printf("%f %f\n", model->block_scores->m[k]->data[i][j], mread->block_scores->m[k]->data[i][j]);
                    abort();
                }
            }
        }
    }

    for(k = 0; k < model->block_loadings->order; k++){
        for(j = 0; j < model->block_loadings->m[k]->col; j++){
            for(i = 0; i < model->block_loadings->m[k]->row; i++){
                if(FLOAT_EQ(model->block_loadings->m[k]->data[i][j], mread->block_loadings->m[k]->data[i][j], ACCEPTABILITY)){
                    continue;   
                }
                else{
                    printf("%f %f\n", model->block_loadings->m[k]->data[i][j], mread->block_loadings->m[k]->data[i][j]);
                    abort();
                }
            }
        }
    }

    for(j = 0; j < model->super_scores->col; j++){
        for(i = 0; i < model->super_scores->row; i++){
            if(FLOAT_EQ(model->super_scores->data[i][j], mread->super_scores->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", model->super_scores->data[i][j], mread->super_scores->data[i][j]);
                abort();
            }
        }
    }

    for(j = 0; j < model->super_weights->col; j++){
        for(i = 0; i < model->super_weights->row; i++){
            if(FLOAT_EQ(model->super_weights->data[i][j], mread->super_weights->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", model->super_weights->data[i][j], mread->super_weights->data[i][j]);
                abort();
            }
        }
    }

    
    for(j = 0; j < model->colaverage->size; j++){
        for(i = 0; i < model->colaverage->d[j]->size; i++){
            if(FLOAT_EQ(model->colaverage->d[j]->data[i], mread->colaverage->d[j]->data[i], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", model->colaverage->d[j]->data[i], mread->colaverage->d[j]->data[i]);
                abort();
            }
        }
    }

    for(j = 0; j < model->colscaling->size; j++){
        for(i = 0; i < model->colscaling->d[j]->size; i++){
            if(FLOAT_EQ(model->colscaling->d[j]->data[i], mread->colscaling->d[j]->data[i], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", model->colscaling->d[j]->data[i], mread->colscaling->d[j]->data[i]);
                abort();
            }
        }
    }

    for(j = 0; j < model->block_expvar->size; j++){
        for(i = 0; i < model->block_expvar->d[j]->size; i++){
            if(FLOAT_EQ(model->block_expvar->d[j]->data[i], mread->block_expvar->d[j]->data[i], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", model->block_expvar->d[j]->data[i], mread->block_expvar->d[j]->data[i]);
                abort();
            }
        }
    }

    for(i = 0; i < model->scaling_factor->size; i++){
        if(FLOAT_EQ(model->scaling_factor->data[i], mread->scaling_factor->data[i], ACCEPTABILITY)){
            continue;   
        }
        else{
            printf("%f %f\n", model->scaling_factor->data[i], mread->scaling_factor->data[i]);
            abort();
        }
    }

    for(i = 0; i < model->total_expvar->size; i++){
        if(FLOAT_EQ(model->total_expvar->data[i], mread->total_expvar->data[i], ACCEPTABILITY)){
            continue;   
        }
        else{
            printf("%f %f\n", model->total_expvar->data[i], mread->total_expvar->data[i]);
            abort();
        }
    }

    DelCPCAModel(&model);
    DelCPCAModel(&mread);
    DelTensor(&t);
    DelUIVector(&colsizes);
}

void TestPCA2()
{
    printf("Test 2 - PCA Model Saving/Reading\n");
    size_t i, j;
    matrix *m; /* Data matrix */
    PCAMODEL *model, *mread;
    int run = SIGSCIENTIFICRUN;

    NewMatrix(&m, 14, 7);


    m->data[0][0] = 4.0000;  m->data[0][1] = 4.0000;  m->data[0][2] = 1.0000;  m->data[0][3] = 84.1400;  m->data[0][4] = 1.0500;  m->data[0][5] = 235.1500;  m->data[0][6] = 357.1500;
    m->data[1][0] = 5.0000;  m->data[1][1] = 5.0000;  m->data[1][2] = 1.0000;  m->data[1][3] = 79.1000;  m->data[1][4] = 0.9780;  m->data[1][5] = 1.5090;  m->data[1][6] = 231.0000;
    m->data[2][0] = 4.0000;  m->data[2][1] = 5.0000;  m->data[2][2] = 1.0000;  m->data[2][3] = 67.0900;  m->data[2][4] = 0.9700;  m->data[2][5] = 249.0000;  m->data[2][6] = 403.0000;
    m->data[3][0] = 4.0000;  m->data[3][1] = 4.0000;  m->data[3][2] = 1.0000;  m->data[3][3] = 68.0700;  m->data[3][4] = 0.9360;  m->data[3][5] = 187.3500;  m->data[3][6] = 304.5500;
    m->data[4][0] = 3.0000;  m->data[4][1] = 4.0000;  m->data[4][2] = 2.0000;  m->data[4][3] = 68.0800;  m->data[4][4] = 1.0300;  m->data[4][5] = 363.0000;  m->data[4][6] = 529.0000;
    m->data[5][0] = 9.0000;  m->data[5][1] = 7.0000;  m->data[5][2] = 1.0000;  m->data[5][3] = 129.1600;  m->data[5][4] = 1.0900;  m->data[5][5] = 258.0000;  m->data[5][6] = 510.0000;
    m->data[6][0] = 10.0000;  m->data[6][1] = 8.0000;  m->data[6][2] = 0.0000;  m->data[6][3] = 128.1600;  m->data[6][4] = 1.1500;  m->data[6][5] = 352.0000;  m->data[6][6] = 491.0000;
    m->data[7][0] = 6.0000;  m->data[7][1] = 6.0000;  m->data[7][2] = 0.0000;  m->data[7][3] = 78.1118;  m->data[7][4] = 0.8765;  m->data[7][5] = 278.6400;  m->data[7][6] = 353.3000;
    m->data[8][0] = 16.0000;  m->data[8][1] = 10.0000;  m->data[8][2] = 0.0000;  m->data[8][3] = 202.2550;  m->data[8][4] = 1.2710;  m->data[8][5] = 429.1500;  m->data[8][6] = 666.6500;
    m->data[9][0] = 6.0000;  m->data[9][1] = 12.0000;  m->data[9][2] = 0.0000;  m->data[9][3] = 84.1600;  m->data[9][4] = 0.7800;  m->data[9][5] = 279.0000;  m->data[9][6] = 354.0000;
    m->data[10][0] = 4.0000;  m->data[10][1] = 8.0000;  m->data[10][2] = 1.0000;  m->data[10][3] = 72.1100;  m->data[10][4] = 0.8900;  m->data[10][5] = 164.5000;  m->data[10][6] = 339.0000;
    m->data[11][0] = 4.0000;  m->data[11][1] = 9.0000;  m->data[11][2] = 1.0000;  m->data[11][3] = 71.1100;  m->data[11][4] = 0.8660;  m->data[11][5] = 210.0000;  m->data[11][6] = 360.0000;
    m->data[12][0] = 5.0000;  m->data[12][1] = 11.0000;  m->data[12][2] = 1.0000;  m->data[12][3] = 85.1500;  m->data[12][4] = 0.8620;  m->data[12][5] = 266.0000;  m->data[12][6] = 379.0000;
    m->data[13][0] = 5.0000;  m->data[13][1] = 10.0000;  m->data[13][2] = 1.0000;  m->data[13][3] = 86.1300;  m->data[13][4] = 0.8800;  m->data[13][5] = 228.0000;  m->data[13][6] = 361.0000;

    NewPCAModel(&model);

    PCA(m, 1, 5, model, &run);
    WritePCA("pca.sqlite3", model);
    NewPCAModel(&mread);
    ReadPCA("pca.sqlite3", mread);

    /* check original results vs readed results */
    for(j = 0; j < model->scores->col; j++){
        for(i = 0; i < model->scores->row; i++){
            if(FLOAT_EQ(model->scores->data[i][j], mread->scores->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", model->scores->data[i][j], mread->scores->data[i][j]);
                abort();
            }
        }
    }

    for(j = 0; j < model->loadings->col; j++){
        for(i = 0; i < model->loadings->row; i++){
            if(FLOAT_EQ(model->loadings->data[i][j], mread->loadings->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", model->loadings->data[i][j], mread->loadings->data[i][j]);
                abort();
            }
        }
    }


    for(i = 0; i < model->colaverage->size; i++){
        if(FLOAT_EQ(model->colaverage->data[i], mread->colaverage->data[i], ACCEPTABILITY)){
            continue;   
        }
        else{
            printf("%f %f\n", model->colaverage->data[i], mread->colaverage->data[i]);
            abort();
        }
    }

    for(i = 0; i < model->colscaling->size; i++){
        if(FLOAT_EQ(model->colscaling->data[i], mread->colscaling->data[i], ACCEPTABILITY)){
            continue;   
        }
        else{
            printf("%f %f\n", model->colscaling->data[i], mread->colscaling->data[i]);
            abort();
        }
    }

    DelPCAModel(&mread);
    DelPCAModel(&model);
    DelMatrix(&m);
}


void TestPLS1()
{
    printf("Test 1 - PLS Model Saving/Reading\n");
    size_t i, j, k;
    matrix *x, *y; /* Data matrix */
    PLSMODEL *m, *mread;

    NewMatrix(&x, 14, 6);
    NewMatrix(&y, 14, 1);

    x->data[0][0] = 4.0000;  x->data[0][1] = 4.0000;  x->data[0][2] = 1.0000;  x->data[0][3] = 84.1400;  x->data[0][4] = 1.0500;  x->data[0][5] = 235.1500;
    x->data[1][0] = 5.0000;  x->data[1][1] = 5.0000;  x->data[1][2] = 1.0000;  x->data[1][3] = 79.1000;  x->data[1][4] = 0.9780;  x->data[1][5] = 231;
    x->data[2][0] = 4.0000;  x->data[2][1] = 5.0000;  x->data[2][2] = 1.0000;  x->data[2][3] = 67.0900;  x->data[2][4] = 0.9700;  x->data[2][5] = 249.0000;
    x->data[3][0] = 4.0000;  x->data[3][1] = 4.0000;  x->data[3][2] = 1.0000;  x->data[3][3] = 68.0700;  x->data[3][4] = 0.9360;  x->data[3][5] = 187.3500;
    x->data[4][0] = 3.0000;  x->data[4][1] = 4.0000;  x->data[4][2] = 2.0000;  x->data[4][3] = 68.0800;  x->data[4][4] = 1.0300;  x->data[4][5] = 363.0000;
    x->data[5][0] = 9.0000;  x->data[5][1] = 7.0000;  x->data[5][2] = 1.0000;  x->data[5][3] = 129.1600;  x->data[5][4] = 1.0900;  x->data[5][5] = 258.0000;
    x->data[6][0] = 10.0000;  x->data[6][1] = 8.0000;  x->data[6][2] = 0.0000;  x->data[6][3] = 128.1600;  x->data[6][4] = 1.1500;  x->data[6][5] = 352.0000;
    x->data[7][0] = 6.0000;  x->data[7][1] = 6.0000;  x->data[7][2] = 0.0000;  x->data[7][3] = 78.1118;  x->data[7][4] = 0.8765;  x->data[7][5] = 278.6400;
    x->data[8][0] = 16.0000;  x->data[8][1] = 10.0000;  x->data[8][2] = 0.0000;  x->data[8][3] = 202.2550;  x->data[8][4] = 1.2710;  x->data[8][5] = 429.1500;
    x->data[9][0] = 6.0000;  x->data[9][1] = 12.0000;  x->data[9][2] = 0.0000;  x->data[9][3] = 84.1600;  x->data[9][4] = 0.7800;  x->data[9][5] = 279.0000;
    x->data[10][0] = 4.0000;  x->data[10][1] = 8.0000;  x->data[10][2] = 1.0000;  x->data[10][3] = 72.1100;  x->data[10][4] = 0.8900;  x->data[10][5] = 164.5000;
    x->data[11][0] = 4.0000;  x->data[11][1] = 9.0000;  x->data[11][2] = 1.0000;  x->data[11][3] = 71.1100;  x->data[11][4] = 0.8660;  x->data[11][5] = 210.0000;
    x->data[12][0] = 5.0000;  x->data[12][1] = 11.0000;  x->data[12][2] = 1.0000;  x->data[12][3] = 85.1500;  x->data[12][4] = 0.8620;  x->data[12][5] = 266.0000;
    x->data[13][0] = 5.0000;  x->data[13][1] = 10.0000;  x->data[13][2] = 1.0000;  x->data[13][3] = 86.1300;  x->data[13][4] = 0.8800;  x->data[13][5] = 228.0000;


    y->data[0][0] = 357.1500;
    y->data[1][0] = 388.0000;
    y->data[2][0] = 403.0000;
    y->data[3][0] = 304.5500;
    y->data[4][0] = 529.0000;
    y->data[5][0] = 510.0000;
    y->data[6][0] = 491.0000;
    y->data[7][0] = 353.3000;
    y->data[8][0] = 666.6500;
    y->data[9][0] = 354.0000;
    y->data[10][0] = 339.0000;
    y->data[11][0] = 360.0000;
    y->data[12][0] = 379.0000;
    y->data[13][0] = 361.0000;

    /*Allocate the final output*/
    NewPLSModel(&m);

    PLS(x, y, 3, 1, 0, m, NULL);

    /*VALIDATE THE MODEL */
    MODELINPUT minpt = initModelInput();
    minpt.mx = x;
    minpt.my = y;
    minpt.nlv = 3;
    minpt.xautoscaling = 0;
    minpt.yautoscaling = 0;

    BootstrapRandomGroupsCV(&minpt, 3, 100, _PLS_, m->predicted_y, m->pred_residuals, 1, NULL, 0);
    //LeaveOneOut(&minpt, _PLS_, m->predicted_y, m->pred_residuals, 4, NULL, 0);
    PLSRegressionStatistics(y, m->predicted_y, m->q2y, m->sdep, m->bias);

    /*ValidationArg varg = initValidationArg();
    varg.vtype = BootstrapRGCV;
    YScrambling(&minpt, _PLS_, varg, 2, m->yscrambling, 1, NULL);
    */
   
    /* test write pls model*/

    WritePLS("pls.sqlite3", m);
    /* test read pls model*/
    NewPLSModel(&mread);
    ReadPLS("pls.sqlite3", mread);

    /* check original results vs readed results */
    for(j = 0; j < m->xscores->col; j++){
        for(i = 0; i < m->xscores->row; i++){
            if(FLOAT_EQ(m->xscores->data[i][j], mread->xscores->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->xscores->data[i][j], mread->xscores->data[i][j]);
                abort();
            }
        }
    }

    for(j = 0; j < m->xloadings->col; j++){
        for(i = 0; i < m->xloadings->row; i++){
            if(FLOAT_EQ(m->xloadings->data[i][j], mread->xloadings->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->xloadings->data[i][j], mread->xloadings->data[i][j]);
                abort();
            }
        }
    }

    for(j = 0; j < m->xweights->col; j++){
        for(i = 0; i < m->xweights->row; i++){
            if(FLOAT_EQ(m->xweights->data[i][j], mread->xweights->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->xweights->data[i][j], mread->xweights->data[i][j]);
                abort();
            }
        }
    }

    for(j = 0; j < m->yscores->col; j++){
        for(i = 0; i < m->yscores->row; i++){
            if(FLOAT_EQ(m->yscores->data[i][j], mread->yscores->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->yscores->data[i][j], mread->yscores->data[i][j]);
                abort();
            }
        }
    }

    for(j = 0; j < m->yloadings->col; j++){
        for(i = 0; i < m->yloadings->row; i++){
            if(FLOAT_EQ(m->yloadings->data[i][j], mread->yloadings->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->yloadings->data[i][j], mread->yloadings->data[i][j]);
                abort();
            }
        }
    }


    for(j = 0; j < m->recalculated_y->col; j++){
        for(i = 0; i < m->recalculated_y->row; i++){
            if(FLOAT_EQ(m->recalculated_y->data[i][j], mread->recalculated_y->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->recalculated_y->data[i][j], mread->recalculated_y->data[i][j]);
                abort();
            }
        }
    }

    for(j = 0; j < m->recalc_residuals->col; j++){
        for(i = 0; i < m->recalc_residuals->row; i++){
            if(FLOAT_EQ(m->recalc_residuals->data[i][j], mread->recalc_residuals->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->recalc_residuals->data[i][j], mread->recalc_residuals->data[i][j]);
                abort();
            }
        }
    }


    for(j = 0; j < m->predicted_y->col; j++){
        for(i = 0; i < m->predicted_y->row; i++){
            if(FLOAT_EQ(m->predicted_y->data[i][j], mread->predicted_y->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->predicted_y->data[i][j], mread->predicted_y->data[i][j]);
                abort();
            }
        }
    }

    for(j = 0; j < m->pred_residuals->col; j++){
        for(i = 0; i < m->pred_residuals->row; i++){
            if(FLOAT_EQ(m->pred_residuals->data[i][j], mread->pred_residuals->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->pred_residuals->data[i][j], mread->pred_residuals->data[i][j]);
                abort();
            }
        }
    }


    for(j = 0; j < m->r2y_recalculated->col; j++){
        for(i = 0; i < m->r2y_recalculated->row; i++){
            if(FLOAT_EQ(m->r2y_recalculated->data[i][j], mread->r2y_recalculated->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->r2y_recalculated->data[i][j], mread->r2y_recalculated->data[i][j]);
                abort();
            }
        }
    }

    for(j = 0; j < m->r2y_validation->col; j++){
        for(i = 0; i < m->r2y_validation->row; i++){
            if(FLOAT_EQ(m->r2y_validation->data[i][j], mread->r2y_validation->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->r2y_validation->data[i][j], mread->r2y_validation->data[i][j]);
                abort();
            }
        }
    }

    for(j = 0; j < m->q2y->col; j++){
        for(i = 0; i < m->q2y->row; i++){
            if(FLOAT_EQ(m->q2y->data[i][j], mread->q2y->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->q2y->data[i][j], mread->q2y->data[i][j]);
                abort();
            }
        }
    }

    for(j = 0; j < m->sdep->col; j++){
        for(i = 0; i < m->sdep->row; i++){
            if(FLOAT_EQ(m->sdep->data[i][j], mread->sdep->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->sdep->data[i][j], mread->sdep->data[i][j]);
                abort();
            }
        }
    }

    for(j = 0; j < m->sdec->col; j++){
        for(i = 0; i < m->sdec->row; i++){
            if(FLOAT_EQ(m->sdec->data[i][j], mread->sdec->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->sdec->data[i][j], mread->sdec->data[i][j]);
                abort();
            }
        }
    }

    for(j = 0; j < m->bias->col; j++){
        for(i = 0; i < m->bias->row; i++){
            if(FLOAT_EQ(m->bias->data[i][j], mread->bias->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->bias->data[i][j], mread->bias->data[i][j]);
                abort();
            }
        }
    }

    for(j = 0; j < m->roc_auc_recalculated->col; j++){
        for(i = 0; i < m->roc_auc_recalculated->row; i++){
            if(FLOAT_EQ(m->roc_auc_recalculated->data[i][j], mread->roc_auc_recalculated->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->roc_auc_recalculated->data[i][j], mread->roc_auc_recalculated->data[i][j]);
                abort();
            }
        }
    }

    for(j = 0; j < m->roc_auc_validation->col; j++){
        for(i = 0; i < m->roc_auc_validation->row; i++){
            if(FLOAT_EQ(m->roc_auc_validation->data[i][j], mread->roc_auc_validation->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->roc_auc_validation->data[i][j], mread->roc_auc_validation->data[i][j]);
                abort();
            }
        }
    }

    for(j = 0; j < m->precision_recall_ap_recalculated->col; j++){
        for(i = 0; i < m->precision_recall_ap_recalculated->row; i++){
            if(FLOAT_EQ(m->precision_recall_ap_recalculated->data[i][j], mread->precision_recall_ap_recalculated->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->precision_recall_ap_recalculated->data[i][j], mread->precision_recall_ap_recalculated->data[i][j]);
                abort();
            }
        }
    }

    for(j = 0; j < m->precision_recall_ap_validation->col; j++){
        for(i = 0; i < m->precision_recall_ap_validation->row; i++){
            if(FLOAT_EQ(m->precision_recall_ap_validation->data[i][j], mread->precision_recall_ap_validation->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->precision_recall_ap_validation->data[i][j], mread->precision_recall_ap_validation->data[i][j]);
                abort();
            }
        }
    }

    for(j = 0; j < m->yscrambling->col; j++){
        for(i = 0; i < m->yscrambling->row; i++){
            if(FLOAT_EQ(m->yscrambling->data[i][j], mread->yscrambling->data[i][j], ACCEPTABILITY)){
                continue;   
            }
            else{
                printf("%f %f\n", m->yscrambling->data[i][j], mread->yscrambling->data[i][j]);
                abort();
            }
        }
    }

    for(k = 0; k < m->roc_recalculated->order; k++){
        for(j = 0; j < m->roc_recalculated->m[k]->col; j++){
            for(i = 0; i < m->roc_recalculated->m[k]->row; i++){
                if(FLOAT_EQ(m->roc_recalculated->m[k]->data[i][j], mread->roc_recalculated->m[k]->data[i][j], ACCEPTABILITY)){
                    continue;   
                }
                else{
                    printf("%f %f\n", m->roc_recalculated->m[k]->data[i][j], mread->roc_recalculated->m[k]->data[i][j]);
                    abort();
                }
            }
        }
    }

    for(k = 0; k < m->roc_validation->order; k++){
        for(j = 0; j < m->roc_validation->m[k]->col; j++){
            for(i = 0; i < m->roc_validation->m[k]->row; i++){
                if(FLOAT_EQ(m->roc_validation->m[k]->data[i][j], mread->roc_validation->m[k]->data[i][j], ACCEPTABILITY)){
                    continue;   
                }
                else{
                    printf("%f %f\n", m->roc_validation->m[k]->data[i][j], mread->roc_validation->m[k]->data[i][j]);
                    abort();
                }
            }
        }
    }

    for(k = 0; k < m->precision_recall_recalculated->order; k++){
        for(j = 0; j < m->precision_recall_recalculated->m[k]->col; j++){
            for(i = 0; i < m->precision_recall_recalculated->m[k]->row; i++){
                if(FLOAT_EQ(m->precision_recall_recalculated->m[k]->data[i][j], mread->precision_recall_recalculated->m[k]->data[i][j], ACCEPTABILITY)){
                    continue;   
                }
                else{
                    printf("%f %f\n", m->precision_recall_recalculated->m[k]->data[i][j], mread->precision_recall_recalculated->m[k]->data[i][j]);
                    abort();
                }
            }
        }
    }

    for(k = 0; k < m->precision_recall_validation->order; k++){
        for(j = 0; j < m->precision_recall_validation->m[k]->col; j++){
            for(i = 0; i < m->precision_recall_validation->m[k]->row; i++){
                if(FLOAT_EQ(m->precision_recall_validation->m[k]->data[i][j], mread->precision_recall_validation->m[k]->data[i][j], ACCEPTABILITY)){
                    continue;   
                }
                else{
                    printf("%f %f\n", m->precision_recall_validation->m[k]->data[i][j], mread->precision_recall_validation->m[k]->data[i][j]);
                    abort();
                }
            }
        }
    }

    
    for(i = 0; i < m->b->size; i++){
        if(FLOAT_EQ(m->b->data[i], mread->b->data[i], ACCEPTABILITY)){
            continue;   
        }
        else{
            printf("%f %f\n", m->b->data[i], mread->b->data[i]);
            abort();
        }
    }

    for(i = 0; i < m->xvarexp->size; i++){
        if(FLOAT_EQ(m->xvarexp->data[i], mread->xvarexp->data[i], ACCEPTABILITY)){
            continue;   
        }
        else{
            printf("%f %f\n", m->xvarexp->data[i], mread->xvarexp->data[i]);
            abort();
        }
    }

    for(i = 0; i < m->xcolaverage->size; i++){
        if(FLOAT_EQ(m->xcolaverage->data[i], mread->xcolaverage->data[i], ACCEPTABILITY)){
            continue;   
        }
        else{
            printf("%f %f\n", m->xcolaverage->data[i], mread->xcolaverage->data[i]);
            abort();
        }
    }

    for(i = 0; i < m->xcolscaling->size; i++){
        if(FLOAT_EQ(m->xcolscaling->data[i], mread->xcolscaling->data[i], ACCEPTABILITY)){
            continue;   
        }
        else{
            printf("%f %f\n", m->xcolscaling->data[i], mread->xcolscaling->data[i]);
            abort();
        }
    }

    for(i = 0; i < m->ycolaverage->size; i++){
        if(FLOAT_EQ(m->ycolaverage->data[i], mread->ycolaverage->data[i], ACCEPTABILITY)){
            continue;   
        }
        else{
            printf("%f %f\n", m->ycolaverage->data[i], mread->ycolaverage->data[i]);
            abort();
        }
    }

    for(i = 0; i < m->ycolscaling->size; i++){
        if(FLOAT_EQ(m->ycolscaling->data[i], mread->ycolscaling->data[i], ACCEPTABILITY)){
            continue;   
        }
        else{
            printf("%f %f\n", m->ycolscaling->data[i], mread->ycolscaling->data[i]);
            abort();
        }
    }

    DelPLSModel(&mread);
    DelPLSModel(&m);
    DelMatrix(&x);
    DelMatrix(&y);
}

int main(void)
{
  /*test 1- 5*/
  TestPLS1();
  TestPCA2();
  TestCPCA3();
 
  return 0;
}
