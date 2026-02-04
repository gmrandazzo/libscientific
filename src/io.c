/* Handles input/output operations and file formats.
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
#include <sqlite3.h>

#include "memwrapper.h"
#include "vector.h"
#include "matrix.h"
#include "pca.h"
#include "pls.h"
#include "cpca.h"
#include "mlr.h"
#include "lda.h"
#include "ica.h"
#include "upca.h"
#include "upls.h"

static void serialize_uivector(uivector *v, dvector *serialized)
{
    size_t i;
    DVectorResize(serialized, v->size + 1);
    serialized->data[0] = (double)v->size;
    for(i = 0; i < v->size; i++) {
        serialized->data[i+1] = (double)v->data[i];
    }
}

static void deserialize_uivector(dvector *serialized, uivector *v)
{
    size_t i;
    size_t size = (size_t)serialized->data[0];
    UIVectorResize(v, size);
    for(i = 0; i < size; i++) {
        v->data[i] = (size_t)serialized->data[i+1];
    }
}

static void serialize_tensor(tensor *tx, dvector *tensor_serialized)
{
    size_t i, j, k, c;
    size_t tot_size = 1;
    for(k = 0; k < tx->order; k++){
        tot_size += 2+(tx->m[k]->row*tx->m[k]->col);        
    }
    DVectorResize(tensor_serialized, tot_size);
    c = 0;
    tensor_serialized->data[c++] = tx->order;
    for(k = 0; k < tx->order; k++){
        tensor_serialized->data[c++] = tx->m[k]->row;
        tensor_serialized->data[c++] = tx->m[k]->col;
        for(i = 0; i < tx->m[k]->row; i++){
            for(j = 0; j < tx->m[k]->col; j++){
                tensor_serialized->data[c++] = tx->m[k]->data[i][j];
            }
        }
    }
}

static void deserialize_tensor(dvector *tensor_serialized, tensor *tx)
{
    size_t i, j, k, c = 0;
    size_t order = tensor_serialized->data[c++];
    size_t nrow, ncol;
    for(k = 0; k < order; k++){
        nrow = tensor_serialized->data[c++];
        ncol = tensor_serialized->data[c++];
        AddTensorMatrix(tx, nrow, ncol);
        for(i = 0; i < nrow; i++){
            for(j = 0; j < ncol; j++){
                tx->m[k]->data[i][j] = tensor_serialized->data[c++];
            }
        }
    }
}

static void serialize_matrix(matrix *mx, dvector *matrix_serialized)
{
    size_t i, j, c = 0;
    DVectorResize(matrix_serialized, 2+mx->row*mx->col);
    matrix_serialized->data[c++]= mx->row;
    matrix_serialized->data[c++]= mx->col;
    for(i = 0; i < mx->row; i++){
        for(j = 0; j < mx->col; j++){
            matrix_serialized->data[c++] = mx->data[i][j];
        }
    }
}

static void deserialize_matrix(dvector *matrix_serialized, matrix *mx)
{
    size_t i, j, c = 0;
    size_t nrow = (size_t)matrix_serialized->data[c++];
    size_t ncol = (size_t)matrix_serialized->data[c++]; 
    ResizeMatrix(mx, nrow, ncol);
    for(i = 0; i < mx->row; i++){
        for(j = 0; j < mx->col; j++){
            mx->data[i][j] = matrix_serialized->data[c++];
        }
    }
}

static void serialize_dvectorlist(dvectorlist *dvlst, dvector *dlst_serialized)
{
    size_t i, j, c = 0;
    size_t tot_size = 0;
    for(i = 0; i < dvlst->size; i++){
        tot_size += dvlst->d[i]->size+1;
    }
    DVectorResize(dlst_serialized, tot_size);
    for(i = 0; i < dvlst->size; i++){
        dlst_serialized->data[c++] = dvlst->d[i]->size;
        for(j = 0; j < dvlst->d[i]->size; j++){
            dlst_serialized->data[c++] = dvlst->d[i]->data[j];
        }
    }
}


static void deserialize_dvectorlist(dvector *dlst_serialized, dvectorlist *dlst)
{
    size_t j, c = 0;
    dvector *v;
    while(c < dlst_serialized->size){
        NewDVector(&v, (size_t)dlst_serialized->data[c++]);
        for(j = 0; j < v->size; j++){
            v->data[j] = dlst_serialized->data[c++];
        }
        DVectorListAppend(dlst, v);
        DelDVector(&v);
    }
}

static int callback(void *data, int argc, char **argv, char **azColName){
    #ifdef DEBUG
    int i;
    fprintf(stderr, "%s: ", (const char*)data);

    for(i = 0; i < argc; i++){
        printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
    }

    printf("\n");
    #endif
    return 0;

}

static void write_vector_into_sqltable(sqlite3 *db, char *tabname, dvector *vect)
{
    size_t i;
    int rc;
    int length;
    char *err_msg = 0;
    char *sql;
    sqlite3_stmt* stmt;

    /* Create SQL statement */
    length = snprintf(NULL, 0, "CREATE TABLE IF NOT EXISTS %s (id INTEGER PRIMARY KEY AUTOINCREMENT, value REAL);", tabname);
    sql = xmalloc(length+1);
    snprintf(sql, length+1, "CREATE TABLE IF NOT EXISTS %s (id INTEGER PRIMARY KEY AUTOINCREMENT, value REAL);", tabname);
    #ifdef DEBUG
    printf("%s\n", sql);
    #endif

    /* Execute SQL statement */
    rc = sqlite3_exec(db, sql, callback, 0, &err_msg);
    if(rc != SQLITE_OK){
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
        xfree(sql);
        return;
    }
    #ifdef DEBUG
    else{
        fprintf(stdout, "Table created successfully\n");
    }
    #endif
    xfree(sql);

    /* Prepare INSERT statement once */
    length = snprintf(NULL, 0, "INSERT INTO %s (value) VALUES (?);", tabname);
    sql = xmalloc(length+1);
    snprintf(sql, length+1, "INSERT INTO %s (value) VALUES (?);", tabname);
    
    rc = sqlite3_prepare_v2(db, sql, -1, &stmt, 0);
    xfree(sql);

    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQL error (Prepare): %s\n", sqlite3_errmsg(db));
        return;
    }

    for (i = 0; i < vect->size; i++){
        sqlite3_bind_double(stmt, 1, vect->data[i]);
        rc = sqlite3_step(stmt);
        if (rc != SQLITE_DONE) {
            fprintf(stderr, "SQL error (Step): %s\n", sqlite3_errmsg(db));
        }
        sqlite3_reset(stmt);
    }
    sqlite3_finalize(stmt);
}

static void read_vector(sqlite3 *db, char *tabname, dvector *vect)
{
 /* SQL query to select the column data */
    int rc;
    char *sql;
    int lenght;
    sqlite3_stmt* stmt;

    lenght = snprintf(NULL, 0, "SELECT value FROM %s;", tabname);
    sql = xmalloc(lenght+1);
    snprintf(sql,lenght+1, "SELECT value FROM %s;", tabname);
    #ifdef DEBUG
    printf("%s\n", sql);
    #endif
    rc = sqlite3_prepare_v2(db, sql, -1, &stmt, 0);
    if (rc != SQLITE_OK){
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        xfree(sql);
        abort();
    }

    /* Fetch and store the data in the dynamic array */
    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
        double value = sqlite3_column_double(stmt, 0);
        DVectorAppend(vect, value);
    }

    if (rc != SQLITE_DONE) {
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
        sqlite3_finalize(stmt);
        sqlite3_close(db);
    }
    xfree(sql);
}

static void OpenDB(char *dbpath, sqlite3 **db)
{
    int rc;
    /* Open the SQLite database */
    rc = sqlite3_open(dbpath, db);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg((*db)));
        sqlite3_close((*db));
        return;
    }
    #ifdef DEBUG
    else{
        fprintf(stdout, "Opened database successfully\n");
    }
    #endif
}

static void DropAllTables(sqlite3 *db)
{
    int rc;
    char *err_msg = 0;
    const char *dropAllObjectsSQL = "SELECT 'DROP TABLE IF EXISTS ' || name || ';' FROM sqlite_master WHERE type = 'table';";
    /* Execute SQL statement */
    rc = sqlite3_exec(db, dropAllObjectsSQL, 0, 0, &err_msg);
    if(rc != SQLITE_OK){
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
    }
    #ifdef DEBUG
    else{
        fprintf(stdout, "Table created successfully\n");
    }
    #endif
}

static void CloseDB(sqlite3 *db)
{
    sqlite3_close(db);
}

void WritePCA(char *dbpath, PCAMODEL *pca)
{
    dvector *serial_vect;
    sqlite3 *db = NULL;
    char *err_msg = 0;

    OpenDB(dbpath, &db);
    sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, &err_msg);
    DropAllTables(db);
    /**
     * PCA model data structure.
     * 
     * - **scores** matrix of scores
     * - **loadings** matrix of loadings
     * - **varexp** vector of explained variance by every component 
     * - **colaverage** input matrix column average
     * - **colscaling** input matrix column scaling
     */

    write_vector_into_sqltable(db, "colaverage", pca->colaverage);
    write_vector_into_sqltable(db, "colscaling", pca->colscaling);
    write_vector_into_sqltable(db, "varexp", pca->varexp);
    
    initDVector(&serial_vect);
    serialize_matrix(pca->scores, serial_vect);
    write_vector_into_sqltable(db, "scores", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(pca->loadings, serial_vect);
    write_vector_into_sqltable(db, "loadings", serial_vect);
    DelDVector(&serial_vect);
    
    sqlite3_exec(db, "COMMIT;", NULL, NULL, &err_msg);
    CloseDB(db);
}

void ReadPCA(char *dbpath, PCAMODEL *pca)
{
    dvector *serial_vect;
    sqlite3* db = NULL;

    OpenDB(dbpath, &db);

    read_vector(db, "colaverage", pca->colaverage);
    read_vector(db, "colscaling", pca->colscaling);
    read_vector(db, "varexp", pca->varexp);

    initDVector(&serial_vect);
    read_vector(db, "scores", serial_vect);
    deserialize_matrix(serial_vect, pca->scores);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "loadings", serial_vect);
    deserialize_matrix(serial_vect, pca->loadings);
    DelDVector(&serial_vect);

    CloseDB(db);
}

void WriteCPCA(char *dbpath, CPCAMODEL *cpca)
{
    dvector *serial_vect;
    sqlite3 *db = NULL;
    char *err_msg = 0;

    OpenDB(dbpath, &db);
    sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, &err_msg);
    DropAllTables(db);

    /*
     * CPCA model data structure.
     * 
     * - **block_scores** matrix of scores
     * - **block_loadings** matrix of loadings
     * - **super_scores** matrix of super scores
     * - **super_weights** matrix of super weigths
     * - **scaling_factor** dvector of scaling factors
     * - **total_expvar** dvector of total explained variance
     * - **block_expvar** dvector list of block explained variance
     * - **colaverage** dvector list of column average
     * - **colscaling** dvector list of column scaling
     */
    write_vector_into_sqltable(db, "scaling_factor", cpca->scaling_factor);
    write_vector_into_sqltable(db, "total_expvar", cpca->total_expvar);
    
    initDVector(&serial_vect);
    serialize_tensor(cpca->block_scores, serial_vect);
    write_vector_into_sqltable(db, "block_scores", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(cpca->block_loadings, serial_vect);
    write_vector_into_sqltable(db, "block_loadings", serial_vect);
    DelDVector(&serial_vect);
    
    initDVector(&serial_vect);
    serialize_matrix(cpca->super_scores, serial_vect);
    write_vector_into_sqltable(db, "super_scores", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(cpca->super_weights, serial_vect);
    write_vector_into_sqltable(db, "super_weights", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_dvectorlist(cpca->block_expvar, serial_vect);
    write_vector_into_sqltable(db, "block_expvar", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_dvectorlist(cpca->colaverage, serial_vect);
    write_vector_into_sqltable(db, "colaverage", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_dvectorlist(cpca->colscaling, serial_vect);
    write_vector_into_sqltable(db, "colscaling", serial_vect);
    DelDVector(&serial_vect);

    sqlite3_exec(db, "COMMIT;", NULL, NULL, &err_msg);
    CloseDB(db);
}

void ReadCPCA(char *dbpath, CPCAMODEL *cpca)
{
    dvector *serial_vect;
    sqlite3* db = NULL;

    OpenDB(dbpath, &db);

    read_vector(db, "scaling_factor", cpca->scaling_factor);
    read_vector(db, "total_expvar", cpca->total_expvar);

    initDVector(&serial_vect);
    read_vector(db, "block_scores", serial_vect);
    deserialize_tensor(serial_vect, cpca->block_scores);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "block_loadings", serial_vect);
    deserialize_tensor(serial_vect, cpca->block_loadings);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "super_scores", serial_vect);
    deserialize_matrix(serial_vect, cpca->super_scores);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "super_weights", serial_vect);
    deserialize_matrix(serial_vect, cpca->super_weights);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "block_expvar", serial_vect);
    deserialize_dvectorlist(serial_vect, cpca->block_expvar);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "colaverage", serial_vect);
    deserialize_dvectorlist(serial_vect, cpca->colaverage);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "colscaling", serial_vect);
    deserialize_dvectorlist(serial_vect, cpca->colscaling);
    DelDVector(&serial_vect);

    CloseDB(db);
}


void WritePLS(char *dbpath, PLSMODEL *pls)
{
    dvector *serial_vect;
    sqlite3 *db = NULL;
    char *err_msg = 0;

    OpenDB(dbpath, &db);
    sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, &err_msg);
    DropAllTables(db);
    /*
     * PLS model data structure to write
     * 
     * - **xscores** x space scores
     * - **xloadings** x space loadings
     * - **xweights** x space weights
     * - **yscores** y space scores
     * - **yloadings** y space loadings
     * - **b** pls regression coefficients
     * - **xvarexp** variance explained in the x space
     * - **xcolaverage** x independent variable column average
     * - **xcolscaling** x independent variable column scaling
     * - **ycolaverage** y independent variable column average
     * - **ycolscaling** y independent variable column scaling
     * - **recalculated_y** y recalculated
     * - **recalc_residuals** y recalculated residuals
     * - **predicted_y** y predicted 
     * - **pred_residuals** y predicted residuals 
     * - **r2y_recalculated** r squared using y recalculated values
     * - **r2y_validation**
     * - **q2y** q squared using y predicted values
     * - **sdep** standard deviation over prediction using y predictions
     * - **sdec** standard deviation over recalculation using y recalculated
     * - **bias** bias
     * - **roc_recalculated** receiver operating characteristic using y recalculated 
     * - **roc_validation** receiver operating characteristic using y predicted 
     * - **roc_auc_recalculated** receiver operating characteristic area under the curve using y recalculated
     * - **roc_auc_validation** receiver operating characteristic area under the curve using y predicted
     * - **precision_recall_recalculated** precision-recall curve using y recalculated
     * - **precision_recall_validation** precision-recall curve using y predicted
     * - **precision_recall_ap_recalculated** precision-recall area under the curve using y recalculated
     * - **precision_recall_ap_validation** precision-recall area under the curve using y predicted
     * - **yscrambling** y-scrambling r-squared and q-squared
     */

    write_vector_into_sqltable(db, "xcolscaling", pls->xcolscaling);
    write_vector_into_sqltable(db, "xcolaverage", pls->xcolaverage);
    write_vector_into_sqltable(db, "ycolscaling", pls->ycolscaling);
    write_vector_into_sqltable(db, "ycolaverage", pls->ycolaverage);
    write_vector_into_sqltable(db, "xvarexp", pls->xvarexp);
    write_vector_into_sqltable(db, "b", pls->b);
    initDVector(&serial_vect);
    serialize_matrix(pls->xscores, serial_vect);
    write_vector_into_sqltable(db, "xscores", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(pls->xloadings, serial_vect);
    write_vector_into_sqltable(db, "xloadings", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(pls->xweights, serial_vect);
    write_vector_into_sqltable(db, "xweights", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(pls->yscores, serial_vect);
    write_vector_into_sqltable(db, "yscores", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(pls->yloadings, serial_vect);
    write_vector_into_sqltable(db, "yloadings", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(pls->recalculated_y, serial_vect);
    write_vector_into_sqltable(db, "recalculated_y", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(pls->recalc_residuals, serial_vect);
    write_vector_into_sqltable(db, "recalc_residuals", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(pls->predicted_y, serial_vect);
    write_vector_into_sqltable(db, "predicted_y", serial_vect);
    DelDVector(&serial_vect);
    
    initDVector(&serial_vect);
    serialize_matrix(pls->pred_residuals, serial_vect);
    write_vector_into_sqltable(db, "pred_residuals", serial_vect);
    DelDVector(&serial_vect);
    
    initDVector(&serial_vect);
    serialize_matrix(pls->r2y_validation, serial_vect);
    write_vector_into_sqltable(db, "r2y_validation", serial_vect);
    DelDVector(&serial_vect);
    
    initDVector(&serial_vect);
    serialize_matrix(pls->r2y_recalculated, serial_vect);
    write_vector_into_sqltable(db, "r2y_recalculated", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(pls->q2y, serial_vect);
    write_vector_into_sqltable(db, "q2y", serial_vect);
    DelDVector(&serial_vect);
    
    initDVector(&serial_vect);
    serialize_matrix(pls->sdep, serial_vect);
    write_vector_into_sqltable(db, "sdep", serial_vect);
    DelDVector(&serial_vect);
    
    initDVector(&serial_vect);
    serialize_matrix(pls->sdec, serial_vect);
    write_vector_into_sqltable(db, "sdec", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(pls->bias, serial_vect);
    write_vector_into_sqltable(db, "bias", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(pls->roc_auc_recalculated, serial_vect);
    write_vector_into_sqltable(db, "roc_auc_recalculated", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(pls->roc_auc_validation, serial_vect);
    write_vector_into_sqltable(db, "roc_auc_validation", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(pls->precision_recall_ap_recalculated, serial_vect);
    write_vector_into_sqltable(db, "precision_recall_ap_recalculated", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(pls->precision_recall_ap_validation, serial_vect);
    write_vector_into_sqltable(db, "precision_recall_ap_validation", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(pls->yscrambling, serial_vect);
    write_vector_into_sqltable(db, "yscrambling", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(pls->roc_recalculated, serial_vect);
    write_vector_into_sqltable(db, "roc_recalculated", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(pls->roc_validation, serial_vect);
    write_vector_into_sqltable(db, "roc_validation", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(pls->precision_recall_recalculated, serial_vect);
    write_vector_into_sqltable(db, "precision_recall_recalculated", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(pls->precision_recall_validation, serial_vect);
    write_vector_into_sqltable(db, "precision_recall_validation", serial_vect);
    DelDVector(&serial_vect);
    
    sqlite3_exec(db, "COMMIT;", NULL, NULL, &err_msg);
    CloseDB(db);
}

void ReadPLS(char *dbpath, PLSMODEL *pls)
{
    dvector *serial_vect;
    sqlite3* db = NULL;
    OpenDB(dbpath, &db);

    read_vector(db, "xcolscaling", pls->xcolscaling);
    read_vector(db, "xcolaverage", pls->xcolaverage);
    read_vector(db, "ycolaverage", pls->ycolaverage);
    read_vector(db, "ycolscaling", pls->ycolscaling);
    read_vector(db, "xvarexp", pls->xvarexp);
    read_vector(db, "b", pls->b);

    initDVector(&serial_vect);
    read_vector(db, "xscores", serial_vect);
    deserialize_matrix(serial_vect, pls->xscores);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "xloadings", serial_vect);
    deserialize_matrix(serial_vect, pls->xloadings);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "xweights", serial_vect);
    deserialize_matrix(serial_vect, pls->xweights);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "yscores", serial_vect);
    deserialize_matrix(serial_vect, pls->yscores);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "yloadings", serial_vect);
    deserialize_matrix(serial_vect, pls->yloadings);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "recalculated_y", serial_vect);
    deserialize_matrix(serial_vect, pls->recalculated_y);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "recalc_residuals", serial_vect);
    deserialize_matrix(serial_vect, pls->recalc_residuals);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "predicted_y", serial_vect);
    deserialize_matrix(serial_vect, pls->predicted_y);
    DelDVector(&serial_vect);
    
    initDVector(&serial_vect);
    read_vector(db, "pred_residuals", serial_vect);
    deserialize_matrix(serial_vect, pls->pred_residuals);
    DelDVector(&serial_vect);
    
    initDVector(&serial_vect);
    read_vector(db, "r2y_validation", serial_vect);
    deserialize_matrix(serial_vect, pls->r2y_validation);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "r2y_recalculated", serial_vect);
    deserialize_matrix(serial_vect, pls->r2y_recalculated);
    DelDVector(&serial_vect);
    
    initDVector(&serial_vect);
    read_vector(db, "q2y", serial_vect);
    deserialize_matrix(serial_vect, pls->q2y);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "sdep", serial_vect);
    deserialize_matrix(serial_vect, pls->sdep);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "sdec", serial_vect);
    deserialize_matrix(serial_vect, pls->sdec);
    DelDVector(&serial_vect);
    
    initDVector(&serial_vect);
    read_vector(db, "bias", serial_vect);
    deserialize_matrix(serial_vect, pls->bias);
    DelDVector(&serial_vect);
    
    initDVector(&serial_vect);
    read_vector(db, "roc_auc_recalculated", serial_vect);
    deserialize_matrix(serial_vect, pls->roc_auc_recalculated);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "roc_auc_validation", serial_vect);
    deserialize_matrix(serial_vect, pls->roc_auc_validation);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "precision_recall_ap_recalculated", serial_vect);
    deserialize_matrix(serial_vect, pls->precision_recall_ap_recalculated);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "precision_recall_ap_validation", serial_vect);
    deserialize_matrix(serial_vect, pls->precision_recall_ap_validation);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "yscrambling", serial_vect);
    deserialize_matrix(serial_vect, pls->yscrambling);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "roc_recalculated", serial_vect);
    deserialize_tensor(serial_vect, pls->roc_recalculated);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "roc_validation", serial_vect);
    deserialize_tensor(serial_vect, pls->roc_validation);
    DelDVector(&serial_vect);
    
    initDVector(&serial_vect);
    read_vector(db, "precision_recall_recalculated", serial_vect);
    deserialize_tensor(serial_vect, pls->precision_recall_recalculated);
    DelDVector(&serial_vect);
 
    initDVector(&serial_vect);
    read_vector(db, "precision_recall_validation", serial_vect);
    deserialize_tensor(serial_vect, pls->precision_recall_validation);
    DelDVector(&serial_vect);

    CloseDB(db);
}

void WriteMatrixCSV(char *path, matrix *m)
{
  FILE *f = fopen(path, "w");
  if(!f){
    fprintf(stderr, "Error: Unable to open file %s for writing\n", path);
    return;
  }
  size_t i, j;
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      fprintf(f, "%f%s", m->data[i][j], (j == m->col - 1) ? "" : ",");
    }
    fprintf(f, "\n");
  }
  fclose(f);
}

void WriteMLR(char *dbpath, MLRMODEL *mlr)
{
    dvector *serial_vect;
    sqlite3 *db = NULL;
    char *err_msg = 0;
    
    OpenDB(dbpath, &db);
    sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, &err_msg);
    DropAllTables(db);

    write_vector_into_sqltable(db, "ymean", mlr->ymean);
    write_vector_into_sqltable(db, "r2y_model", mlr->r2y_model);
    write_vector_into_sqltable(db, "q2y", mlr->q2y);
    write_vector_into_sqltable(db, "sdep", mlr->sdep);
    write_vector_into_sqltable(db, "sdec", mlr->sdec);
    write_vector_into_sqltable(db, "bias", mlr->bias);

    initDVector(&serial_vect);
    serialize_matrix(mlr->b, serial_vect);
    write_vector_into_sqltable(db, "b", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(mlr->recalculated_y, serial_vect);
    write_vector_into_sqltable(db, "recalculated_y", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(mlr->predicted_y, serial_vect);
    write_vector_into_sqltable(db, "predicted_y", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(mlr->recalc_residuals, serial_vect);
    write_vector_into_sqltable(db, "recalc_residuals", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(mlr->pred_residuals, serial_vect);
    write_vector_into_sqltable(db, "pred_residuals", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(mlr->r2q2scrambling, serial_vect);
    write_vector_into_sqltable(db, "r2q2scrambling", serial_vect);
    DelDVector(&serial_vect);

    sqlite3_exec(db, "COMMIT;", NULL, NULL, &err_msg);
    CloseDB(db);
}

void ReadMLR(char *dbpath, MLRMODEL *mlr)
{
    dvector *serial_vect;
    sqlite3* db = NULL;
    OpenDB(dbpath, &db);

    read_vector(db, "ymean", mlr->ymean);
    read_vector(db, "r2y_model", mlr->r2y_model);
    read_vector(db, "q2y", mlr->q2y);
    read_vector(db, "sdep", mlr->sdep);
    read_vector(db, "sdec", mlr->sdec);
    read_vector(db, "bias", mlr->bias);

    initDVector(&serial_vect);
    read_vector(db, "b", serial_vect);
    deserialize_matrix(serial_vect, mlr->b);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "recalculated_y", serial_vect);
    deserialize_matrix(serial_vect, mlr->recalculated_y);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "predicted_y", serial_vect);
    deserialize_matrix(serial_vect, mlr->predicted_y);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "recalc_residuals", serial_vect);
    deserialize_matrix(serial_vect, mlr->recalc_residuals);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "pred_residuals", serial_vect);
    deserialize_matrix(serial_vect, mlr->pred_residuals);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "r2q2scrambling", serial_vect);
    deserialize_matrix(serial_vect, mlr->r2q2scrambling);
    DelDVector(&serial_vect);

    CloseDB(db);
}

void WriteLDA(char *dbpath, LDAMODEL *lda)
{
    dvector *serial_vect;
    sqlite3 *db = NULL;
    char *err_msg = 0;
    
    OpenDB(dbpath, &db);
    sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, &err_msg);
    DropAllTables(db);

    write_vector_into_sqltable(db, "eval", lda->eval);
    write_vector_into_sqltable(db, "pprob", lda->pprob);
    write_vector_into_sqltable(db, "roc_aucs", lda->roc_aucs);
    write_vector_into_sqltable(db, "pr_aucs", lda->pr_aucs);

    /* Scalars stored as vectors of size 1 */
    dvector *scalars;
    NewDVector(&scalars, 2);
    scalars->data[0] = (double)lda->nclass;
    scalars->data[1] = (double)lda->class_start;
    write_vector_into_sqltable(db, "scalars", scalars);
    DelDVector(&scalars);

    initDVector(&serial_vect);
    serialize_uivector(lda->classid, serial_vect);
    write_vector_into_sqltable(db, "classid", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(lda->inv_cov, serial_vect);
    write_vector_into_sqltable(db, "inv_cov", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(lda->features, serial_vect);
    write_vector_into_sqltable(db, "features", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(lda->mnpdf, serial_vect);
    write_vector_into_sqltable(db, "mnpdf", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(lda->evect, serial_vect);
    write_vector_into_sqltable(db, "evect", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(lda->mu, serial_vect);
    write_vector_into_sqltable(db, "mu", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(lda->fmean, serial_vect);
    write_vector_into_sqltable(db, "fmean", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(lda->fsdev, serial_vect);
    write_vector_into_sqltable(db, "fsdev", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(lda->yscrambling, serial_vect);
    write_vector_into_sqltable(db, "yscrambling", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(lda->recalculated_y, serial_vect);
    write_vector_into_sqltable(db, "recalculated_y", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(lda->recalculated_residuals, serial_vect);
    write_vector_into_sqltable(db, "recalculated_residuals", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(lda->predicted_y, serial_vect);
    write_vector_into_sqltable(db, "predicted_y", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(lda->predicted_residuals, serial_vect);
    write_vector_into_sqltable(db, "predicted_residuals", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(lda->roc, serial_vect);
    write_vector_into_sqltable(db, "roc", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(lda->pr, serial_vect);
    write_vector_into_sqltable(db, "pr", serial_vect);
    DelDVector(&serial_vect);

    sqlite3_exec(db, "COMMIT;", NULL, NULL, &err_msg);
    CloseDB(db);
}

void ReadLDA(char *dbpath, LDAMODEL *lda)
{
    dvector *serial_vect;
    sqlite3* db = NULL;
    OpenDB(dbpath, &db);

    read_vector(db, "eval", lda->eval);
    read_vector(db, "pprob", lda->pprob);
    read_vector(db, "roc_aucs", lda->roc_aucs);
    read_vector(db, "pr_aucs", lda->pr_aucs);

    dvector *scalars;
    initDVector(&scalars);
    read_vector(db, "scalars", scalars);
    if(scalars->size >= 2) {
        lda->nclass = (size_t)scalars->data[0];
        lda->class_start = (size_t)scalars->data[1];
    }
    DelDVector(&scalars);

    initDVector(&serial_vect);
    read_vector(db, "classid", serial_vect);
    deserialize_uivector(serial_vect, lda->classid);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "inv_cov", serial_vect);
    deserialize_matrix(serial_vect, lda->inv_cov);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "features", serial_vect);
    deserialize_tensor(serial_vect, lda->features);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "mnpdf", serial_vect);
    deserialize_tensor(serial_vect, lda->mnpdf);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "evect", serial_vect);
    deserialize_matrix(serial_vect, lda->evect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "mu", serial_vect);
    deserialize_matrix(serial_vect, lda->mu);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "fmean", serial_vect);
    deserialize_matrix(serial_vect, lda->fmean);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "fsdev", serial_vect);
    deserialize_matrix(serial_vect, lda->fsdev);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "yscrambling", serial_vect);
    deserialize_matrix(serial_vect, lda->yscrambling);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "recalculated_y", serial_vect);
    deserialize_matrix(serial_vect, lda->recalculated_y);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "recalculated_residuals", serial_vect);
    deserialize_matrix(serial_vect, lda->recalculated_residuals);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "predicted_y", serial_vect);
    deserialize_matrix(serial_vect, lda->predicted_y);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "predicted_residuals", serial_vect);
    deserialize_matrix(serial_vect, lda->predicted_residuals);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "roc", serial_vect);
    deserialize_tensor(serial_vect, lda->roc);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "pr", serial_vect);
    deserialize_tensor(serial_vect, lda->pr);
    DelDVector(&serial_vect);

    CloseDB(db);
}

void WriteICA(char *dbpath, ICAMODEL *ica)
{
    dvector *serial_vect;
    sqlite3 *db = NULL;
    char *err_msg = 0;
    
    OpenDB(dbpath, &db);
    sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, &err_msg);
    DropAllTables(db);

    write_vector_into_sqltable(db, "colaverage", ica->colaverage);
    write_vector_into_sqltable(db, "colscaling", ica->colscaling);

    initDVector(&serial_vect);
    serialize_matrix(ica->S, serial_vect);
    write_vector_into_sqltable(db, "S", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(ica->W, serial_vect);
    write_vector_into_sqltable(db, "W", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(ica->whitening_matrix, serial_vect);
    write_vector_into_sqltable(db, "whitening_matrix", serial_vect);
    DelDVector(&serial_vect);

    sqlite3_exec(db, "COMMIT;", NULL, NULL, &err_msg);
    CloseDB(db);
}

void ReadICA(char *dbpath, ICAMODEL *ica)
{
    dvector *serial_vect;
    sqlite3* db = NULL;
    OpenDB(dbpath, &db);

    read_vector(db, "colaverage", ica->colaverage);
    read_vector(db, "colscaling", ica->colscaling);

    initDVector(&serial_vect);
    read_vector(db, "S", serial_vect);
    deserialize_matrix(serial_vect, ica->S);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "W", serial_vect);
    deserialize_matrix(serial_vect, ica->W);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "whitening_matrix", serial_vect);
    deserialize_matrix(serial_vect, ica->whitening_matrix);
    DelDVector(&serial_vect);

    CloseDB(db);
}

void WriteUPCA(char *dbpath, UPCAMODEL *upca)
{
    dvector *serial_vect;
    sqlite3 *db = NULL;
    char *err_msg = 0;
    
    OpenDB(dbpath, &db);
    sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, &err_msg);
    DropAllTables(db);

    write_vector_into_sqltable(db, "varexp", upca->varexp);

    initDVector(&serial_vect);
    serialize_matrix(upca->scores, serial_vect);
    write_vector_into_sqltable(db, "scores", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(upca->loadings, serial_vect);
    write_vector_into_sqltable(db, "loadings", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_dvectorlist(upca->colaverage, serial_vect);
    write_vector_into_sqltable(db, "colaverage", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_dvectorlist(upca->colscaling, serial_vect);
    write_vector_into_sqltable(db, "colscaling", serial_vect);
    DelDVector(&serial_vect);

    sqlite3_exec(db, "COMMIT;", NULL, NULL, &err_msg);
    CloseDB(db);
}

void ReadUPCA(char *dbpath, UPCAMODEL *upca)
{
    dvector *serial_vect;
    sqlite3* db = NULL;
    OpenDB(dbpath, &db);

    read_vector(db, "varexp", upca->varexp);

    initDVector(&serial_vect);
    read_vector(db, "scores", serial_vect);
    deserialize_matrix(serial_vect, upca->scores);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "loadings", serial_vect);
    deserialize_tensor(serial_vect, upca->loadings);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "colaverage", serial_vect);
    deserialize_dvectorlist(serial_vect, upca->colaverage);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "colscaling", serial_vect);
    deserialize_dvectorlist(serial_vect, upca->colscaling);
    DelDVector(&serial_vect);

    CloseDB(db);
}

void WriteUPLS(char *dbpath, UPLSMODEL *upls)
{
    dvector *serial_vect;
    sqlite3 *db = NULL;
    char *err_msg = 0;
    
    OpenDB(dbpath, &db);
    sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, &err_msg);
    DropAllTables(db);

    write_vector_into_sqltable(db, "xvarexp", upls->xvarexp);
    write_vector_into_sqltable(db, "b", upls->b);
    write_vector_into_sqltable(db, "r2x_model", upls->r2x_model);
    write_vector_into_sqltable(db, "r2x_validation", upls->r2x_validation);

    initDVector(&serial_vect);
    serialize_matrix(upls->xscores, serial_vect);
    write_vector_into_sqltable(db, "xscores", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(upls->xloadings, serial_vect);
    write_vector_into_sqltable(db, "xloadings", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(upls->xweights, serial_vect);
    write_vector_into_sqltable(db, "xweights", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_matrix(upls->yscores, serial_vect);
    write_vector_into_sqltable(db, "yscores", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(upls->yloadings, serial_vect);
    write_vector_into_sqltable(db, "yloadings", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_dvectorlist(upls->xcolaverage, serial_vect);
    write_vector_into_sqltable(db, "xcolaverage", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_dvectorlist(upls->xcolscaling, serial_vect);
    write_vector_into_sqltable(db, "xcolscaling", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_dvectorlist(upls->ycolaverage, serial_vect);
    write_vector_into_sqltable(db, "ycolaverage", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_dvectorlist(upls->ycolscaling, serial_vect);
    write_vector_into_sqltable(db, "ycolscaling", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(upls->r2y_model, serial_vect);
    write_vector_into_sqltable(db, "r2y_model", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(upls->r2y_validation, serial_vect);
    write_vector_into_sqltable(db, "r2y_validation", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(upls->q2y, serial_vect);
    write_vector_into_sqltable(db, "q2y", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(upls->sdep, serial_vect);
    write_vector_into_sqltable(db, "sdep", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(upls->sdec, serial_vect);
    write_vector_into_sqltable(db, "sdec", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(upls->recalculated_y, serial_vect);
    write_vector_into_sqltable(db, "recalculated_y", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(upls->predicted_y, serial_vect);
    write_vector_into_sqltable(db, "predicted_y", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(upls->recalc_residuals, serial_vect);
    write_vector_into_sqltable(db, "recalc_residuals", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(upls->pred_residuals, serial_vect);
    write_vector_into_sqltable(db, "pred_residuals", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(upls->q2y_yscrambling, serial_vect);
    write_vector_into_sqltable(db, "q2y_yscrambling", serial_vect);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    serialize_tensor(upls->sdep_yscrambling, serial_vect);
    write_vector_into_sqltable(db, "sdep_yscrambling", serial_vect);
    DelDVector(&serial_vect);

    sqlite3_exec(db, "COMMIT;", NULL, NULL, &err_msg);
    CloseDB(db);
}

void ReadUPLS(char *dbpath, UPLSMODEL *upls)
{
    dvector *serial_vect;
    sqlite3* db = NULL;
    OpenDB(dbpath, &db);

    read_vector(db, "xvarexp", upls->xvarexp);
    read_vector(db, "b", upls->b);
    read_vector(db, "r2x_model", upls->r2x_model);
    read_vector(db, "r2x_validation", upls->r2x_validation);

    initDVector(&serial_vect);
    read_vector(db, "xscores", serial_vect);
    deserialize_matrix(serial_vect, upls->xscores);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "xloadings", serial_vect);
    deserialize_tensor(serial_vect, upls->xloadings);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "xweights", serial_vect);
    deserialize_tensor(serial_vect, upls->xweights);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "yscores", serial_vect);
    deserialize_matrix(serial_vect, upls->yscores);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "yloadings", serial_vect);
    deserialize_tensor(serial_vect, upls->yloadings);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "xcolaverage", serial_vect);
    deserialize_dvectorlist(serial_vect, upls->xcolaverage);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "xcolscaling", serial_vect);
    deserialize_dvectorlist(serial_vect, upls->xcolscaling);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "ycolaverage", serial_vect);
    deserialize_dvectorlist(serial_vect, upls->ycolaverage);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "ycolscaling", serial_vect);
    deserialize_dvectorlist(serial_vect, upls->ycolscaling);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "r2y_model", serial_vect);
    deserialize_tensor(serial_vect, upls->r2y_model);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "r2y_validation", serial_vect);
    deserialize_tensor(serial_vect, upls->r2y_validation);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "q2y", serial_vect);
    deserialize_tensor(serial_vect, upls->q2y);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "sdep", serial_vect);
    deserialize_tensor(serial_vect, upls->sdep);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "sdec", serial_vect);
    deserialize_tensor(serial_vect, upls->sdec);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "recalculated_y", serial_vect);
    deserialize_tensor(serial_vect, upls->recalculated_y);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "predicted_y", serial_vect);
    deserialize_tensor(serial_vect, upls->predicted_y);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "recalc_residuals", serial_vect);
    deserialize_tensor(serial_vect, upls->recalc_residuals);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "pred_residuals", serial_vect);
    deserialize_tensor(serial_vect, upls->pred_residuals);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "q2y_yscrambling", serial_vect);
    deserialize_tensor(serial_vect, upls->q2y_yscrambling);
    DelDVector(&serial_vect);

    initDVector(&serial_vect);
    read_vector(db, "sdep_yscrambling", serial_vect);
    deserialize_tensor(serial_vect, upls->sdep_yscrambling);
    DelDVector(&serial_vect);

    CloseDB(db);
}