#include <stdio.h>
#include <sqlite3.h>

#include "vector.h"
#include "matrix.h"
#include "pca.h"
#include "pls.h"
#include "cpca.h"

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
    size_t i, j, c = 0;
    dvector *v;
    for(i = 0; i < dlst_serialized->size; i++){
        NewDVector(&v, (int)dlst_serialized->data[c++]);
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
    int lenght;
    char *zErrMsg = 0;
    char *sql;

    /* Create SQL statement */
    lenght = snprintf(NULL, 0, "CREATE TABLE IF NOT EXISTS %s (id INTEGER PRIMARY KEY AUTOINCREMENT, value REAL);", tabname);
    sql = malloc(lenght+1);
    snprintf(sql, lenght+1, "CREATE TABLE IF NOT EXISTS %s (id INTEGER PRIMARY KEY AUTOINCREMENT, value REAL);", tabname);
    #ifdef DEBUG
    printf("%s\n", sql);
    #endif

    /* Execute SQL statement */
    rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
    if(rc != SQLITE_OK){
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
    }
    #ifdef DEBUG
    else{
        fprintf(stdout, "Table created successfully\n");
    }
    #endif
    free(sql);

    for (i = 0; i < vect->size; i++){
        lenght = snprintf(NULL, 0, "INSERT INTO %s (value) VALUES (%.18f);", tabname, vect->data[i]);
        sql = malloc(lenght+1);
        snprintf(sql, lenght+1, "INSERT INTO %s (value) VALUES (%.18f);", tabname, vect->data[i]);
        #ifdef DEBUG
        printf("%s\n", sql);
        #endif
        sqlite3_stmt* stmt;
        rc = sqlite3_prepare_v2(db, sql, -1, &stmt, 0);
        if (rc != SQLITE_OK) {
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
            sqlite3_close(db);
            return;
        }
        sqlite3_bind_double(stmt, 1, vect->data[i]);
        rc = sqlite3_step(stmt);
        if (rc != SQLITE_DONE) {
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
        }
        sqlite3_finalize(stmt);
        free(sql);
    }
}

static void read_vector(sqlite3 *db, char *tabname, dvector *vect)
{
 /* SQL query to select the column data */
    int rc;
    char *sql;
    int lenght;
    sqlite3_stmt* stmt;

    lenght = snprintf(NULL, 0, "SELECT value FROM %s;", tabname);
    sql = malloc(lenght+1);
    snprintf(sql, lenght+1, "SELECT value FROM %s;", tabname);
    #ifdef DEBUG
    printf("%s\n", sql);
    #endif
    rc = sqlite3_prepare_v2(db, sql, -1, &stmt, 0);
    if (rc != SQLITE_OK){
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        return;
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
        return;
    }
}

static void OpenDB(char *dbpath, sqlite3 **db)
{
    int rc;
    /* Open the SQLite database */
    rc = sqlite3_open(dbpath, db);
    if (rc) {
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg((*db)));
        return;
    }
    #ifdef DEBUG
    else{
        fprintf(stdout, "Opened database successfully\n");
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
    OpenDB(dbpath, &db);
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
    OpenDB(dbpath, &db);
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
    OpenDB(dbpath, &db);

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
     * - **roc_auc_validation** eceiver operating characteristic area under the curve using y predicted
     * - **precision_recall_recalculated** precision-recall curve using y recalculated
     * - **precision_recall_validation** precision-recall curve using y predicted
     * - **precision_recall_ap_recalculated** precision-recall aread under the curve using y recalculated
     * - **precision_recall_ap_validation** precision-recall aread under the curve using y predicted
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
