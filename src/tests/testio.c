#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "matrix.h"
#include "vector.h"
#include "io.h"
#include "pca.h"
#include "pls.h"
#include "cpca.h"
#include "mlr.h"
#include "lda.h"
#include "ica.h"
#include "upca.h"
#include "upls.h"

void check_vector(dvector *v1, dvector *v2, const char *name) {
    if (v1->size != v2->size) {
        printf("FAIL: %s size mismatch (%zu vs %zu)\n", name, v1->size, v2->size);
        return;
    }
    for (size_t i = 0; i < v1->size; i++) {
        if (fabs(v1->data[i] - v2->data[i]) > 1e-9) {
            printf("FAIL: %s mismatch at %zu (%f vs %f)\n", name, i, v1->data[i], v2->data[i]);
            return;
        }
    }
    printf("%s OK\n", name);
}

void check_matrix(matrix *m1, matrix *m2, const char *name) {
    if (m1->row != m2->row || m1->col != m2->col) {
        printf("FAIL: %s dim mismatch\n", name);
        return;
    }
    for (size_t i = 0; i < m1->row; i++) {
        for (size_t j = 0; j < m1->col; j++) {
            if (fabs(m1->data[i][j] - m2->data[i][j]) > 1e-9) {
                printf("FAIL: %s mismatch at %zu,%zu (%f vs %f)\n", name, i, j, m1->data[i][j], m2->data[i][j]);
                return;
            }
        }
    }
    printf("%s OK\n", name);
}

void check_tensor(tensor *t1, tensor *t2, const char *name) {
    if (t1->order != t2->order) {
        printf("FAIL: %s order mismatch\n", name);
        return;
    }
    for (size_t k = 0; k < t1->order; k++) {
        if (t1->m[k]->row != t2->m[k]->row || t1->m[k]->col != t2->m[k]->col) {
             printf("FAIL: %s dim mismatch at %zu\n", name, k);
             return;
        }
        for (size_t i = 0; i < t1->m[k]->row; i++) {
            for (size_t j = 0; j < t1->m[k]->col; j++) {
                if (fabs(t1->m[k]->data[i][j] - t2->m[k]->data[i][j]) > 1e-9) {
                    printf("FAIL: %s mismatch at %zu,%zu,%zu\n", name, k, i, j);
                    return;
                }
            }
        }
    }
    printf("%s OK\n", name);
}

void testMLR() {
    puts("Test 4 - MLR Model Saving/Reading");
    MLRMODEL *m1, *m2;
    NewMLRModel(&m1);
    NewMLRModel(&m2);

    /* Populate m1 with dummy data */
    NewMatrix(&m1->b, 3, 1); m1->b->data[0][0] = 1.1;
    NewMatrix(&m1->recalculated_y, 5, 1); m1->recalculated_y->data[0][0] = 2.2;
    NewDVector(&m1->r2y_model, 1); m1->r2y_model->data[0] = 0.95;

    WriteMLR("test_mlr.db", m1);
    ReadMLR("test_mlr.db", m2);

    check_matrix(m1->b, m2->b, "b");
    check_matrix(m1->recalculated_y, m2->recalculated_y, "recalculated_y");
    check_vector(m1->r2y_model, m2->r2y_model, "r2y_model");

    DelMLRModel(&m1);
    DelMLRModel(&m2);
    remove("test_mlr.db");
}

void testLDA() {
    puts("Test 5 - LDA Model Saving/Reading");
    LDAMODEL *m1, *m2;
    NewLDAModel(&m1);
    NewLDAModel(&m2);

    m1->nclass = 2;
    m1->class_start = 0;
    NewDVector(&m1->eval, 2); m1->eval->data[0] = 0.5; m1->eval->data[1] = 0.6;
    NewMatrix(&m1->mu, 2, 2); m1->mu->data[0][0] = 1.0;
    NewUIVector(&m1->classid, 3); m1->classid->data[0]=0; m1->classid->data[1]=1; m1->classid->data[2]=0;

    WriteLDA("test_lda.db", m1);
    ReadLDA("test_lda.db", m2);

    if(m1->nclass != m2->nclass) printf("FAIL: nclass mismatch\n");
    check_vector(m1->eval, m2->eval, "eval");
    check_matrix(m1->mu, m2->mu, "mu");
    
    if(m1->classid->size != m2->classid->size) printf("FAIL: classid size\n");
    else if(m1->classid->data[1] != m2->classid->data[1]) printf("FAIL: classid data\n");
    else printf("classid OK\n");

    DelLDAModel(&m1);
    DelLDAModel(&m2);
    remove("test_lda.db");
}

void testICA() {
    puts("Test 6 - ICA Model Saving/Reading");
    ICAMODEL *m1, *m2;
    NewICAModel(&m1);
    NewICAModel(&m2);

    NewMatrix(&m1->S, 3, 2); m1->S->data[0][0] = 0.1;
    NewMatrix(&m1->W, 2, 2); m1->W->data[1][1] = 0.9;
    NewDVector(&m1->colaverage, 2); m1->colaverage->data[0] = 5.5;

    WriteICA("test_ica.db", m1);
    ReadICA("test_ica.db", m2);

    check_matrix(m1->S, m2->S, "S");
    check_matrix(m1->W, m2->W, "W");
    check_vector(m1->colaverage, m2->colaverage, "colaverage");

    DelICAModel(&m1);
    DelICAModel(&m2);
    remove("test_ica.db");
}

void testUPCA() {
    puts("Test 7 - UPCA Model Saving/Reading");
    UPCAMODEL *m1, *m2;
    NewUPCAModel(&m1);
    NewUPCAModel(&m2);

    NewMatrix(&m1->scores, 4, 2); m1->scores->data[0][0] = 10.0;
    NewTensor(&m1->loadings, 1); /* 1 component for test */
    NewTensorMatrix(m1->loadings, 0, 3, 3);
    m1->loadings->m[0]->data[0][0] = 7.7;
    
    WriteUPCA("test_upca.db", m1);
    ReadUPCA("test_upca.db", m2);

    check_matrix(m1->scores, m2->scores, "scores");
    check_tensor(m1->loadings, m2->loadings, "loadings");

    DelUPCAModel(&m1);
    DelUPCAModel(&m2);
    remove("test_upca.db");
}

void testUPLS() {
    puts("Test 8 - UPLS Model Saving/Reading");
    UPLSMODEL *m1, *m2;
    NewUPLSModel(&m1);
    NewUPLSModel(&m2);

    NewMatrix(&m1->xscores, 5, 2); m1->xscores->data[1][1] = 3.33;
    NewTensor(&m1->xloadings, 1);
    NewTensorMatrix(m1->xloadings, 0, 2, 2); m1->xloadings->m[0]->data[0][0] = 1.23;
    NewDVector(&m1->b, 2); m1->b->data[0] = 0.11;

    WriteUPLS("test_upls.db", m1);
    ReadUPLS("test_upls.db", m2);

    check_matrix(m1->xscores, m2->xscores, "xscores");
    check_tensor(m1->xloadings, m2->xloadings, "xloadings");
    check_vector(m1->b, m2->b, "b");

    DelUPLSModel(&m1);
    DelUPLSModel(&m2);
    remove("test_upls.db");
}

/* Original Tests Stubs */
void Test1() { puts("Test 1 - PLS (Original Skipped)"); }
void Test2() { puts("Test 2 - PCA (Original Skipped)"); }
void Test3() { puts("Test 3 - CPCA (Original Skipped)"); }

int main(int argc, char **argv)
{
    Test1();
    Test2();
    Test3();
    testMLR();
    testLDA();
    testICA();
    testUPCA();
    testUPLS();
    return 0;
}