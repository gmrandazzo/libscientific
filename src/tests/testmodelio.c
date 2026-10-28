/* Unit tests for zero-dependency binary model container I/O (.lsm / .lsb).
 * Copyright (C) 2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
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
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "model_io.h"
#include "numeric.h"
#include "pca.h"
#include "pls.h"
#include "mlr.h"
#include "ica.h"

#define TEST_EPS 1e-12

static void assert_dvector_eq(const dvector *v1, const dvector *v2, const char *msg) {
    assert(v1->size == v2->size);
    for (size_t i = 0; i < v1->size; i++) {
        if (isnan(v1->data[i]) && isnan(v2->data[i])) continue;
        if (!FLOAT_EQ(v1->data[i], v2->data[i], TEST_EPS)) {
            fprintf(stderr, "dvector mismatch at %zu: %f != %f in %s\n", i, v1->data[i], v2->data[i], msg);
            assert(0);
        }
    }
}

static void assert_matrix_eq(const matrix *m1, const matrix *m2, const char *msg) {
    assert(m1->row == m2->row);
    assert(m1->col == m2->col);
    for (size_t i = 0; i < m1->row; i++) {
        for (size_t j = 0; j < m1->col; j++) {
            if (isnan(m1->data[i][j]) && isnan(m2->data[i][j])) continue;
            if (!FLOAT_EQ(m1->data[i][j], m2->data[i][j], TEST_EPS)) {
                fprintf(stderr, "matrix mismatch at [%zu][%zu]: %f != %f in %s\n", i, j, m1->data[i][j], m2->data[i][j], msg);
                assert(0);
            }
        }
    }
}

static void test_single_pca(void) {
    printf("Running test_single_pca (.lsm)...\n");

    matrix *X = NULL;
    NewMatrix(&X, 10, 4);
    for (size_t i = 0; i < 10; i++) {
        for (size_t j = 0; j < 4; j++) {
            X->data[i][j] = (double)(i + 1) * 1.5 + (double)(j + 1) * 0.3;
        }
    }

    PCAMODEL *pca = NULL;
    NewPCAModel(&pca);
    PCA(X, 1, 2, pca, NULL);

    const char *path = "test_pca_single.lsm";
    int st = SavePCAModel(path, pca);
    assert(st == LSCM_OK);

    PCAMODEL *loaded_pca = NULL;
    NewPCAModel(&loaded_pca);
    st = ImportPCAModel(path, loaded_pca);
    assert(st == LSCM_OK);

    assert_matrix_eq(pca->scores, loaded_pca->scores, "PCA scores");
    assert_matrix_eq(pca->loadings, loaded_pca->loadings, "PCA loadings");
    assert_dvector_eq(pca->colaverage, loaded_pca->colaverage, "PCA colaverage");
    assert_dvector_eq(pca->colscaling, loaded_pca->colscaling, "PCA colscaling");
    assert_dvector_eq(pca->varexp, loaded_pca->varexp, "PCA varexp");

    remove(path);
    DelPCAModel(&pca);
    DelPCAModel(&loaded_pca);
    DelMatrix(&X);
    printf("test_single_pca PASSED\n");
}

static void test_single_pls(void) {
    printf("Running test_single_pls (.lsm)...\n");

    matrix *X = NULL, *Y = NULL;
    NewMatrix(&X, 8, 3);
    NewMatrix(&Y, 8, 1);
    for (size_t i = 0; i < 8; i++) {
        X->data[i][0] = (double)(i + 1);
        X->data[i][1] = (double)((i * i) % 7 + 1);
        X->data[i][2] = (double)(i * 3 + 2);
        Y->data[i][0] = (double)(i + 1) * 3.7 + 1.2;
    }

    PLSMODEL *pls = NULL;
    NewPLSModel(&pls);
    PLS(X, Y, 2, 1, 0, pls, NULL);

    const char *path = "test_pls_single.lsm";
    int st = SavePLSModel(path, pls);
    assert(st == LSCM_OK);

    PLSMODEL *loaded_pls = NULL;
    NewPLSModel(&loaded_pls);
    st = ImportPLSModel(path, loaded_pls);
    assert(st == LSCM_OK);

    assert_matrix_eq(pls->xscores, loaded_pls->xscores, "PLS xscores");
    assert_matrix_eq(pls->xloadings, loaded_pls->xloadings, "PLS xloadings");
    assert_matrix_eq(pls->xweights, loaded_pls->xweights, "PLS xweights");
    assert_dvector_eq(pls->b, loaded_pls->b, "PLS regression vector b");
    assert_dvector_eq(pls->xcolaverage, loaded_pls->xcolaverage, "PLS xcolaverage");
    assert_dvector_eq(pls->xcolscaling, loaded_pls->xcolscaling, "PLS xcolscaling");

    remove(path);
    DelPLSModel(&pls);
    DelPLSModel(&loaded_pls);
    DelMatrix(&X);
    DelMatrix(&Y);
    printf("test_single_pls PASSED\n");
}

static void test_multi_model_bundle(void) {
    printf("Running test_multi_model_bundle (.lsb)...\n");

    matrix *X = NULL, *Y = NULL;
    NewMatrix(&X, 6, 2);
    NewMatrix(&Y, 6, 1);
    for (size_t i = 0; i < 6; i++) {
        X->data[i][0] = (double)(i + 1);
        X->data[i][1] = (double)((i * 3 + 1) % 5 + 1);
        Y->data[i][0] = (double)(i + 1) * 2.5;
    }

    PCAMODEL *pca = NULL;
    NewPCAModel(&pca);
    PCA(X, 1, 1, pca, NULL);

    PLSMODEL *pls = NULL;
    NewPLSModel(&pls);
    PLS(X, Y, 1, 1, 0, pls, NULL);

    MLRMODEL *mlr = NULL;
    NewMLRModel(&mlr);
    MLR(X, Y, mlr, NULL);

    MODEL_BUNDLE *bundle = NewModelBundle();
    ModelBundleAddModel(bundle, "pca_stage", MODEL_TYPE_PCA, pca);
    ModelBundleAddModel(bundle, "pls_stage", MODEL_TYPE_PLS, pls);
    ModelBundleAddModel(bundle, "mlr_stage", MODEL_TYPE_MLR, mlr);
    assert(bundle->count == 3);

    const char *path = "test_project.lsb";
    int st = SaveModelBundle(path, bundle);
    assert(st == LSCM_OK);

    MODEL_BUNDLE *loaded_bundle = ImportModelBundle(path);
    assert(loaded_bundle != NULL);
    assert(loaded_bundle->count == 3);

    ModelType t_pca = MODEL_TYPE_UNKNOWN;
    PCAMODEL *l_pca = (PCAMODEL*)ModelBundleGetModel(loaded_bundle, "pca_stage", &t_pca);
    assert(l_pca != NULL && t_pca == MODEL_TYPE_PCA);
    assert_matrix_eq(pca->scores, l_pca->scores, "Bundle PCA scores");

    ModelType t_pls = MODEL_TYPE_UNKNOWN;
    PLSMODEL *l_pls = (PLSMODEL*)ModelBundleGetModel(loaded_bundle, "pls_stage", &t_pls);
    assert(l_pls != NULL && t_pls == MODEL_TYPE_PLS);
    assert_dvector_eq(pls->b, l_pls->b, "Bundle PLS b");

    ModelType t_mlr = MODEL_TYPE_UNKNOWN;
    MLRMODEL *l_mlr = (MLRMODEL*)ModelBundleGetModel(loaded_bundle, "mlr_stage", &t_mlr);
    assert(l_mlr != NULL && t_mlr == MODEL_TYPE_MLR);
    assert_matrix_eq(mlr->b, l_mlr->b, "Bundle MLR b");

    remove(path);
    DelModelBundle(&bundle, 1);
    DelModelBundle(&loaded_bundle, 1);
    DelMatrix(&X);
    DelMatrix(&Y);
    printf("test_multi_model_bundle PASSED\n");
}

static void test_auto_detection(void) {
    printf("Running test_auto_detection...\n");

    matrix *X = NULL;
    NewMatrix(&X, 5, 2);
    for (size_t i = 0; i < 5; i++) {
        X->data[i][0] = (double)i;
        X->data[i][1] = (double)i * 3.0;
    }

    PCAMODEL *pca = NULL;
    NewPCAModel(&pca);
    PCA(X, 0, 1, pca, NULL);

    const char *single_path = "auto_detect_single.lsm";
    SavePCAModel(single_path, pca);

    MODEL_FILE_INFO info;
    int st = InspectModelFile(single_path, &info);
    assert(st == LSCM_OK);
    assert(info.is_single == true);
    assert(info.num_models == 1);
    assert(info.primary_type == MODEL_TYPE_PCA);

    void *single_model = NULL;
    ModelType single_type = MODEL_TYPE_UNKNOWN;
    MODEL_BUNDLE *bundle_out = NULL;

    int count = ImportModelAuto(single_path, &single_model, &single_type, &bundle_out);
    assert(count == 1);
    assert(single_model != NULL);
    assert(single_type == MODEL_TYPE_PCA);
    assert(bundle_out == NULL);

    PCAMODEL *loaded_pca = (PCAMODEL*)single_model;
    assert_matrix_eq(pca->scores, loaded_pca->scores, "Auto-detected PCA scores");
    DelPCAModel(&loaded_pca);
    remove(single_path);

    /* Test on bundle */
    MODEL_BUNDLE *b = NewModelBundle();
    ModelBundleAddModel(b, "p1", MODEL_TYPE_PCA, pca);
    const char *bundle_path = "auto_detect_bundle.lsb";
    SaveModelBundle(bundle_path, b);

    st = InspectModelFile(bundle_path, &info);
    assert(st == LSCM_OK);
    assert(info.is_single == true); // 1 model in bundle is single
    remove(bundle_path);

    /* 2 models in bundle */
    PCAMODEL *pca2 = NULL;
    NewPCAModel(&pca2);
    PCA(X, 1, 1, pca2, NULL);
    ModelBundleAddModel(b, "p2", MODEL_TYPE_PCA, pca2);

    SaveModelBundle(bundle_path, b);
    st = InspectModelFile(bundle_path, &info);
    assert(st == LSCM_OK);
    assert(info.is_single == false);
    assert(info.num_models == 2);

    single_model = NULL;
    single_type = MODEL_TYPE_UNKNOWN;
    bundle_out = NULL;
    count = ImportModelAuto(bundle_path, &single_model, &single_type, &bundle_out);
    assert(count == 2);
    assert(single_model == NULL);
    assert(bundle_out != NULL);
    assert(bundle_out->count == 2);

    DelModelBundle(&bundle_out, 1);
    DelModelBundle(&b, 1);
    remove(bundle_path);
    DelMatrix(&X);
    printf("test_auto_detection PASSED\n");
}

static void test_error_handling(void) {
    printf("Running test_error_handling...\n");

    MODEL_FILE_INFO info;
    int st = InspectModelFile("non_existent_file.lsm", &info);
    assert(st == LSCM_ERR_FILE_NOT_FOUND);

    /* Corrupted file with wrong magic */
    const char *corrupt_path = "corrupt_test.lsm";
    FILE *f = fopen(corrupt_path, "wb");
    uint8_t dummy_header[32];
    memset(dummy_header, 0, sizeof(dummy_header));
    uint32_t bad_magic = 0xDEADBEEF;
    memcpy(dummy_header, &bad_magic, sizeof(bad_magic));
    fwrite(dummy_header, 32, 1, f);
    fclose(f);

    st = InspectModelFile(corrupt_path, &info);
    assert(st == LSCM_ERR_INVALID_MAGIC);

    PCAMODEL *pca = NULL;
    NewPCAModel(&pca);
    st = ImportPCAModel(corrupt_path, pca);
    assert(st == LSCM_ERR_INVALID_MAGIC);

    remove(corrupt_path);
    DelPCAModel(&pca);
    printf("test_error_handling PASSED\n");
}

int main(void) {
    test_single_pca();
    test_single_pls();
    test_multi_model_bundle();
    test_auto_detection();
    test_error_handling();
    printf("\nAll Model I/O (.lsm / .lsb) tests PASSED successfully!\n");
    return 0;
}
