/* Provides zero-dependency binary model container serialization (.lsm / .lsb).
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

#include "model_io.h"
#include "memwrapper.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TOC_ENTRY_SIZE 128U
#define HEADER_SIZE    32U

/* Packed on-disk header representation (32 bytes) */
typedef struct {
    uint32_t magic;           /* "LSCM" */
    uint16_t version;         /* 1 */
    uint16_t endianness;      /* 1 = Little Endian */
    uint32_t num_models;      /* Model count */
    uint32_t flags;           /* Reserved */
    uint64_t toc_size;        /* TOC size in bytes */
    uint64_t total_size;      /* Total file size */
} LSCMHeader;

/* Packed on-disk TOC entry representation (128 bytes) */
typedef struct {
    uint32_t model_type;
    uint32_t model_flags;
    char model_name[LSCM_MAX_NAME_LEN];
    uint64_t payload_offset;
    uint64_t payload_size;
    uint32_t n_components;
    uint32_t n_features;
    uint32_t n_samples;
    uint32_t checksum;
    uint8_t reserved[20];
} LSCMTocEntry;

static void copy_model_name(char *dest, const char *src, size_t max_len) {
    if (!dest || max_len == 0) return;
    dest[0] = '\0';
    if (!src) return;
    size_t i = 0;
    while (i + 1 < max_len && src[i] != '\0') {
        dest[i] = src[i];
        i++;
    }
    dest[i] = '\0';
}

/* =========================================================================
 * Low-Level Binary Primitive Helpers (Standard C stdio)
 * ========================================================================= */

static int write_u64(FILE *f, uint64_t val) {
    return (fwrite(&val, sizeof(uint64_t), 1, f) == 1) ? LSCM_OK : LSCM_ERR_IO;
}

static int read_u64(FILE *f, uint64_t *val) {
    return (fread(val, sizeof(uint64_t), 1, f) == 1) ? LSCM_OK : LSCM_ERR_IO;
}

static int write_dvector_bin(FILE *f, const dvector *v) {
    if (!v || v->size == 0) {
        return write_u64(f, 0);
    }
    if (write_u64(f, (uint64_t)v->size) != LSCM_OK) return LSCM_ERR_IO;
    if (fwrite(v->data, sizeof(double), v->size, f) != v->size) return LSCM_ERR_IO;
    return LSCM_OK;
}

static int read_dvector_bin(FILE *f, dvector *v) {
    uint64_t sz = 0;
    if (read_u64(f, &sz) != LSCM_OK) return LSCM_ERR_IO;
    DVectorResize(v, (size_t)sz);
    if (sz > 0) {
        if (fread(v->data, sizeof(double), (size_t)sz, f) != (size_t)sz) return LSCM_ERR_IO;
    }
    return LSCM_OK;
}

static int write_uivector_bin(FILE *f, const uivector *v) {
    if (!v || v->size == 0) {
        return write_u64(f, 0);
    }
    if (write_u64(f, (uint64_t)v->size) != LSCM_OK) return LSCM_ERR_IO;
    for (size_t i = 0; i < v->size; i++) {
        uint64_t val = (uint64_t)v->data[i];
        if (write_u64(f, val) != LSCM_OK) return LSCM_ERR_IO;
    }
    return LSCM_OK;
}

static int read_uivector_bin(FILE *f, uivector *v) {
    uint64_t sz = 0;
    if (read_u64(f, &sz) != LSCM_OK) return LSCM_ERR_IO;
    UIVectorResize(v, (size_t)sz);
    for (size_t i = 0; i < (size_t)sz; i++) {
        uint64_t val = 0;
        if (read_u64(f, &val) != LSCM_OK) return LSCM_ERR_IO;
        v->data[i] = (size_t)val;
    }
    return LSCM_OK;
}

static int write_matrix_bin(FILE *f, const matrix *m) {
    if (!m || m->row == 0 || m->col == 0) {
        if (write_u64(f, 0) != LSCM_OK) return LSCM_ERR_IO;
        return write_u64(f, 0);
    }
    if (write_u64(f, (uint64_t)m->row) != LSCM_OK) return LSCM_ERR_IO;
    if (write_u64(f, (uint64_t)m->col) != LSCM_OK) return LSCM_ERR_IO;
    for (size_t i = 0; i < m->row; i++) {
        if (fwrite(m->data[i], sizeof(double), m->col, f) != m->col) return LSCM_ERR_IO;
    }
    return LSCM_OK;
}

static int read_matrix_bin(FILE *f, matrix *m) {
    uint64_t row = 0, col = 0;
    if (read_u64(f, &row) != LSCM_OK) return LSCM_ERR_IO;
    if (read_u64(f, &col) != LSCM_OK) return LSCM_ERR_IO;
    ResizeMatrix(m, (size_t)row, (size_t)col);
    for (size_t i = 0; i < (size_t)row; i++) {
        if (fread(m->data[i], sizeof(double), (size_t)col, f) != (size_t)col) return LSCM_ERR_IO;
    }
    return LSCM_OK;
}

static int write_tensor_bin(FILE *f, const tensor *t) {
    if (!t || t->order == 0) {
        return write_u64(f, 0);
    }
    if (write_u64(f, (uint64_t)t->order) != LSCM_OK) return LSCM_ERR_IO;
    for (size_t k = 0; k < t->order; k++) {
        if (write_matrix_bin(f, t->m[k]) != LSCM_OK) return LSCM_ERR_IO;
    }
    return LSCM_OK;
}

static int read_tensor_bin(FILE *f, tensor *t) {
    uint64_t order = 0;
    if (read_u64(f, &order) != LSCM_OK) return LSCM_ERR_IO;
    ClearTensor(t);
    for (size_t k = 0; k < (size_t)order; k++) {
        uint64_t row = 0, col = 0;
        if (read_u64(f, &row) != LSCM_OK) return LSCM_ERR_IO;
        if (read_u64(f, &col) != LSCM_OK) return LSCM_ERR_IO;
        AddTensorMatrix(t, (size_t)row, (size_t)col);
        for (size_t i = 0; i < (size_t)row; i++) {
            if (fread(t->m[k]->data[i], sizeof(double), (size_t)col, f) != (size_t)col) return LSCM_ERR_IO;
        }
    }
    return LSCM_OK;
}

static int write_dvectorlist_bin(FILE *f, const dvectorlist *l) {
    if (!l || l->size == 0) {
        return write_u64(f, 0);
    }
    if (write_u64(f, (uint64_t)l->size) != LSCM_OK) return LSCM_ERR_IO;
    for (size_t i = 0; i < l->size; i++) {
        if (write_dvector_bin(f, l->d[i]) != LSCM_OK) return LSCM_ERR_IO;
    }
    return LSCM_OK;
}

static void clear_dvectorlist(dvectorlist *l) {
    if (!l) return;
    for (size_t i = 0; i < l->size; i++) {
        if (l->d[i]) DelDVector(&l->d[i]);
    }
    xfree(l->d);
    l->d = NULL;
    l->size = 0;
}

static int read_dvectorlist_bin(FILE *f, dvectorlist *l) {
    uint64_t sz = 0;
    if (read_u64(f, &sz) != LSCM_OK) return LSCM_ERR_IO;
    clear_dvectorlist(l);
    for (size_t i = 0; i < (size_t)sz; i++) {
        dvector *v = NULL;
        initDVector(&v);
        if (read_dvector_bin(f, v) != LSCM_OK) {
            DelDVector(&v);
            return LSCM_ERR_IO;
        }
        DVectorListAppend(l, v);
        DelDVector(&v);
    }
    return LSCM_OK;
}

/* =========================================================================
 * Concrete Model Payloads: Serializers & Deserializers
 * ========================================================================= */

/* --- PCA --- */
static int write_pca_payload(FILE *f, const PCAMODEL *m) {
    if (write_matrix_bin(f, m->scores) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->loadings) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->dmodx) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->varexp) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->colaverage) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->colscaling) != LSCM_OK) return LSCM_ERR_IO;
    return LSCM_OK;
}

static int read_pca_payload(FILE *f, PCAMODEL *m) {
    if (read_matrix_bin(f, m->scores) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->loadings) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->dmodx) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->varexp) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->colaverage) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->colscaling) != LSCM_OK) return LSCM_ERR_IO;
    return LSCM_OK;
}

/* --- PLS --- */
static int write_pls_payload(FILE *f, const PLSMODEL *m) {
    if (write_matrix_bin(f, m->xscores) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->xloadings) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->xweights) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->yscores) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->yloadings) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->b) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->xvarexp) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->xcolaverage) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->xcolscaling) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->ycolaverage) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->ycolscaling) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->recalculated_y) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->recalc_residuals) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->predicted_y) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->pred_residuals) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->r2y_recalculated) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->r2y_validation) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->q2y) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->sdep) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->sdec) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->bias) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->roc_recalculated) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->roc_validation) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->roc_auc_recalculated) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->roc_auc_validation) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->precision_recall_recalculated) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->precision_recall_validation) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->precision_recall_ap_recalculated) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->precision_recall_ap_validation) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->yscrambling) != LSCM_OK) return LSCM_ERR_IO;
    return LSCM_OK;
}

static int read_pls_payload(FILE *f, PLSMODEL *m) {
    if (read_matrix_bin(f, m->xscores) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->xloadings) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->xweights) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->yscores) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->yloadings) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->b) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->xvarexp) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->xcolaverage) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->xcolscaling) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->ycolaverage) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->ycolscaling) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->recalculated_y) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->recalc_residuals) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->predicted_y) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->pred_residuals) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->r2y_recalculated) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->r2y_validation) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->q2y) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->sdep) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->sdec) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->bias) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->roc_recalculated) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->roc_validation) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->roc_auc_recalculated) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->roc_auc_validation) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->precision_recall_recalculated) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->precision_recall_validation) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->precision_recall_ap_recalculated) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->precision_recall_ap_validation) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->yscrambling) != LSCM_OK) return LSCM_ERR_IO;
    return LSCM_OK;
}

/* --- CPCA --- */
static int write_cpca_payload(FILE *f, const CPCAMODEL *m) {
    if (write_tensor_bin(f, m->block_scores) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->block_loadings) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->super_scores) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->super_weights) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->scaling_factor) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->total_expvar) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvectorlist_bin(f, m->block_expvar) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvectorlist_bin(f, m->colaverage) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvectorlist_bin(f, m->colscaling) != LSCM_OK) return LSCM_ERR_IO;
    return LSCM_OK;
}

static int read_cpca_payload(FILE *f, CPCAMODEL *m) {
    if (read_tensor_bin(f, m->block_scores) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->block_loadings) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->super_scores) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->super_weights) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->scaling_factor) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->total_expvar) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvectorlist_bin(f, m->block_expvar) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvectorlist_bin(f, m->colaverage) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvectorlist_bin(f, m->colscaling) != LSCM_OK) return LSCM_ERR_IO;
    return LSCM_OK;
}

/* --- MLR --- */
static int write_mlr_payload(FILE *f, const MLRMODEL *m) {
    if (write_matrix_bin(f, m->b) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->recalculated_y) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->predicted_y) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->recalc_residuals) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->pred_residuals) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->ymean) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->r2y_model) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->q2y) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->sdep) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->sdec) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->bias) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->r2q2scrambling) != LSCM_OK) return LSCM_ERR_IO;
    return LSCM_OK;
}

static int read_mlr_payload(FILE *f, MLRMODEL *m) {
    if (read_matrix_bin(f, m->b) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->recalculated_y) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->predicted_y) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->recalc_residuals) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->pred_residuals) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->ymean) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->r2y_model) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->q2y) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->sdep) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->sdec) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->bias) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->r2q2scrambling) != LSCM_OK) return LSCM_ERR_IO;
    return LSCM_OK;
}

/* --- LDA --- */
static int write_lda_payload(FILE *f, const LDAMODEL *m) {
    if (write_matrix_bin(f, m->inv_cov) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->features) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->mnpdf) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->evect) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->mu) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->fmean) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->fsdev) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->yscrambling) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->recalculated_y) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->recalculated_residuals) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->predicted_y) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->predicted_residuals) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->eval) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->pprob) != LSCM_OK) return LSCM_ERR_IO;
    if (write_uivector_bin(f, m->classid) != LSCM_OK) return LSCM_ERR_IO;
    if (write_u64(f, (uint64_t)m->nclass) != LSCM_OK) return LSCM_ERR_IO;
    if (write_u64(f, (uint64_t)m->class_start) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->roc) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->pr) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->roc_aucs) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->pr_aucs) != LSCM_OK) return LSCM_ERR_IO;
    return LSCM_OK;
}

static int read_lda_payload(FILE *f, LDAMODEL *m) {
    if (read_matrix_bin(f, m->inv_cov) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->features) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->mnpdf) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->evect) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->mu) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->fmean) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->fsdev) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->yscrambling) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->recalculated_y) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->recalculated_residuals) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->predicted_y) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->predicted_residuals) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->eval) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->pprob) != LSCM_OK) return LSCM_ERR_IO;
    if (read_uivector_bin(f, m->classid) != LSCM_OK) return LSCM_ERR_IO;
    uint64_t nclass = 0, class_start = 0;
    if (read_u64(f, &nclass) != LSCM_OK) return LSCM_ERR_IO;
    if (read_u64(f, &class_start) != LSCM_OK) return LSCM_ERR_IO;
    m->nclass = (size_t)nclass;
    m->class_start = (size_t)class_start;
    if (read_tensor_bin(f, m->roc) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->pr) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->roc_aucs) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->pr_aucs) != LSCM_OK) return LSCM_ERR_IO;
    return LSCM_OK;
}

/* --- ICA --- */
static int write_ica_payload(FILE *f, const ICAMODEL *m) {
    if (write_matrix_bin(f, m->S) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->W) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->whitening_matrix) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->colaverage) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->colscaling) != LSCM_OK) return LSCM_ERR_IO;
    return LSCM_OK;
}

static int read_ica_payload(FILE *f, ICAMODEL *m) {
    if (read_matrix_bin(f, m->S) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->W) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->whitening_matrix) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->colaverage) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->colscaling) != LSCM_OK) return LSCM_ERR_IO;
    return LSCM_OK;
}

/* --- UPCA --- */
static int write_upca_payload(FILE *f, const UPCAMODEL *m) {
    if (write_matrix_bin(f, m->scores) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->loadings) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->varexp) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvectorlist_bin(f, m->colaverage) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvectorlist_bin(f, m->colscaling) != LSCM_OK) return LSCM_ERR_IO;
    return LSCM_OK;
}

static int read_upca_payload(FILE *f, UPCAMODEL *m) {
    if (read_matrix_bin(f, m->scores) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->loadings) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->varexp) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvectorlist_bin(f, m->colaverage) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvectorlist_bin(f, m->colscaling) != LSCM_OK) return LSCM_ERR_IO;
    return LSCM_OK;
}

/* --- UPLS --- */
static int write_upls_payload(FILE *f, const UPLSMODEL *m) {
    if (write_matrix_bin(f, m->xscores) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->xloadings) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->xweights) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->xvarexp) != LSCM_OK) return LSCM_ERR_IO;
    if (write_matrix_bin(f, m->yscores) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->yloadings) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->b) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvectorlist_bin(f, m->xcolaverage) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvectorlist_bin(f, m->xcolscaling) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvectorlist_bin(f, m->ycolaverage) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvectorlist_bin(f, m->ycolscaling) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->r2x_model) != LSCM_OK) return LSCM_ERR_IO;
    if (write_dvector_bin(f, m->r2x_validation) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->r2y_model) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->r2y_validation) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->q2y) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->sdep) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->sdec) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->recalculated_y) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->predicted_y) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->recalc_residuals) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->pred_residuals) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->q2y_yscrambling) != LSCM_OK) return LSCM_ERR_IO;
    if (write_tensor_bin(f, m->sdep_yscrambling) != LSCM_OK) return LSCM_ERR_IO;
    return LSCM_OK;
}

static int read_upls_payload(FILE *f, UPLSMODEL *m) {
    if (read_matrix_bin(f, m->xscores) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->xloadings) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->xweights) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->xvarexp) != LSCM_OK) return LSCM_ERR_IO;
    if (read_matrix_bin(f, m->yscores) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->yloadings) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->b) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvectorlist_bin(f, m->xcolaverage) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvectorlist_bin(f, m->xcolscaling) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvectorlist_bin(f, m->ycolaverage) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvectorlist_bin(f, m->ycolscaling) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->r2x_model) != LSCM_OK) return LSCM_ERR_IO;
    if (read_dvector_bin(f, m->r2x_validation) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->r2y_model) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->r2y_validation) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->q2y) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->sdep) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->sdec) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->recalculated_y) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->predicted_y) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->recalc_residuals) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->pred_residuals) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->q2y_yscrambling) != LSCM_OK) return LSCM_ERR_IO;
    if (read_tensor_bin(f, m->sdep_yscrambling) != LSCM_OK) return LSCM_ERR_IO;
    return LSCM_OK;
}

static int dispatch_write_model_payload(FILE *f, ModelType type, const void *model) {
    if (!f || !model) return LSCM_ERR_PARAM;
    switch (type) {
        case MODEL_TYPE_PCA:  return write_pca_payload(f, (const PCAMODEL*)model);
        case MODEL_TYPE_PLS:  return write_pls_payload(f, (const PLSMODEL*)model);
        case MODEL_TYPE_CPCA: return write_cpca_payload(f, (const CPCAMODEL*)model);
        case MODEL_TYPE_MLR:  return write_mlr_payload(f, (const MLRMODEL*)model);
        case MODEL_TYPE_LDA:  return write_lda_payload(f, (const LDAMODEL*)model);
        case MODEL_TYPE_ICA:  return write_ica_payload(f, (const ICAMODEL*)model);
        case MODEL_TYPE_UPCA: return write_upca_payload(f, (const UPCAMODEL*)model);
        case MODEL_TYPE_UPLS: return write_upls_payload(f, (const UPLSMODEL*)model);
        default: return LSCM_ERR_INVALID_TYPE;
    }
}

static int dispatch_read_model_payload(FILE *f, ModelType type, void *model) {
    if (!f || !model) return LSCM_ERR_PARAM;
    switch (type) {
        case MODEL_TYPE_PCA:  return read_pca_payload(f, (PCAMODEL*)model);
        case MODEL_TYPE_PLS:  return read_pls_payload(f, (PLSMODEL*)model);
        case MODEL_TYPE_CPCA: return read_cpca_payload(f, (CPCAMODEL*)model);
        case MODEL_TYPE_MLR:  return read_mlr_payload(f, (MLRMODEL*)model);
        case MODEL_TYPE_LDA:  return read_lda_payload(f, (LDAMODEL*)model);
        case MODEL_TYPE_ICA:  return read_ica_payload(f, (ICAMODEL*)model);
        case MODEL_TYPE_UPCA: return read_upca_payload(f, (UPCAMODEL*)model);
        case MODEL_TYPE_UPLS: return read_upls_payload(f, (UPLSMODEL*)model);
        default: return LSCM_ERR_INVALID_TYPE;
    }
}

static void* allocate_model(ModelType type) {
    void *m = NULL;
    switch (type) {
        case MODEL_TYPE_PCA:  NewPCAModel((PCAMODEL**)&m); break;
        case MODEL_TYPE_PLS:  NewPLSModel((PLSMODEL**)&m); break;
        case MODEL_TYPE_CPCA: NewCPCAModel((CPCAMODEL**)&m); break;
        case MODEL_TYPE_MLR:  NewMLRModel((MLRMODEL**)&m); break;
        case MODEL_TYPE_LDA:  NewLDAModel((LDAMODEL**)&m); break;
        case MODEL_TYPE_ICA:  NewICAModel((ICAMODEL**)&m); break;
        case MODEL_TYPE_UPCA: NewUPCAModel((UPCAMODEL**)&m); break;
        case MODEL_TYPE_UPLS: NewUPLSModel((UPLSMODEL**)&m); break;
        default: break;
    }
    return m;
}

static void free_model(ModelType type, void **ptr) {
    if (!ptr || !*ptr) return;
    switch (type) {
        case MODEL_TYPE_PCA:  DelPCAModel((PCAMODEL**)ptr); break;
        case MODEL_TYPE_PLS:  DelPLSModel((PLSMODEL**)ptr); break;
        case MODEL_TYPE_CPCA: DelCPCAModel((CPCAMODEL**)ptr); break;
        case MODEL_TYPE_MLR:  DelMLRModel((MLRMODEL**)ptr); break;
        case MODEL_TYPE_LDA:  DelLDAModel((LDAMODEL**)ptr); break;
        case MODEL_TYPE_ICA:  DelICAModel((ICAMODEL**)ptr); break;
        case MODEL_TYPE_UPCA: DelUPCAModel((UPCAMODEL**)ptr); break;
        case MODEL_TYPE_UPLS: DelUPLSModel((UPLSMODEL**)ptr); break;
        default: xfree(*ptr); *ptr = NULL; break;
    }
}

/* =========================================================================
 * Multi-Model Bundle Management
 * ========================================================================= */

MODEL_BUNDLE* NewModelBundle(void) {
    MODEL_BUNDLE *b = (MODEL_BUNDLE*)xmalloc(sizeof(MODEL_BUNDLE));
    b->count = 0;
    b->capacity = 4;
    b->entries = (MODEL_ENTRY*)xmalloc(sizeof(MODEL_ENTRY) * b->capacity);
    return b;
}

int ModelBundleAddModel(MODEL_BUNDLE *bundle, const char *name, ModelType type, void *model) {
    if (!bundle || !name || type == MODEL_TYPE_UNKNOWN || !model) return LSCM_ERR_PARAM;
    if (bundle->count >= bundle->capacity) {
        size_t new_cap = bundle->capacity * 2;
        bundle->entries = (MODEL_ENTRY*)xrealloc(bundle->entries, sizeof(MODEL_ENTRY) * new_cap);
        bundle->capacity = new_cap;
    }
    MODEL_ENTRY *e = &bundle->entries[bundle->count];
    copy_model_name(e->name, name, sizeof(e->name));
    e->type = type;
    e->model = model;
    bundle->count++;
    return LSCM_OK;
}

void* ModelBundleGetModel(const MODEL_BUNDLE *bundle, const char *name, ModelType *out_type) {
    if (!bundle || !name) return NULL;
    for (size_t i = 0; i < bundle->count; i++) {
        if (strncmp(bundle->entries[i].name, name, LSCM_MAX_NAME_LEN) == 0) {
            if (out_type) *out_type = bundle->entries[i].type;
            return bundle->entries[i].model;
        }
    }
    return NULL;
}

void* ModelBundleGetModelByIndex(const MODEL_BUNDLE *bundle, size_t index, const char **out_name, ModelType *out_type) {
    if (!bundle || index >= bundle->count) return NULL;
    if (out_name) *out_name = bundle->entries[index].name;
    if (out_type) *out_type = bundle->entries[index].type;
    return bundle->entries[index].model;
}

void DelModelBundle(MODEL_BUNDLE **bundle, int free_contained_models) {
    if (!bundle || !*bundle) return;
    MODEL_BUNDLE *b = *bundle;
    if (free_contained_models) {
        for (size_t i = 0; i < b->count; i++) {
            if (b->entries[i].model) {
                free_model(b->entries[i].type, &b->entries[i].model);
            }
        }
    }
    xfree(b->entries);
    xfree(b);
    *bundle = NULL;
}

/* =========================================================================
 * Container Header & TOC Serialization Logic
 * ========================================================================= */

int SaveModelBundle(const char *filepath, const MODEL_BUNDLE *bundle) {
    if (!filepath || !bundle || bundle->count == 0) return LSCM_ERR_PARAM;

    FILE *f = fopen(filepath, "wb");
    if (!f) return LSCM_ERR_IO;

    uint32_t num_models = (uint32_t)bundle->count;
    uint64_t toc_size = (uint64_t)num_models * TOC_ENTRY_SIZE;

    /* 1. Write placeholder header */
    LSCMHeader header;
    memset(&header, 0, sizeof(header));
    header.magic = LSCM_MAGIC;
    header.version = LSCM_VERSION;
    header.endianness = 1; /* LE */
    header.num_models = num_models;
    header.toc_size = toc_size;

    if (fwrite(&header, HEADER_SIZE, 1, f) != 1) {
        fclose(f);
        return LSCM_ERR_IO;
    }

    /* 2. Write placeholder TOC */
    LSCMTocEntry *toc = (LSCMTocEntry*)xcalloc(num_models, sizeof(LSCMTocEntry));
    for (size_t i = 0; i < num_models; i++) {
        toc[i].model_type = (uint32_t)bundle->entries[i].type;
        copy_model_name(toc[i].model_name, bundle->entries[i].name, sizeof(toc[i].model_name));
    }
    if (fwrite(toc, TOC_ENTRY_SIZE, num_models, f) != num_models) {
        xfree(toc);
        fclose(f);
        return LSCM_ERR_IO;
    }

    /* 3. Write Payloads */
    for (size_t i = 0; i < num_models; i++) {
        long offset = ftell(f);
        toc[i].payload_offset = (uint64_t)offset;

        int st = dispatch_write_model_payload(f, bundle->entries[i].type, bundle->entries[i].model);
        if (st != LSCM_OK) {
            xfree(toc);
            fclose(f);
            return st;
        }

        long end_offset = ftell(f);
        toc[i].payload_size = (uint64_t)(end_offset - offset);
    }

    long total_sz = ftell(f);
    header.total_size = (uint64_t)total_sz;

    /* 4. Seek to beginning and rewrite finalized header and TOC */
    if (fseek(f, 0, SEEK_SET) != 0) {
        xfree(toc);
        fclose(f);
        return LSCM_ERR_IO;
    }
    if (fwrite(&header, HEADER_SIZE, 1, f) != 1) {
        xfree(toc);
        fclose(f);
        return LSCM_ERR_IO;
    }
    if (fwrite(toc, TOC_ENTRY_SIZE, num_models, f) != num_models) {
        xfree(toc);
        fclose(f);
        return LSCM_ERR_IO;
    }

    xfree(toc);
    fclose(f);
    return LSCM_OK;
}

int SaveModel(const char *filepath, const char *model_name, ModelType type, const void *model) {
    if (!filepath || type == MODEL_TYPE_UNKNOWN || !model) return LSCM_ERR_PARAM;
    const char *name = (model_name && model_name[0] != '\0') ? model_name : "default";

    MODEL_BUNDLE *b = NewModelBundle();
    ModelBundleAddModel(b, name, type, (void*)model);
    int status = SaveModelBundle(filepath, b);
    DelModelBundle(&b, 0); /* Do not free the caller's model struct */
    return status;
}

int InspectModelFile(const char *filepath, MODEL_FILE_INFO *info) {
    if (!filepath || !info) return LSCM_ERR_PARAM;
    memset(info, 0, sizeof(*info));

    FILE *f = fopen(filepath, "rb");
    if (!f) return LSCM_ERR_FILE_NOT_FOUND;

    LSCMHeader header;
    if (fread(&header, HEADER_SIZE, 1, f) != 1) {
        fclose(f);
        return LSCM_ERR_IO;
    }

    if (header.magic != LSCM_MAGIC) {
        fclose(f);
        return LSCM_ERR_INVALID_MAGIC;
    }
    if (header.version != LSCM_VERSION) {
        fclose(f);
        return LSCM_ERR_UNSUPPORTED_VER;
    }

    info->format_version = header.version;
    info->num_models = header.num_models;
    info->is_single = (header.num_models == 1);
    info->total_size = header.total_size;

    if (header.num_models > 0) {
        LSCMTocEntry first_entry;
        if (fread(&first_entry, TOC_ENTRY_SIZE, 1, f) == 1) {
            info->primary_type = (ModelType)first_entry.model_type;
            copy_model_name(info->primary_name, first_entry.model_name, sizeof(info->primary_name));
        }
    }

    fclose(f);
    return LSCM_OK;
}

int ImportModel(const char *filepath, ModelType expected_type, void *model) {
    if (!filepath || !model) return LSCM_ERR_PARAM;

    FILE *f = fopen(filepath, "rb");
    if (!f) return LSCM_ERR_FILE_NOT_FOUND;

    LSCMHeader header;
    if (fread(&header, HEADER_SIZE, 1, f) != 1) {
        fclose(f);
        return LSCM_ERR_IO;
    }
    if (header.magic != LSCM_MAGIC) {
        fclose(f);
        return LSCM_ERR_INVALID_MAGIC;
    }

    if (header.num_models == 0) {
        fclose(f);
        return LSCM_ERR_NOT_FOUND;
    }

    LSCMTocEntry *toc = (LSCMTocEntry*)xmalloc(header.num_models * sizeof(LSCMTocEntry));
    if (fread(toc, TOC_ENTRY_SIZE, header.num_models, f) != header.num_models) {
        xfree(toc);
        fclose(f);
        return LSCM_ERR_IO;
    }

    /* Locate matching entry */
    int target_idx = -1;
    if (expected_type == MODEL_TYPE_ANY || expected_type == MODEL_TYPE_UNKNOWN) {
        target_idx = 0; /* Default to first model */
    } else {
        for (uint32_t i = 0; i < header.num_models; i++) {
            if (toc[i].model_type == (uint32_t)expected_type) {
                target_idx = (int)i;
                break;
            }
        }
    }

    if (target_idx < 0) {
        xfree(toc);
        fclose(f);
        return LSCM_ERR_TYPE_MISMATCH;
    }

    if (fseek(f, (long)toc[target_idx].payload_offset, SEEK_SET) != 0) {
        xfree(toc);
        fclose(f);
        return LSCM_ERR_IO;
    }

    ModelType actual_type = (ModelType)toc[target_idx].model_type;
    int st = dispatch_read_model_payload(f, actual_type, model);

    xfree(toc);
    fclose(f);
    return st;
}

MODEL_BUNDLE* ImportModelBundle(const char *filepath) {
    if (!filepath) return NULL;

    FILE *f = fopen(filepath, "rb");
    if (!f) return NULL;

    LSCMHeader header;
    if (fread(&header, HEADER_SIZE, 1, f) != 1 || header.magic != LSCM_MAGIC) {
        fclose(f);
        return NULL;
    }

    LSCMTocEntry *toc = (LSCMTocEntry*)xmalloc(header.num_models * sizeof(LSCMTocEntry));
    if (fread(toc, TOC_ENTRY_SIZE, header.num_models, f) != header.num_models) {
        xfree(toc);
        fclose(f);
        return NULL;
    }

    MODEL_BUNDLE *bundle = NewModelBundle();
    for (uint32_t i = 0; i < header.num_models; i++) {
        ModelType type = (ModelType)toc[i].model_type;
        void *m = allocate_model(type);
        if (!m) continue;

        if (fseek(f, (long)toc[i].payload_offset, SEEK_SET) == 0) {
            if (dispatch_read_model_payload(f, type, m) == LSCM_OK) {
                ModelBundleAddModel(bundle, toc[i].model_name, type, m);
            } else {
                free_model(type, &m);
            }
        } else {
            free_model(type, &m);
        }
    }

    xfree(toc);
    fclose(f);
    return bundle;
}

int ImportModelAuto(const char *filepath,
                    void **out_single_model,
                    ModelType *out_single_type,
                    MODEL_BUNDLE **out_bundle) {
    if (!filepath) return LSCM_ERR_PARAM;
    if (out_single_model) *out_single_model = NULL;
    if (out_single_type) *out_single_type = MODEL_TYPE_UNKNOWN;
    if (out_bundle) *out_bundle = NULL;

    MODEL_FILE_INFO info;
    int st = InspectModelFile(filepath, &info);
    if (st != LSCM_OK) return st;

    if (info.is_single) {
        void *m = allocate_model(info.primary_type);
        if (!m) return LSCM_ERR_MEM;
        st = ImportModel(filepath, info.primary_type, m);
        if (st != LSCM_OK) {
            free_model(info.primary_type, &m);
            return st;
        }
        if (out_single_model) *out_single_model = m;
        if (out_single_type) *out_single_type = info.primary_type;
        return 1;
    } else {
        MODEL_BUNDLE *b = ImportModelBundle(filepath);
        if (!b) return LSCM_ERR_IO;
        if (out_bundle) *out_bundle = b;
        return (int)b->count;
    }
}

/* =========================================================================
 * Convenience Wrappers for Individual Models
 * ========================================================================= */

int SavePCAModel(const char *filepath, const PCAMODEL *pca) {
    return SaveModel(filepath, "pca", MODEL_TYPE_PCA, pca);
}
int ImportPCAModel(const char *filepath, PCAMODEL *pca) {
    return ImportModel(filepath, MODEL_TYPE_PCA, pca);
}

int SavePLSModel(const char *filepath, const PLSMODEL *pls) {
    return SaveModel(filepath, "pls", MODEL_TYPE_PLS, pls);
}
int ImportPLSModel(const char *filepath, PLSMODEL *pls) {
    return ImportModel(filepath, MODEL_TYPE_PLS, pls);
}

int SaveCPCAModel(const char *filepath, const CPCAMODEL *cpca) {
    return SaveModel(filepath, "cpca", MODEL_TYPE_CPCA, cpca);
}
int ImportCPCAModel(const char *filepath, CPCAMODEL *cpca) {
    return ImportModel(filepath, MODEL_TYPE_CPCA, cpca);
}

int SaveMLRModel(const char *filepath, const MLRMODEL *mlr) {
    return SaveModel(filepath, "mlr", MODEL_TYPE_MLR, mlr);
}
int ImportMLRModel(const char *filepath, MLRMODEL *mlr) {
    return ImportModel(filepath, MODEL_TYPE_MLR, mlr);
}

int SaveLDAModel(const char *filepath, const LDAMODEL *lda) {
    return SaveModel(filepath, "lda", MODEL_TYPE_LDA, lda);
}
int ImportLDAModel(const char *filepath, LDAMODEL *lda) {
    return ImportModel(filepath, MODEL_TYPE_LDA, lda);
}

int SaveICAModel(const char *filepath, const ICAMODEL *ica) {
    return SaveModel(filepath, "ica", MODEL_TYPE_ICA, ica);
}
int ImportICAModel(const char *filepath, ICAMODEL *ica) {
    return ImportModel(filepath, MODEL_TYPE_ICA, ica);
}

int SaveUPCAModel(const char *filepath, const UPCAMODEL *upca) {
    return SaveModel(filepath, "upca", MODEL_TYPE_UPCA, upca);
}
int ImportUPCAModel(const char *filepath, UPCAMODEL *upca) {
    return ImportModel(filepath, MODEL_TYPE_UPCA, upca);
}

int SaveUPLSModel(const char *filepath, const UPLSMODEL *upls) {
    return SaveModel(filepath, "upls", MODEL_TYPE_UPLS, upls);
}
int ImportUPLSModel(const char *filepath, UPLSMODEL *upls) {
    return ImportModel(filepath, MODEL_TYPE_UPLS, upls);
}

const char* ModelTypeToString(ModelType type) {
    switch (type) {
        case MODEL_TYPE_PCA:  return "PCA";
        case MODEL_TYPE_PLS:  return "PLS";
        case MODEL_TYPE_CPCA: return "CPCA";
        case MODEL_TYPE_MLR:  return "MLR";
        case MODEL_TYPE_LDA:  return "LDA";
        case MODEL_TYPE_ICA:  return "ICA";
        case MODEL_TYPE_UPCA: return "UPCA";
        case MODEL_TYPE_UPLS: return "UPLS";
        default: return "UNKNOWN";
    }
}

ModelType ModelTypeFromString(const char *str) {
    if (!str) return MODEL_TYPE_UNKNOWN;
    if (strcmp(str, "PCA") == 0 || strcmp(str, "pca") == 0) return MODEL_TYPE_PCA;
    if (strcmp(str, "PLS") == 0 || strcmp(str, "pls") == 0) return MODEL_TYPE_PLS;
    if (strcmp(str, "CPCA") == 0 || strcmp(str, "cpca") == 0) return MODEL_TYPE_CPCA;
    if (strcmp(str, "MLR") == 0 || strcmp(str, "mlr") == 0) return MODEL_TYPE_MLR;
    if (strcmp(str, "LDA") == 0 || strcmp(str, "lda") == 0) return MODEL_TYPE_LDA;
    if (strcmp(str, "ICA") == 0 || strcmp(str, "ica") == 0) return MODEL_TYPE_ICA;
    if (strcmp(str, "UPCA") == 0 || strcmp(str, "upca") == 0) return MODEL_TYPE_UPCA;
    if (strcmp(str, "UPLS") == 0 || strcmp(str, "upls") == 0) return MODEL_TYPE_UPLS;
    return MODEL_TYPE_UNKNOWN;
}
