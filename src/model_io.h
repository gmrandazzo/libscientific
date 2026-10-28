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

#ifndef MODEL_IO_H
#define MODEL_IO_H

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#include "pca.h"
#include "pls.h"
#include "cpca.h"
#include "mlr.h"
#include "lda.h"
#include "ica.h"
#include "upca.h"
#include "upls.h"

#ifdef __cplusplus
extern "C" {
#endif

#define LSCM_MAGIC 0x4D43534CU /* "LSCM" in little-endian */
#define LSCM_VERSION 1U
#define LSCM_MAX_NAME_LEN 64

/**
 * Enumeration of supported model types in the LSCM container format.
 */
typedef enum {
    MODEL_TYPE_UNKNOWN = 0,
    MODEL_TYPE_PCA     = 1,
    MODEL_TYPE_PLS     = 2,
    MODEL_TYPE_CPCA    = 3,
    MODEL_TYPE_MLR     = 4,
    MODEL_TYPE_LDA     = 5,
    MODEL_TYPE_ICA     = 6,
    MODEL_TYPE_UPCA    = 7,
    MODEL_TYPE_UPLS    = 8,
    MODEL_TYPE_ANY     = 255
} ModelType;

/**
 * Error codes returned by model I/O routines.
 */
typedef enum {
    LSCM_OK                 =  0,
    LSCM_ERR_FILE_NOT_FOUND = -1,
    LSCM_ERR_IO             = -2,
    LSCM_ERR_INVALID_MAGIC  = -3,
    LSCM_ERR_UNSUPPORTED_VER= -4,
    LSCM_ERR_INVALID_TYPE   = -5,
    LSCM_ERR_TYPE_MISMATCH  = -6,
    LSCM_ERR_NOT_FOUND      = -7,
    LSCM_ERR_MEM            = -8,
    LSCM_ERR_PARAM          = -9,
    LSCM_ERR_CORRUPT        = -10
} LSCMStatus;

/**
 * Metadata descriptor for inspecting files without full payload loading.
 */
typedef struct {
    uint32_t num_models;
    bool is_single;
    ModelType primary_type;
    char primary_name[LSCM_MAX_NAME_LEN];
    uint32_t format_version;
    uint64_t total_size;
} MODEL_FILE_INFO;

/**
 * Single entry within a multi-model bundle.
 */
typedef struct {
    char name[LSCM_MAX_NAME_LEN];
    ModelType type;
    void *model;
} MODEL_ENTRY;

/**
 * Multi-model bundle container representing a project or pipeline.
 */
typedef struct {
    size_t count;
    size_t capacity;
    MODEL_ENTRY *entries;
} MODEL_BUNDLE;

/* =========================================================================
 * Core Single-Model APIs (.lsm)
 * ========================================================================= */

/**
 * Save a single model into a standalone .lsm file.
 *
 * @param [in] filepath Destination file path (typically with .lsm extension)
 * @param [in] model_name Identifier tag for the model (e.g. "qsar_pls" or NULL for "default")
 * @param [in] type Model type (MODEL_TYPE_PCA, MODEL_TYPE_PLS, etc.)
 * @param [in] model Pointer to the populated model struct
 * @return LSCM_OK (0) on success, or a negative LSCMStatus error code
 */
int SaveModel(const char *filepath, const char *model_name, ModelType type, const void *model);

/**
 * Import a single model from a file (.lsm or .lsb with 1 model).
 *
 * @param [in] filepath Source file path
 * @param [in] expected_type Expected ModelType (or MODEL_TYPE_ANY)
 * @param [out] model Pointer to pre-allocated or initialized target model struct
 * @return LSCM_OK (0) on success, or a negative LSCMStatus error code
 */
int ImportModel(const char *filepath, ModelType expected_type, void *model);

/* =========================================================================
 * Multi-Model Project Bundle APIs (.lsb)
 * ========================================================================= */

/**
 * Create a new empty model bundle in memory.
 */
MODEL_BUNDLE* NewModelBundle(void);

/**
 * Add a model into a bundle. The bundle takes custody of the model pointer.
 *
 * @param [in,out] bundle Target bundle
 * @param [in] name Unique model identifier name within the project
 * @param [in] type Model type
 * @param [in] model Pointer to the concrete model struct
 * @return LSCM_OK (0) on success, or a negative LSCMStatus error code
 */
int ModelBundleAddModel(MODEL_BUNDLE *bundle, const char *name, ModelType type, void *model);

/**
 * Retrieve a model from a bundle by its unique name.
 *
 * @param [in] bundle Source bundle
 * @param [in] name Unique model name
 * @param [out] out_type Optional output parameter for model type
 * @return Pointer to the model struct, or NULL if not found
 */
void* ModelBundleGetModel(const MODEL_BUNDLE *bundle, const char *name, ModelType *out_type);

/**
 * Retrieve a model from a bundle by index.
 *
 * @param [in] bundle Source bundle
 * @param [in] index 0-based index
 * @param [out] out_name Optional output buffer/pointer for model name
 * @param [out] out_type Optional output parameter for model type
 * @return Pointer to the model struct, or NULL if index is out of bounds
 */
void* ModelBundleGetModelByIndex(const MODEL_BUNDLE *bundle, size_t index, const char **out_name, ModelType *out_type);

/**
 * Free a model bundle container.
 *
 * @param [in,out] bundle Pointer to bundle pointer
 * @param [in] free_contained_models If non-zero, frees each contained model struct appropriately
 */
void DelModelBundle(MODEL_BUNDLE **bundle, int free_contained_models);

/**
 * Save multiple models from a bundle into a single .lsb file.
 *
 * @param [in] filepath Destination file path (typically with .lsb extension)
 * @param [in] bundle Model bundle containing 1 or more models
 * @return LSCM_OK (0) on success, or a negative LSCMStatus error code
 */
int SaveModelBundle(const char *filepath, const MODEL_BUNDLE *bundle);

/**
 * Import a multi-model project bundle from an .lsb file.
 *
 * @param [in] filepath Source file path
 * @return Pointer to loaded MODEL_BUNDLE, or NULL on error
 */
MODEL_BUNDLE* ImportModelBundle(const char *filepath);

/* =========================================================================
 * Inspection & Auto-Detection APIs
 * ========================================================================= */

/**
 * Inspect a model file header without loading full matrix/tensor payloads.
 *
 * @param [in] filepath Source file path (.lsm or .lsb)
 * @param [out] info Output structure with model count, types, and sizes
 * @return LSCM_OK (0) on success, or a negative LSCMStatus error code
 */
int InspectModelFile(const char *filepath, MODEL_FILE_INFO *info);

/**
 * Automatically import a file, detecting whether it is a single model or a bundle.
 *
 * If single model:
 *   - `*out_bundle` will be NULL.
 *   - `*out_single_model` receives the newly allocated model struct.
 *   - `*out_single_type` receives the model type.
 *   - Returns 1.
 *
 * If multi-model bundle:
 *   - `*out_single_model` will be NULL.
 *   - `*out_bundle` receives the loaded MODEL_BUNDLE.
 *   - Returns the total count of models in the bundle (> 1).
 *
 * @param [in] filepath Source file path (.lsm, .lsb, or generic)
 * @param [out] out_single_model Output pointer for single model
 * @param [out] out_single_type Output pointer for single model type
 * @param [out] out_bundle Output pointer for bundle
 * @return Number of models imported on success (> 0), or a negative error code
 */
int ImportModelAuto(const char *filepath,
                    void **out_single_model,
                    ModelType *out_single_type,
                    MODEL_BUNDLE **out_bundle);

/* =========================================================================
 * Type-Specific Convenience Helpers
 * ========================================================================= */

int SavePCAModel(const char *filepath, const PCAMODEL *pca);
int ImportPCAModel(const char *filepath, PCAMODEL *pca);

int SavePLSModel(const char *filepath, const PLSMODEL *pls);
int ImportPLSModel(const char *filepath, PLSMODEL *pls);

int SaveCPCAModel(const char *filepath, const CPCAMODEL *cpca);
int ImportCPCAModel(const char *filepath, CPCAMODEL *cpca);

int SaveMLRModel(const char *filepath, const MLRMODEL *mlr);
int ImportMLRModel(const char *filepath, MLRMODEL *mlr);

int SaveLDAModel(const char *filepath, const LDAMODEL *lda);
int ImportLDAModel(const char *filepath, LDAMODEL *lda);

int SaveICAModel(const char *filepath, const ICAMODEL *ica);
int ImportICAModel(const char *filepath, ICAMODEL *ica);

int SaveUPCAModel(const char *filepath, const UPCAMODEL *upca);
int ImportUPCAModel(const char *filepath, UPCAMODEL *upca);

int SaveUPLSModel(const char *filepath, const UPLSMODEL *upls);
int ImportUPLSModel(const char *filepath, UPLSMODEL *upls);

/**
 * Helper to convert ModelType enum to human-readable string.
 */
const char* ModelTypeToString(ModelType type);

/**
 * Helper to parse human-readable string to ModelType enum.
 */
ModelType ModelTypeFromString(const char *str);

#ifdef __cplusplus
}
#endif

#endif /* MODEL_IO_H */
