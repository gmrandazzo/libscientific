/* Implements Linear Discriminant Analysis (LDA).
 * Copyright (C) 2016-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
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

#ifndef LDA_H
#define LDA_H
#include <stdio.h>
#include "matrix.h"
#include "vector.h"
#include "tensor.h"
#include "scientificinfo.h"


/**
 * Linear Discriminant Analysis model (Fisher LDA).
 *
 * Create it with :c:func:`NewLDAModel`, fill it with :c:func:`LDA` and free
 * it with :c:func:`DelLDAModel`. With ``n`` objects, ``f`` features and
 * ``K`` classes (labels ``0..K-1`` or ``1..K``):
 *
 * Members:
 *   inv_cov: Inverse (or pseudo-inverse) of the pooled within-class scatter
 *     matrix, ``f x f``.
 *   features: Training objects projected on the discriminant directions:
 *     one ``n_k x f`` matrix per class.
 *   mnpdf: Normal density of each projected training object in its class,
 *     same layout as ``features``.
 *   evect: Discriminant directions, ``f x f`` (one per column).
 *   mu: Class means, ``K x f``.
 *   fmean: Mean of the projected features per class, ``K x f``.
 *   fsdev: Standard deviation of the projected features per class, ``K x f``.
 *   yscrambling: Y-scrambling results (validation functions).
 *   recalculated_y: Fitted classes (validation functions).
 *   recalculated_residuals: Fitted-class residuals (validation functions).
 *   predicted_y: Cross-validated classes (validation functions).
 *   predicted_residuals: Cross-validated residuals (validation functions).
 *   eval: Eigenvalue of each discriminant direction (size ``f``).
 *   pprob: Prior probability of each class (fraction of training objects).
 *   classid: Row indexes of the training objects, grouped by class.
 *   nclass: Number of classes ``K``.
 *   class_start: First class label, 0 or 1.
 *   roc: ROC curves (validation functions).
 *   pr: Precision-recall curves (validation functions).
 *   roc_aucs: ROC areas under the curve (validation functions).
 *   pr_aucs: Precision-recall areas under the curve (validation functions).
 */
typedef struct{
  matrix *inv_cov;
  tensor *features;
  tensor *mnpdf;
  matrix *evect;
  matrix *mu;
  matrix *fmean;
  matrix *fsdev;
  matrix *yscrambling;
  matrix *recalculated_y;
  matrix *recalculated_residuals;
  matrix *predicted_y;
  matrix *predicted_residuals;
  dvector *eval;
  dvector *pprob;
  uivector *classid;
  size_t nclass;
  size_t class_start;
  tensor *roc;
  tensor *pr;
  dvector *roc_aucs;
  dvector *pr_aucs;
} LDAMODEL;

/**
 * Allocate an empty LDA model.
 *
 * Args:
 *   m: [out] Address of the pointer that receives the new model.
 *
 * Note:
 *   Release the model with :c:func:`DelLDAModel`.
 */
void NewLDAModel(LDAMODEL **m);

/**
 * Free an LDA model and all its members, and set the pointer to ``NULL``.
 *
 * Args:
 *   m: [in,out] Address of the model pointer created with :c:func:`NewLDAModel`.
 */
void DelLDAModel(LDAMODEL **m);

/**
 * Print the content of an LDA model to standard output.
 *
 * Args:
 *   m: [in] The model to print.
 */
void PrintLDAModel(LDAMODEL *m);

/**
 * Fit a Fisher linear discriminant analysis model.
 *
 * Computes the class means and priors, the within-class (``Sw``) and
 * between-class (``Sb``) scatter matrices and the discriminant directions
 * (eigenvectors of ``Sw^-1 Sb``). When ``Sw`` is singular its
 * pseudo-inverse is used (Liu et al., Int. J. Patt. Recogn. Artif. Intell.
 * 21, 1265, 2007). Fitting an already fitted model replaces its results.
 *
 * Args:
 *   mx: [in] Training objects, ``n x f``.
 *   my: [in] Class labels, ``n x 1``: integers ``0..K-1`` or ``1..K``, with
 *     every class present.
 *   lda: [out] Model created with :c:func:`NewLDAModel`; filled with the
 *     results (see :c:type:`LDAMODEL`).
 *
 * Returns:
 *   Nothing, the results are written into ``lda``.
 */
void LDA(matrix *mx, matrix *my, LDAMODEL *lda);

/**
 * Classify new objects with an LDA model.
 *
 * The discriminant score of class ``k`` is
 * ``mu_k' S^-1 x - 0.5 mu_k' S^-1 mu_k + ln(prior_k)``; the predicted class
 * is the one with the highest score.
 *
 * Args:
 *   mx: [in] Objects, ``m x f``.
 *   lda: [in] Model fitted with :c:func:`LDA`.
 *   pfeatures: [out] **Empty** matrix created with ``initMatrix``: receives
 *     the projections on the discriminant directions, ``m x f``.
 *   probability: [out] Matrix created with ``initMatrix``; resized to
 *     ``m x K``: discriminant score of each class (log scale, not a
 *     probability).
 *   mnpdf: [out] Matrix created with ``initMatrix``; resized to ``m x f``:
 *     normal density of each projection in the predicted class.
 *   prediction: [out] Matrix created with ``initMatrix``; resized to
 *     ``m x 1``: predicted labels, in the numbering of the training labels.
 *
 * Returns:
 *   Nothing, the results are written into the outputs.
 */
void LDAPrediction(matrix *mx,
                   LDAMODEL *lda,
                   matrix *pfeatures,
                   matrix *probability,
                   matrix *mnpdf,
                   matrix *prediction);

/**
 * Compute ROC and precision-recall curves of a binary classification.
 *
 * Args:
 *   y_true: [in] True classes (0/1), size ``n``.
 *   y_pred: [in] Predicted class scores, size ``n``.
 *   roc: [out] Matrix created with ``initMatrix``; receives the ROC curve
 *     (columns: false and true positive rate).
 *   roc_auc: [out] Area under the ROC curve.
 *   precision_recal: [out] Matrix created with ``initMatrix``; receives the
 *     precision-recall curve (columns: recall and precision).
 *   pr_auc: [out] Area under the precision-recall curve.
 *
 * Returns:
 *   Nothing, the results are written into the outputs.
 */
void LDAStatistics(dvector *y_true,
                   dvector *y_pred,
                   matrix *roc,
                   double *roc_auc,
                   matrix *precision_recal,
                   double *pr_auc);

/**
 * Compute one-vs-rest ROC and precision-recall curves from class labels.
 *
 * For every class, the objects of that class are the positives and the
 * predicted label is used as a 0/1 score. With two classes a single curve
 * (first class vs the other) is computed.
 *
 * Args:
 *   y_true: [in] True labels, ``n x 1`` (``0..K-1`` or ``1..K``).
 *   y_pred: [in] Predicted labels, ``n x 1``, e.g. the ``prediction``
 *     output of :c:func:`LDAPrediction`.
 *   roc: [out] Tensor created with ``initTensor``; one ROC curve per class
 *     is appended. ``NULL`` to skip.
 *   roc_aucs: [out] Vector created with ``initDVector``; one ROC AUC per
 *     class is appended.
 *   precision_recals: [out] Tensor; one precision-recall curve per class is
 *     appended. ``NULL`` to skip.
 *   pr_aucs: [out] Vector; one precision-recall AUC per class is appended.
 *
 * Returns:
 *   Nothing, the results are written into the outputs.
 */
void LDAMulticlassStatistics(matrix *y_true,
                             matrix *y_pred,
                             tensor *roc,
                             dvector *roc_aucs,
                             tensor *precision_recals,
                             dvector *pr_aucs);

/**
 * Count the distinct class labels in a label matrix.
 *
 * Args:
 *   my: [in] Class labels, ``n x 1``.
 *
 * Returns:
 *   The number of distinct labels.
 */
int getNClasses(matrix *my);

#endif
