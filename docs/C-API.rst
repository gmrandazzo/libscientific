***************
C API Reference
***************

Double Vector Functions
-----------------------

.. autocstruct:: vector.h::dvector

.. autocfunction:: vector.h::NewDVector

.. autocfunction:: vector.h::DelDVector

.. autocfunction:: vector.h::initDVector

.. autocfunction:: vector.h::DVectorAppend

.. autocfunction:: vector.h::DVectorRemoveAt

.. autocfunction:: vector.h::PrintDVector


Unsigned int Vector Functions
-----------------------------

.. autocstruct:: vector.h::uivector

.. autocfunction:: vector.h::NewUIVector

.. autocfunction:: vector.h::DelUIVector

.. autocfunction:: vector.h::initUIVector

.. autocfunction:: vector.h::UIVectorAppend

.. autocfunction:: vector.h::UIVectorRemoveAt

.. autocfunction:: vector.h::PrintUIVector


Int Vector Functions
--------------------
.. autocstruct:: vector.h::ivector

.. autocfunction:: vector.h::NewIVector

.. autocfunction:: vector.h::DelIVector

.. autocfunction:: vector.h::initIVector

.. autocfunction:: vector.h::IVectorAppend

.. autocfunction:: vector.h::IVectorRemoveAt

.. autocfunction:: vector.h::PrintIVector


String Vector Functions
-----------------------

.. autocstruct:: vector.h::strvector

.. autocfunction:: vector.h::NewStrVector

.. autocfunction:: vector.h::DelStrVector

.. autocfunction:: vector.h::initStrVector

.. autocfunction:: vector.h::StrVectorAppend

.. autocfunction:: vector.h::PrintStrVector

Matrix Functions
----------------

.. autocstruct:: matrix.h::matrix

.. autocfunction:: matrix.h::NewMatrix

.. autocfunction:: matrix.h::DelMatrix

.. autocfunction:: matrix.h::initMatrix

.. autocfunction:: matrix.h::MatrixAppendRow

.. autocfunction:: matrix.h::MatrixAppendCol

.. autocfunction:: matrix.h::MatrixDeleteRowAt

.. autocfunction:: matrix.h::MatrixDeleteColAt

.. autocfunction:: matrix.h::PrintMatrix

Tensor Functions
----------------

.. autocstruct:: tensor.h::tensor

.. autocfunction:: tensor.h::NewTensor

.. autocfunction:: tensor.h::NewTensorMatrix

.. autocfunction:: tensor.h::AddTensorMatrix

.. autocfunction:: tensor.h::DelTensor

.. autocfunction:: tensor.h::initTensor

.. autocfunction:: tensor.h::PrintTensor

PCA Functions
-------------

Principal Component Analysis with the NIPALS algorithm (``#include <scientific.h>``).
Typical workflow: :c:func:`NewPCAModel` → :c:func:`PCA` → read the members of
:c:type:`PCAMODEL` or :c:func:`PCAScorePredictor` for new objects →
:c:func:`DelPCAModel`.

.. autocstruct:: pca.h::PCAMODEL
   :members:

.. autocfunction:: pca.h::NewPCAModel

.. autocfunction:: pca.h::DelPCAModel

.. autocfunction:: pca.h::PCA

.. autocfunction:: pca.h::PCAScorePredictor

.. autocfunction:: pca.h::PCAIndVarPredictor

.. autocfunction:: pca.h::PCARSquared

.. autocfunction:: pca.h::PCARankValidation

.. autocfunction:: pca.h::PCATsqContributions

.. autocfunction:: pca.h::GetResidualMatrix

.. autocfunction:: pca.h::PrintPCA

Complete example (``docs/c_code_examples/pcaexample1.c``):

.. literalinclude:: c_code_examples/pcaexample1.c
   :language: c

CPCA Functions
--------------

Consensus (multi-block) Principal Component Analysis. The blocks are the
layers of a tensor, all with the same objects (rows). Workflow:
:c:func:`NewCPCAModel` → :c:func:`CPCA` → :c:func:`CPCAScorePredictor` →
:c:func:`DelCPCAModel`.

.. autocstruct:: cpca.h::CPCAMODEL
   :members:

.. autocfunction:: cpca.h::NewCPCAModel
.. autocfunction:: cpca.h::DelCPCAModel
.. autocfunction:: cpca.h::CPCA
.. autocfunction:: cpca.h::CPCAScorePredictor
.. autocfunction:: cpca.h::PrintCPCA

Complete example (``docs/c_code_examples/cpcaexample1.c``):

.. literalinclude:: c_code_examples/cpcaexample1.c
   :language: c

PLS Functions
-------------

Partial Least Squares regression and discriminant analysis (PLS-DA) with the
NIPALS algorithm. Workflow: :c:func:`NewPLSModel` → :c:func:`PLS` →
:c:func:`PLSYPredictorAllLV` (or :c:func:`PLSScorePredictor` +
:c:func:`PLSYPredictor`) → :c:func:`DelPLSModel`. Predictions for all the
latent variables are laid out as ``n x (ny*nlv)``: column ``ny*lv + j`` is
response ``j`` with ``lv+1`` latent variables.

.. autocstruct:: pls.h::PLSMODEL
   :members:

.. autocfunction:: pls.h::NewPLSModel
.. autocfunction:: pls.h::DelPLSModel
.. autocfunction:: pls.h::PLS
.. autocfunction:: pls.h::PLSBetasCoeff
.. autocfunction:: pls.h::PLSScorePredictor
.. autocfunction:: pls.h::PLSYPredictor
.. autocfunction:: pls.h::PLSYPredictorAllLV
.. autocfunction:: pls.h::PLSRegressionStatistics
.. autocfunction:: pls.h::PLSDiscriminantAnalysisStatistics
.. autocfunction:: pls.h::PLSVIP
.. autocfunction:: pls.h::GetLVCCutoff
.. autocfunction:: pls.h::PrintPLSModel

Complete example (``docs/c_code_examples/plsexample1.c``):

.. literalinclude:: c_code_examples/plsexample1.c
   :language: c

MLR Functions
-------------

Multiple Linear Regression by ordinary least squares. Workflow:
:c:func:`NewMLRModel` → :c:func:`MLR` → :c:func:`MLRPredictY` →
:c:func:`DelMLRModel`.

.. autocstruct:: mlr.h::MLRMODEL
   :members:

.. autocfunction:: mlr.h::NewMLRModel
.. autocfunction:: mlr.h::DelMLRModel
.. autocfunction:: mlr.h::MLR
.. autocfunction:: mlr.h::MLRPredictY
.. autocfunction:: mlr.h::MLRRegressionStatistics
.. autocfunction:: mlr.h::PrintMLR

Complete example (``docs/c_code_examples/mlrexample1.c``):

.. literalinclude:: c_code_examples/mlrexample1.c
   :language: c

LDA Functions
-------------

Fisher Linear Discriminant Analysis for classification. Class labels are
integers ``0..K-1`` or ``1..K``. Workflow: :c:func:`NewLDAModel` →
:c:func:`LDA` → :c:func:`LDAPrediction` → :c:func:`DelLDAModel`.

.. autocstruct:: lda.h::LDAMODEL
   :members:

.. autocfunction:: lda.h::NewLDAModel
.. autocfunction:: lda.h::DelLDAModel
.. autocfunction:: lda.h::LDA
.. autocfunction:: lda.h::LDAPrediction
.. autocfunction:: lda.h::LDAStatistics
.. autocfunction:: lda.h::LDAMulticlassStatistics
.. autocfunction:: lda.h::getNClasses
.. autocfunction:: lda.h::PrintLDAModel

Complete example (``docs/c_code_examples/ldaexample1.c``):

.. literalinclude:: c_code_examples/ldaexample1.c
   :language: c

UPCA Functions
--------------

Multi-way (unfolded) Principal Component Analysis of three-way data: a
tensor of layers, all ``n x f``. Workflow: :c:func:`NewUPCAModel` →
:c:func:`UPCA` → :c:func:`UPCAScorePredictor` / :c:func:`UPCAIndVarPredictor`
→ :c:func:`DelUPCAModel`.

.. autocstruct:: upca.h::UPCAMODEL
   :members:

.. autocfunction:: upca.h::NewUPCAModel
.. autocfunction:: upca.h::DelUPCAModel
.. autocfunction:: upca.h::CheckTensor
.. autocfunction:: upca.h::UPCA
.. autocfunction:: upca.h::UPCAScorePredictor
.. autocfunction:: upca.h::UPCAIndVarPredictor
.. autocfunction:: upca.h::PrintUPCAModel

Complete example (``docs/c_code_examples/upcaexample1.c``):

.. literalinclude:: c_code_examples/upcaexample1.c
   :language: c

UPLS Functions
--------------

Multi-way (unfolded) Partial Least Squares regression. X and Y are tensors
of layers with the same objects. Workflow: :c:func:`NewUPLSModel` →
:c:func:`UPLS` → :c:func:`UPLSScorePredictor` + :c:func:`UPLSYPredictor` →
:c:func:`DelUPLSModel`; :c:func:`UPLSRandomGroupsCV` and
:c:func:`UPLSLOOCV` estimate the predictive ability.

.. autocstruct:: upls.h::UPLSMODEL
   :members:

.. autocfunction:: upls.h::NewUPLSModel
.. autocfunction:: upls.h::DelUPLSModel
.. autocfunction:: upls.h::UPLS
.. autocfunction:: upls.h::UPLSScorePredictor
.. autocfunction:: upls.h::UPLSYPredictor
.. autocfunction:: upls.h::UPLSRSquared
.. autocfunction:: upls.h::UPLSRSquared_SSErr_SStot
.. autocfunction:: upls.h::UPLSRandomGroupsCV
.. autocfunction:: upls.h::UPLSLOOCV
.. autocfunction:: upls.h::UPLSYScrambling
.. autocfunction:: upls.h::UPLSGetPCModelCutOff
.. autocfunction:: upls.h::PrintUPLSModel

Complete example (``docs/c_code_examples/uplsexample1.c``):

.. literalinclude:: c_code_examples/uplsexample1.c
   :language: c

Model Container I/O Functions
-----------------------------

.. autocstruct:: model_io.h::MODEL_FILE_INFO

.. autocstruct:: model_io.h::MODEL_ENTRY

.. autocstruct:: model_io.h::MODEL_BUNDLE

.. autocfunction:: model_io.h::SaveModel

.. autocfunction:: model_io.h::ImportModel

.. autocfunction:: model_io.h::NewModelBundle

.. autocfunction:: model_io.h::ModelBundleAddModel

.. autocfunction:: model_io.h::ModelBundleGetModel

.. autocfunction:: model_io.h::ModelBundleGetModelByIndex

.. autocfunction:: model_io.h::DelModelBundle

.. autocfunction:: model_io.h::SaveModelBundle

.. autocfunction:: model_io.h::ImportModelBundle

.. autocfunction:: model_io.h::InspectModelFile

.. autocfunction:: model_io.h::ImportModelAuto

.. autocfunction:: model_io.h::SavePCAModel

.. autocfunction:: model_io.h::ImportPCAModel

.. autocfunction:: model_io.h::SavePLSModel

.. autocfunction:: model_io.h::ImportPLSModel

.. autocfunction:: model_io.h::SaveCPCAModel

.. autocfunction:: model_io.h::ImportCPCAModel

.. autocfunction:: model_io.h::SaveMLRModel

.. autocfunction:: model_io.h::ImportMLRModel

.. autocfunction:: model_io.h::SaveLDAModel

.. autocfunction:: model_io.h::ImportLDAModel

.. autocfunction:: model_io.h::SaveICAModel

.. autocfunction:: model_io.h::ImportICAModel

.. autocfunction:: model_io.h::SaveUPCAModel

.. autocfunction:: model_io.h::ImportUPCAModel

.. autocfunction:: model_io.h::SaveUPLSModel

.. autocfunction:: model_io.h::ImportUPLSModel

.. autocfunction:: model_io.h::ModelTypeToString

.. autocfunction:: model_io.h::ModelTypeFromString
