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

.. autocstruct:: pca.h::PCAMODEL

.. autocfunction:: pca.h::NewPCAModel

.. autocfunction:: pca.h::DelPCAModel

.. autocfunction:: pca.h::PCA

.. autocfunction:: pca.h::PCAScorePredictor

.. autocfunction:: pca.h::PCAIndVarPredictor

.. autocfunction:: pca.h::GetResidualMatrix

.. autocfunction:: pca.h::PrintPCA

CPCA Functions
--------------

.. autocstruct:: cpca.h::CPCAMODEL

.. autocfunction:: cpca.h::NewCPCAModel

.. autocfunction:: cpca.h::DelCPCAModel

.. autocfunction:: cpca.h::CPCA

.. autocfunction:: cpca.h::CPCAScorePredictor

.. autocfunction:: cpca.h::PrintCPCA
    
PLS Functions
-------------

.. autocstruct:: pls.h::PLSMODEL

.. autocfunction:: pls.h::NewPLSModel

.. autocfunction:: pls.h::DelPLSModel

.. autocfunction:: pls.h::PLS

.. autocfunction:: pls.h::PLSBetasCoeff

.. autocfunction:: pls.h::PLSScorePredictor

.. autocfunction:: pls.h::PLSYPredictor

.. autocfunction:: pls.h::PLSYPredictorAllLV

.. autocfunction:: pls.h::PLSRegressionStatistics

.. autocfunction:: pls.h::PLSDiscriminantAnalysisStatistics

.. autocfunction:: pls.h::GetLVCCutoff

.. autocfunction:: pls.h::PrintPLSModel


MLR Functions
-------------

.. autocstruct:: mlr.h::MLRMODEL

.. autocfunction:: mlr.h::NewMLRModel

.. autocfunction:: mlr.h::DelMLRModel

.. autocfunction:: mlr.h::MLR

.. autocfunction:: mlr.h::MLRPredictY

.. autocfunction:: mlr.h::MLRRegressionStatistics

.. autocfunction:: mlr.h::PrintMLR

LDA Functions
-------------

.. autocstruct:: lda.h::LDAMODEL

.. autocfunction:: lda.h::NewLDAModel

.. autocfunction:: lda.h::DelLDAModel

.. autocfunction:: lda.h::LDA

.. autocfunction:: lda.h::LDAPrediction

.. autocfunction:: lda.h::LDAStatistics

.. autocfunction:: lda.h::LDAMulticlassStatistics   

.. autocfunction:: lda.h::PrintLDAModel


UPCA Functions
--------------

.. autocstruct:: upca.h::UPCAMODEL

.. autocfunction:: upca.h::NewUPCAModel

.. autocfunction:: upca.h::DelUPCAModel

.. autocfunction:: upca.h::UPCA

.. autocfunction:: upca.h::UPCAScorePredictor

.. autocfunction:: upca.h::UPCAIndVarPredictor

.. autocfunction:: upca.h::PrintUPCAModel   


Object Selection Functions
--------------------------
