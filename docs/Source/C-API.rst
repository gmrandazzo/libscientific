***************
C API Reference
***************


Matrix Functions
----------------

.. autocstruct:: matrix.h::matrix

.. autocfunction:: matrix.h::NewMatrix

.. autocfunction:: matrix.h::DelMatrix

.. autocfunction:: matrix.h::initMatrix

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

.. autocfunction:: ldupcaa.h::PrintUPCAModel   
