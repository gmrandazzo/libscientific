.. libscientific documentation master file, created by
   sphinx-quickstart on Tue Jul 12 10:35:29 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


What is libscientific?
======================

Libscientific aims to do numerical/multivariate and statistical analysis. The library is written in C language for almost all the computations except for a few methods such as the singular value decomposition (SVD) and the eigenvector/eigenvalue transformation from a third-party library, the Lapack library. It requires two dependencies: CMake and a c/fortran compiler and supports Windows, Linux, macOS, and embedded systems. Python bindings are available using ctypes to avoid dependencies. The library is distributed under GPLv3 license allowing it to be used for public and commercial purposes.

The significant advantage of libscientific is that it does not require dependencies unless a c/fortran compiler and CMake. Moreover, the size of the library is under 1.5MB.
Libscientific also has been tested with old computers.


What can I do with libscientific?
---------------------------------

* Multivariate analysis:
        * Principal Component Analysis using the NIPALS algorithm
        * Partial Least Squares (PLS) using the NIPALS algorithm (Regression and Classification)
        * Multiple Linear Regression (MLR) using the Ordinary Least Squares algorithm
        * Linear Discriminant Analysis (LDA) using the Fisher algorithm
        * Multi-waw unfolding PCA (UPCA) using the S. Wold, P. Geladi and K. Esbensen algorithm
        * Multi-way unfolding PLS (UPLS) using the S. Wold, P. Geladi and K. Esbensen algorithm

* Matrix/Vector/Tensor computations:
        * Matrix/Vector and Vector/Matrix dot product
        * Matrix/Matrix dot product
        * Vector/Vector dot product
        * Matrix inversion using the Gauss-Jordan algorithm or the LU decomposition
        * Matrix pseudo inversion using the Moore-Penrose formula
        * Matrix transpose
        * Tensor/Tensor dot product
        * Vector/Matrix Kronecker product
        * Tensor/Matrix dot product
        * Tensor transpose

* Statistical analysis:
        * Descriptive statistics of a matrix
        * R squared (R^2)
        * Mean absolute error (MAE)
        * Mean squared error (MSE)
        * Root mean squared error (RMSE)
        * Bias estimation (BIAS)
        * Sensitivity binary classification test
        * Positive predicted values binary classification test
        * Receiver operating characteristic curve (ROC)
        * Precision-Recal curve (PR)

* Design of experiment:
        * Bifactorial matrix expansion used in design of experiments (DOE)
        * Yates analysis to assess the variable effects in a DOE

* Clustering:
        * K-Means with different initialization (Random, ++ and so on) for divisive clustering
        * Hierarchical clustering for aglomerative clustering
        * Hyper grid map which create clusters dividing the hyperspace into a hyper grid and abundance score plot

* Instance/Object selection:
        * Most descriptive compound selection
        * Maximum dissimilarity selection

* Optimization:
        * Downhill optimization using the Nelder–Mead method aka Simplex method

* Interpolation:
        * Natural cubic spline interpolation

* Similarity analysis:
        * Euclidean distance analysis
        * Manhattan distance analysis
        * Cosine distance analysis

* Model validation:
        * Leave-One-Out
        * K-Fold cross validation
        * Repeated K-Fold cross validation
        * Train/Test split

* Variable selection:
        * Genetic agorithm variable selection for PLS only
        * Particle Swarm Optimization variable selection for PLS only
        * Spearman Ranking variable selection for PLS only

History
-------

* 2009: Developed and used during my PhD thesis at the Laboratory of chemometrics and cheminformatics at the University of Perugia
* 2016: Open-source the code release under GPLv3. The development still continue under open-source licence


Integration with other open-source projects
-------------------------------------------


* Libscientific is the engine of `QStudioMetrics <https://github.com/gmrandazzo/QStudioMetrics>`_, a software to develop easy multivariate analysis
* `MolecularFramework <https://github.com/gmrandazzo/MolecularFramework>`_ a c/c++ cheminformatic library to analyze 3D molecules and develop process 3D information for model development


Licensing
---------

This documentation is copiright (C) 2011-2020 by Giuseppe Marco Randazzo

This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative Commons, 543 Howard Street, 5th Floor, San Francisco, California, 94105, USA.


