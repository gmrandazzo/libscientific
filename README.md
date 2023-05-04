libscientific
=============

![Page views](https://visitor-badge.glitch.me/badge?page_id=gmrandazzo.libscientific)
[![status](https://joss.theoj.org/papers/afc8dfc4cdd496f6f51813dbaa5ad310/status.svg)](https://joss.theoj.org/papers/afc8dfc4cdd496f6f51813dbaa5ad310)
[![Licence: GPL v3](https://img.shields.io/github/license/gmrandazzo/libscientific)](https://github.com/gmrandazzo/libscientific/blob/master/LICENSE)
[![Maintainability Rating](https://sonarqube.gmrandazzo.com/api/project_badges/measure?project=libscientific&metric=sqale_rating&token=sqb_426f9a683b14ac981c8e32a1782672fc11c1a789)](https://sonarqube.gmrandazzo.com/dashboard?id=libscientific)
[![Security Rating](https://sonarqube.gmrandazzo.com/api/project_badges/measure?project=libscientific&metric=security_rating&token=sqb_426f9a683b14ac981c8e32a1782672fc11c1a789)](https://sonarqube.gmrandazzo.com/dashboard?id=libscientific)
[![Bugs](https://sonarqube.gmrandazzo.com/api/project_badges/measure?project=libscientific&metric=bugs&token=sqb_426f9a683b14ac981c8e32a1782672fc11c1a789)](https://sonarqube.gmrandazzo.com/dashboard?id=libscientific)
[![Security Hotspots](https://sonarqube.gmrandazzo.com/api/project_badges/measure?project=libscientific&metric=security_hotspots&token=sqb_426f9a683b14ac981c8e32a1782672fc11c1a789)](https://sonarqube.gmrandazzo.com/dashboard?id=libscientific)

Documentation available at http://gmrandazzo.github.io/libscientific/

Libscientific is a C framework for multivariate and other statistical analysis.
It is written in C language using their numerical algorithms for matrix computations.
Only left eigenvectors and eigenvalues are calculated using a third-party library, lapack.

Actually libscientific is able to compute:

  - Multivariate analysis
    - Principal Component Analysis (PCA) NIPALS algorithm) [1]
    - Partial Least Squares (PLS) NIPALS algorithm [1]
    - Consensus PCA (CPCA) NIPALS algorithm [7]
    - Multiple Linear Regression (MLR) Ordinary least squares algorithm
    - Unfold PCA (UPCA) [2]
    - Unfold PLS (UPLS) [2]

  - Pattern recognition
    - Fisher LDA

  - Cluster analysis
    - K-means++ (David Arthur modification) [3]
    - Hierarchical clustering

  - Object selection
    - Most Descriptive Compounds (MDC) [4]
    - Most Dissimilar Compounds  (DIS) [5]

  - Statistical analyisis
    - R2, MSE, MAE, RMSE, BIAS, Sensitivity, Positive Predicted Values 
    - Yates analysis
    - Receiver operating characteristic (ROC)
    - Precision-Recal
    - Matrix-Matrix Euclidean, Manhattan, Cosine and Mahalanobis  distances

  - Numerical analysis
    - Estimate of an integral over a xy region (numerical integration using the trapezoid rule)
    - Natural cubic spline interpolation and prediction
    - Ordinary Least-Squares (OrdinaryLeastSquares)
    - Linear Equation Solver (SolveLSE)
    - Singular value decomposition
  
  - Optimization
     - Nelder-Mead simplex algorithm 

Moreover for some algorithms is possible to run validation methods
with parallel computing to be faster:

- Bootstrap k-fold Cross Validation (RGCV)
- Leave-One-Out
- Y-Scrambling [6]

Usage examples
---------------

* Sampling example on a drug dataset    [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1gp8ppAsGlUbC4qGT-1Frc9ru1PiDvBl1)

* PLS example on the Solubility Dataset [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1eQxLoZOrDMnTkxSkjuyTS_SuAyV3BcSF?usp=sharing)


TODO
----
- Implement Independent Component Analysis ICA
- Exstensive test of some numerical analyisis methods: CholeskyReduction, QR Decomposition, LU Decomposition, HouseholderReduction and so on.

References:
-----------

[1] P. Geladi, B.R. Kowalski
     Partial least-squares regression: a tutorial
     Analytica Chimica Acta Volume 185, 1986, Pages 1-17
     [link](http://dx.doi.org/10.1016/0003-2670(86)80028-9)

[2] S. Wold, P. Geladi, K. Esbensen and J. Öhman
    MULTI-WAY PRINCIPAL COMPONENTSAND PLS-ANALYSIS
    Journal of Chemometrics
    Volume 1, Issue 1, pages 41–56, January 1987
    [link](http://dx.doi.org/10.1002/cem.1180010107)

[3] T. Kanungo, D.M. Mount, N.S. Netanyahu, C.D. Piatko, R. Silverman, A.Y. Wu
    An efficient k-means clustering algorithm: analysis and implementation
    Pattern Analysis and Machine Intelligence, IEEE Transactions on
    Issue Date: Jul 2002
    On page(s): 881 - 892
    [link](http://dx.doi.org/10.1109/TPAMI.2002.1017616)

[4] B.D. Hudson, R.M. Hyde, E. Rahr, J, Wood and J. Osman
    Parameter Based Methods for Compound Selection from Chemical Databases
    Quantitative Structure-Activity Relationships
    Volume 15, Issue 4, pages 285–289, 1996
    [link](http://dx.doi.org/10.1002/qsar.19960150402)

[5] J. Holliday, P. Willett
    Definitions of "Dissimilarity" for Dissimilarity-Based Compound Selection
    Journal of Biomolecular Screening Volume 1, Number 3, 1996 Pages: 145-151
    [link](http://dx.doi.org/10.1177/108705719600100308)

[6] R.D. Clark , P.C. Fox
    Statistical variation in progressive scrambling.
    J Comput Aided Mol Des. 2004 Jul-Sep;18(7-9):563-76.
    [link](http://dx.doi.org/10.1007/s10822-004-4077-z)

[7] J. A. Westerhuis, T. Kourti and J.F. Macgregor
    Analysis of multiblock and hierarchical PCA and PLS models
    Journal of Chemometrics 1998 12, 301-321
    [link](http://dx.doi.org/10.1002/(SICI)1099-128X(199809/10)12:5<301::AID-CEM515>3.0.CO;2-S)


License
============

Libscientific is distributed under GPLv3 license.
To know more in detail how the license work, please read the file "LICENSE" or
go to ["http://www.gnu.org/licenses/gpl-3.0.en.html"](http://www.gnu.org/licenses/gpl-3.0.en.html)


Dependencies
============

The required dependencies to use libscientific are:

- lapack/blas library or a fortran compiler
- c compiler (gcc or clang for osx)
- cmake

Install
=======

Compile from source
-------------------

  cmake -DCMAKE_INSTALL_PREFIX=/usr/local/ ..
  make -j5
  mate test
  sudo make install

Homebrew OSX
------------

brew install --HEAD libscientific


Development
===========
GIT
---

You can check the latest sources with the command::

  git clone https://github.com/gmrandazzo/libscientific.git


Contributing
------------

To contribute, you can fork the project, or if you have already forked the project
update to the latest version of libscientific, make the changes and open a Pull Request.

However, here are some requests.
Before opening a Pull Request:
  * Be sure that your code it's working.
  * No leaks. Run valgrind
  * Comment your code with Parameters, attributes, returns, notes, and References.
  * Test examples are necessary. Tests must prove that
    - the algorithm works correctly
    - the algorithm do not present any memory leak

### How to write a unit test?

Please first read the cmake documentation about [testing with cmake and ctest](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html)

Then write a test for the algorithm you propose and save it in "src/tests" directory.
Run and submit the resulting output in the pull request specifying:
- What the algorithm does
- What the unit tests represent and what they prove.
