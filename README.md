libscientific
=============

Libscientific is a C framework  for multivariate and other statistical analysis.
It is purely writen in C language using their own numerical algorithms
for matrix computations. Only left eigenvectors and eigenvalues are calculated
using a third part library, lapack version 3.5.0.

Actually libscientific is able to compute:

  -Multivariate analysis
    - PCA (NIPALS ALGORITHM) [1]
    - PLS (NIPALS ALGORITHM) [1]
    - MLR (Ordinary least squares algorithm)
    - UPCA [2]
    - UPLS [2]

  -Pattern recognition
    - Fisher LDA

  -Cluster analysis
    - K-means++ (David Arthur modification) [3]
    - Hierarchical clustering

  -Object selection
    - Most Descriptive Compounds (MDC) [4]
    - Most Dissimilar Compounds  (DIS) [5]

Moreover for some algorithms is possible to run validation methods
with parallel computing to be faster:

- Bootstrap k-fold Cross Validation (RGCV)
- Leave-One-Out
- Y-Scrambling [6]


References:

[1] P. Geladi, B.R. Kowalski
     Partial least-squares regression: a tutorial
     Analytica Chimica Acta Volume 185, 1986, Pages 1–17
     DOI:10.1016/0003-2670(86)80028-9

[2] S. Wold, P. Geladi, K. Esbensen and J. Öhman3
    MULTI-WAY PRINCIPAL COMPONENTSAND PLS-ANALYSIS
    Journal of Chemometrics
    Volume 1, Issue 1, pages 41–56, January 1987
    DOI: 10.1002/cem.1180010107

[3] T. Kanungo, D.M. Mount, N.S. Netanyahu, C.D. Piatko, R. Silverman, A.Y. Wu
    An efficient k-means clustering algorithm: analysis and implementation
    Pattern Analysis and Machine Intelligence, IEEE Transactions on
    Issue Date: Jul 2002
    On page(s): 881 - 892
    DOI:10.1109/TPAMI.2002.1017616

[4] B.D. Hudson, R.M. Hyde, E. Rahr, J, Wood and J. Osman
    Parameter Based Methods for Compound Selection from Chemical Databases
    Quantitative Structure-Activity Relationships
    Volume 15, Issue 4, pages 285–289, 1996
    DOI: 10.1002/qsar.19960150402

[5] J. Holliday, P. Willett
    Definitions of "Dissimilarity" for Dissimilarity-Based Compound Selection
    Journal of Biomolecular Screening Volume 1, Number 3, 1996 Pages: 145-151
    DOI:10.1177/108705719600100308

[6] R.D. Clark , P.C. Fox
    Statistical variation in progressive scrambling.
    J Comput Aided Mol Des. 2004 Jul-Sep;18(7-9):563-76.
    DOI 10.1007/s10822-004-4077-z

License
============

Libscientific is distributed under GPLv3 license, this means that:

- you can use this library where you want doing what you want.
- you can modify this library and commit changes.
- you can not use this library inside a commercial software.

To know more in details how the licens work please read the file "LICENSE" or
go to "http://www.gnu.org/licenses/gpl-3.0.en.html"

Libscientific is currently property of Giuseppe Marco Randazzo which is also the
current package maintainer.

Voluntary contributions are welcome. :-)


Dependencies
============

The required dependencies to use libscientific are:

- fortran compiler
- c compiler (gcc or clang for osx)
- cmake

Install
=======

To install for all users on Unix/Linux/OSX/Windows:

  cmake .. -DPREFIX=/usr/local/  
  make -j5
  sudo make install


Development
===========
GIT
~~~

You can check the latest sources with the command::

  git clone https://github.com/zeld/libscientific.git

~~~

Contributing
~~~~~~~~~~~~~

To contribute you can fork the project, or if you have already forked the project
update to the latest version of libscientific, make the changes and open a Pull Request.

However some recommendations before open a Pull Request:
  * Be sure that your code it's working.
  * Comment your code with Parameters, Attribute, Return, Notes and References.
  * A test example is necessary.

Probably your code will be integrated but some quality controls and goals
have to be respected.
