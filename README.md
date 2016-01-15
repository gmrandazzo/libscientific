libscientific
=============

Libscientific is a C framework  for multivariate and other statistical analysis.
It is purely writen in C language using their own numerical algorithms
for matrix computations. Only left eigenvectors and eigenvalues are calculated
using a third part library, lapack version 3.5.0.

Actually libscientific is able to compute:
- Multivariate analysis
  * PCA (NIPALS ALGORITHM)
  * PLS (NIPALS ALGORITHM)
  * MLR (Ordinary least squares algorithm)
  * UPCA
  * UPLS

- Pattern recognition
  * Fisher LDA

- Cluster analysis
  * K-means++ (David Arthur algorithm)
  * Hierarchical clustering

- Object selection
  * Most Descriptive Compounds (MDC)
  * Most Dissimilar Compounds  (DIS)

Moreover for some algorithms this library contains validation methods such as:
- Bootstrap k-fold Cross Validation (RGCV)
- Leave-One-Out
- Y-Scrambling
running in with parallel computing to be faster.


Libscientific is distributed under GPLv3 license, this means that:

- you can use this library where you want doing what you want.
- you can modify this library and commit changes.
- you can not use this library inside a commercial software.

License
============

Libscientific is currently property of Giuseppe Marco Randazzo which is also the
current package maintainer.

Voluntary contributions are welcome. :-)


Dependencies
============

The required dependencies to use libscientific are:

- fortran compiler
- c compiler (gcc or clang)
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
