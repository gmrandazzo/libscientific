.. libscientific documentation master file, created by
   sphinx-quickstart on Tue Jul 12 10:35:29 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Getting started in C
====================


Every data type object in libscientific is stored in the HEAP and then supports 
dynamic memory allocation.

Hence, every data object such as matrix, vectors, tensors, models, in general,
need to be manually allocated/deallocated by the programmer using the library predefined
constructs "NewSOMETHING(&...);" and "DelSOMETHING(&...);".

Hence, every data object is a pointer and needs to pass by reference "&".
To avoid memory fragmentation problems, please, consider deallocating every allocated variable at the end of your program :-)


Compile a program that use libscientific
----------------------------------------

A program that use libscientific requires only one directive as follow:

.. code-block:: c
   :linenos:

   #include <scientific.h>

and to compile the code with a C or C++ compiler linking with `-lscientific`
and specify the right paths using the `-L<library path/of/libscientific` and `-I<include path/of libscientific>`
options.


.. code-block:: c
   :linenos:

   gcc -o example1 -L/usr/local/lib/ -I/usr/local/include -lscientific example1.c



Vector operations
=================

Create/Allocate a vector
------------------------

There are four different types of vectors

* Double vector: dvector
* Integer vector: ivector
* Unsigned integer vector: uivector
* String vector: strvector


Here we show an example on how to allocate/deallocate these four vector types.

.. literalinclude:: c_code_examples/vexample1.c
   :language: c
   :linenos:


Append a value to a given vector
--------------------------------

Here we show an example on how to append a value to a vector.

.. literalinclude:: c_code_examples/vexample2.c
  :language: c
  :linenos:


Work with string vectors
------------------------



Matrix operations
=================


Matrix is a user-defined data type that contains:
- the number of rows
- the number of columns
- the 2D data array, which defines the matrix.

The data array is selected explicitly as a double type to work with an extensive range of numbers.

.. literalinclude:: ../../src/matrix.h
   :language: c
   :lines: 30-33


Create/Allocate a matrix with a specific size
---------------------------------------------
Create a simple matrix of 10 rows and 15 columns and fill it with numbers.

Then print it to the terminal using "PrintMatrix();"

.. literalinclude:: c_code_examples/mxexample1.c
   :language: c
   :linenos:


Initialize an empty matrix and append a row/column to it
--------------------------------------------------------

An empty matrix is an object with rows and cols equal to 0. However in that
matrix object we can add dynamically rows and columns and or resize it later on.

In this example we will initialize an empty matrix and we will add to it several rows.

.. literalinclude:: c_code_examples/mxexample2.c
   :language: c
   :linenos:


Of course the same code can be reused to add a columns using the construct "MatrixAppendCol" instead of "MatrixAppendRow".


Matrix x Column vector dot product
----------------------------------

In this example we illustrate the product between a matrix of sizes M x N and a column double vector of size N.

.. literalinclude:: c_code_examples/mxexample3.c
   :language: c
   :linenos:


Transpose a matrix
------------------

A matrix transpose is an operation that flips a matrix over its diagonal.
Here is an example that shows how to produce a transpose of a given matrix.

.. literalinclude:: c_code_examples/mxexample4.c
   :language: c
   :linenos:



Invert a matrix
---------------

In this example we show how to invert a matrix with libscientific

.. literalinclude:: c_code_examples/mxexample5.c
   :language: c
   :linenos:

Calculate eigenvectors and eigenvalues of a matrix
--------------------------------------------------

This example shows how to calculate eigenvectors and eigenvalues of an N x N real nonsymmetric matrix.
The eigenvector/eigenvalue is computed thanks to the dgeev.f code extracted from the Lapack library.

.. literalinclude:: c_code_examples/mxexample6.c
  :language: c
  :linenos:

Singular Value Decomposition of a square matrix
-----------------------------------------------

In this example we show how to factorize a square matrix using the singular value decomposition (SVD) method

.. literalinclude:: c_code_examples/mxexample7.c
   :language: c
   :lineos:

Tensor operations
=================

TO BE COMPLETED


Multivariate analysis algorithms
================================

In this section, you will find examples of running multivariate analysis algorithms.
In particular, the algorithm described here is extracted from official scientific publications
and is adapted to run in multithreading to speed up the calculation.

* PCA and PLS implements the NIPALS algorithm described in the following publication:

| P. Geladi, B.R. Kowalski
| Partial least-squares regression: a tutorial
| Analytica Chimica Acta Volume 185, 1986, Pages 1-17
| DOI:10.1016/0003-2670(86)80028-9



Principal Component Analysis (PCA)
----------------------------------

Here is an example that shows how to compute a principal component analysis on a matrix.


.. literalinclude:: c_code_examples/pcaexample1.c
   :language: c
   :linenos:


Partial Least Squares (PLS)
---------------------------

A matrix of features or independent vIariables and a matrix of targets or dependent variables is requested to calculate a PLS model.
Here is a simple example that shows how to calculate a PLS model.


.. literalinclude:: c_code_examples/plsexample1.c
   :language: c
   :linenos:
