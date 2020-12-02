.. libscientific documentation master file, created by
   sphinx-quickstart on Tue Jul 12 10:35:29 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Getting started in C
====================

Every data type object in libscientific is stored in the HEAP and then supports
dynamic memory allocation.

This means that every data object such as, matrix, vectors, tensors, models in general need to be manually allocated/deallocated
by the programmer using the library predefined constructs "NewSOMETHING(&...);" and "DelSOMETHING(&...);".

This also means that every data object is a pointer and need to pass by reference "&".

To avoid memory fragmentation problems please, consider that every allocated variable
need to be deallocated at the end of your program :-)


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


Matrix operations
=================

Matrix is an user defined data type which contains informations in regards to
- the number of rows
- the number of columns
- the 2D data array which define the matrix

The data array is specifically selected to be double type
to work with large range of numbers.

.. literalinclude:: ../src/matrix.h
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

A matrix transpose is an operation which flip a matrix over it's diagonal.
Here an example that shows how to produce a transpose of a given matrix.

.. literalinclude:: c_code_examples/mxexample4.c
   :language: c
   :linenos:



Multivariate analysis algorithms
================================

In this section you will find examples regarding how to run multivariate analysis algorithms.
In particular the algorithm described here are extracted from official scientific publications
and are adapted to run in multithreading to speedup the calculation.

* PCA and PLS implements the NIPALS algorithm described in the following publication:

| P. Geladi, B.R. Kowalski
| Partial least-squares regression: a tutorial
| Analytica Chimica Acta Volume 185, 1986, Pages 1–17
| DOI:10.1016/0003-2670(86)80028-9



Principal Component Analysis (PCA)
----------------------------------

Here an example to shows how to compute a principal component analysis on a matrix.


.. literalinclude:: c_code_examples/pcaexample1.c
   :language: c
   :linenos:


Partial Least Squares (PLS)
---------------------------
To calculate a PLS model, a matrix of features or independent variables and a matrix of targets or dependent variables is requested.
Here a simple example that shows how to calculate a PLS model.


.. literalinclude:: c_code_examples/plsexample1.c
   :language: c
   :linenos:
