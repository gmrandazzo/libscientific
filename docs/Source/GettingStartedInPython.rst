.. libscientific documentation master file, created by
   sphinx-quickstart on Tue Jul 12 10:35:29 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Getting started in Python
=========================

Every data type object in libscientific is stored in the HEAP and then supports
dynamic memory allocation.

In python, there is no need to allocate/deallocate matrix/vectors/tensors and models in general
because are automatically handled by the python binding itself.


Use libscientific in python
---------------------------
First you need to install the python package. For that follow the installation process described 
here.

A program that use libscientific requires to import the python binding as follow

.. code-block:: python
   :linenos:

   import scientific
   ...


Vector operations
=================

Create a vector in python
-------------------------

There are four different types of vectors

* Double vector: dvector
* Integer vector: ivector
* Unsigned integer vector: uivector
* String vector: strvector


Here we show an example on how create these four vector types.

.. literalinclude:: python_code_examples/vexample1.py
   :language: c
   :linenos:


Append a value to a given vector
--------------------------------

Here we show an example on how to append a value to a vector.

.. literalinclude:: python_code_examples/vexample2.py
  :language: c
  :linenos:



Matrix operations
=================

Matrix is an user defined data type which contains informations in regards to
- the number of rows
- the number of columns
- the 2D data array which define the matrix

The data array in python use the same implementation of the c language version.
However memory allocation/destruction are carried out directly from the python class.
Hence there is no need to manually free up the memory.


Create a matrix in python
--------------------------

In this example we show how to create a matrix from a list of list (or numpy array)
and we show how to modify its content and convert it again to a list of list.

.. literalinclude:: python_code_examples/mxexample1.py
   :language: c
   :linenos:



Tensor operations
=================

TO BE COMPLETED


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


.. literalinclude:: python_code_examples/pcaexample1.py
   :language: c
   :linenos:


Partial Least Squares (PLS)
---------------------------

To calculate a PLS model, a matrix of features or independent variables and a matrix of targets or dependent variables is requested.

Here a simple example that shows how to calculate a PLS model.

.. literalinclude:: python_code_examples/plsexample1.py
   :language: c
   :linenos:


