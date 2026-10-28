.. libscientific documentation master file, created by
   sphinx-quickstart on Tue Jul 12 10:35:29 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Getting started in Python
=========================

Every data type object in libscientific is stored in the HEAP and then supports dynamic memory allocation.

In python, there is no need to allocate/deallocate matrix/vectors/tensors and models in general
because the python binding itself automatically handles them.


Use libscientific in python
---------------------------
First, you need to install the c library and the python package.
Please follow the process described `here <http://gmrandazzo.github.io/libscientific/Source/_build/html/Install.html>`_.

A program that use libscientific requires to import the python binding as follows

.. code-block:: python
   :linenos:

   import libscientific
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
   :language: python
   :linenos:


Append a value to a given vector
--------------------------------

Here we show an example on how to append a value to a vector.

.. literalinclude:: python_code_examples/vexample2.py
  :language: python
  :linenos:



Matrix operations
=================

Matrix is a user-defined data type that contains information in regards to
- the number of rows
- the number of columns
- the 2D data array which defines the matrix

The data array in python uses the c language implementation.
However, memory allocation/destruction is carried out directly from the python class.
Hence there is no need to free up the memory manually.


Create a matrix in python
--------------------------

In this example, we show how to create a matrix from a list of lists (or numpy array),
modify its content and convert it again to a list of lists.

.. literalinclude:: python_code_examples/mxexample1.py
   :language: python
   :linenos:

Tensor operations
=================

Tensor is a user-defined data type that contains:
- order: the number of matrix
- m: the array the 2D data array, which defines the tensor itself.

The data array in python uses the c language implementation.
However, memory allocation/destruction is carried out directly from the python class.
Hence there is no need to free up the memory manually


Create a tensor in python
--------------------------

In this example, we show how to create a tensor from a list of list of lists (or numpy array),
modify its content and convert it again to a list of lists.

.. literalinclude:: python_code_examples/tnsexample1.py
   :language: python
   :linenos:

Multivariate analysis algorithms
================================

In this section, you will find examples of running multivariate analysis algorithms.
In particular, the algorithm described here is extracted from official libscientific publications and is adapted to run in multithreading to speed up the calculation.


* PCA and PLS implements the NIPALS algorithm described in the following publication:

| P. Geladi, B.R. Kowalski
| Partial least-squares regression: a tutorial
| Analytica Chimica Acta Volume 185, 1986, Pages 1–17
| DOI:10.1016/0003-2670(86)80028-9

* CPCA implements the NIPALS algorithm described in the following publication:

| ANALYSIS OF MULTIBLOCK AND HIERARCHICAL PCA AND PLS MODELS
| JOHAN A. WESTERHUIS, THEODORA KOURTI* AND JOHN F. MACGREGOR
| J. Chemometrics 12, 301–321 (1998)
| DOI:/10.1002/(SICI)1099-128X(199809/10)12:5<301::AID-CEM515>3.0.CO;2-S


Principal Component Analysis (PCA)
----------------------------------


Here is an example that shows to compute a principal component analysis on a matrix.


.. literalinclude:: python_code_examples/pcaexample1.py
   :language: python
   :linenos:

Consensus Principal Component Analysis (CPCA)
---------------------------------------------

Here is an example that shows how to compute a consenus principal component analysis on a tensor.


.. literalinclude:: python_code_examples/cpcaexample1.py
   :language: python
   :linenos:



Partial Least Squares (PLS)
---------------------------

A matrix of features or independent variables and a matrix of targets or dependent variables is requested to calculate a PLS model.

Here is a simple example that shows how to calculate a PLS model.

.. literalinclude:: python_code_examples/plsexample1.py
   :language: python
   :linenos:


Model persistence and serialization (.lsm and .lsb)
===================================================

Libscientific allows saving and loading models in Python with zero external dependencies.

Saving and loading a single model (.lsm)
----------------------------------------

.. code-block:: python
   :linenos:

   from libscientific.pca import PCA
   from libscientific.model_io import save_pca_model, load_pca_model

   # Train PCA model
   model = PCA(matrix, n_components=3)

   # Save model to .lsm container
   save_pca_model("pca_model.lsm", model)

   # Load model back from disk
   loaded_model = load_pca_model("pca_model.lsm")


Working with multi-model bundles (.lsb)
---------------------------------------

Bundles allow grouping multiple related models together into a single file:

.. code-block:: python
   :linenos:

   from libscientific.model_io import ModelBundle

   # Create and populate a project bundle
   bundle = ModelBundle()
   bundle.add_pca("screening_pca", pca_model)
   bundle.add_pls("activity_pls", pls_model)

   # Save bundle to disk
   bundle.save("project.lsb")

   # Load bundle back
   loaded_bundle = ModelBundle.load("project.lsb")
   retrieved_pca = loaded_bundle.get_pca("screening_pca")
   retrieved_pls = loaded_bundle.get_pls("activity_pls")
