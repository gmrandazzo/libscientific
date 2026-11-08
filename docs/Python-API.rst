********************
Python API Reference
********************

Every model follows the same pattern::

    model = PCA(scaling=1, npc=2)     # 1. create with the options
    model.fit(X)                      # 2. fit on training data
    scores = model.predict(X_new)     # 3. use on new data
    model.save("model.lsm")           # 4. save / load

Inputs can be Python lists of rows, numpy arrays or
:class:`~libscientific.matrix.Matrix` objects (see
:data:`~libscientific.matrix.MatrixLike`). Results are Python lists: the
documentation of every method gives their exact shape. Wrap them with
``numpy.array(...)`` if you need arrays.

The examples in the documentation are executed by the test suite, so they
always match the behaviour of the library.

Models
======

.. toctree::
   :maxdepth: 1

   python_api/pca
   python_api/pls
   python_api/cpca
   python_api/mlr
   python_api/lda
   python_api/upca
   python_api/upls
   python_api/clustering
   python_api/interpolate
   python_api/model_io

Data structures and utilities
=============================

.. toctree::
   :maxdepth: 1

   python_api/matrix
   python_api/vector
   python_api/tensor
   python_api/vectlist
   python_api/misc
   python_api/info
   python_api/loadlibrary
