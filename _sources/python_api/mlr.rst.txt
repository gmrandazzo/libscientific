MLR — Multiple Linear Regression
================================

.. automodule:: libscientific.mlr
   :no-members:

.. autoclass:: libscientific.mlr.MLR
   :members:

Low-level API
-------------

Thin ctypes wrappers around the C functions, working on raw C
pointers. They are needed only to interoperate with other low-level code.

.. autoclass:: libscientific.mlr.MLRMODEL
   :no-members:

.. autofunction:: libscientific.mlr.new_mlr_model
.. autofunction:: libscientific.mlr.del_mlr_model
.. autofunction:: libscientific.mlr.mlr_algorithm
.. autofunction:: libscientific.mlr.mlr_predict_y
.. autofunction:: libscientific.mlr.write_mlr
.. autofunction:: libscientific.mlr.read_mlr
