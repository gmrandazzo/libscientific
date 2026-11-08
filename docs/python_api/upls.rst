UPLS — Multi-way PLS
====================

.. automodule:: libscientific.upls
   :no-members:

.. autoclass:: libscientific.upls.UPLS
   :members:

Low-level API
-------------

Thin ctypes wrappers around the C functions, working on raw C
pointers. They are needed only to interoperate with other low-level code.

.. autoclass:: libscientific.upls.UPLSMODEL
   :no-members:

.. autofunction:: libscientific.upls.new_upls_model
.. autofunction:: libscientific.upls.del_upls_model
.. autofunction:: libscientific.upls.upls_algorithm
.. autofunction:: libscientific.upls.upls_score_predictor
.. autofunction:: libscientific.upls.upls_y_predictor
.. autofunction:: libscientific.upls.write_upls
.. autofunction:: libscientific.upls.read_upls
