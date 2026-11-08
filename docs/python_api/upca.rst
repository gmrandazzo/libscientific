UPCA — Multi-way PCA
====================

.. automodule:: libscientific.upca
   :no-members:

.. autoclass:: libscientific.upca.UPCA
   :members:

Functions
---------

.. autofunction:: libscientific.upca.as_c_tensor

Low-level API
-------------

Thin ctypes wrappers around the C functions, working on raw C
pointers. They are needed only to interoperate with other low-level code.

.. autoclass:: libscientific.upca.UPCAMODEL
   :no-members:

.. autofunction:: libscientific.upca.new_upca_model
.. autofunction:: libscientific.upca.del_upca_model
.. autofunction:: libscientific.upca.upca_algorithm
.. autofunction:: libscientific.upca.upca_score_predictor
.. autofunction:: libscientific.upca.upca_ind_var_predictor
.. autofunction:: libscientific.upca.write_upca
.. autofunction:: libscientific.upca.read_upca
