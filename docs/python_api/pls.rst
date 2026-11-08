PLS — Partial Least Squares regression
======================================

.. automodule:: libscientific.pls
   :no-members:

.. autoclass:: libscientific.pls.PLS
   :members:

Functions
---------

.. autofunction:: libscientific.pls.apply_scaling
.. autofunction:: libscientific.pls.revert_scaling
.. autofunction:: libscientific.pls.pls_beta_inference

Low-level API
-------------

Thin ctypes wrappers around the C functions, working on raw C
pointers. They are needed only to interoperate with other low-level code.

.. autoclass:: libscientific.pls.PLSMODEL
   :no-members:

.. autofunction:: libscientific.pls.new_pls_model
.. autofunction:: libscientific.pls.del_pls_model
.. autofunction:: libscientific.pls.pls_algorithm
.. autofunction:: libscientific.pls.pls_beta_coefficients
.. autofunction:: libscientific.pls.pls_vip
.. autofunction:: libscientific.pls.pls_score_predictor
.. autofunction:: libscientific.pls.pls_y_predictor
.. autofunction:: libscientific.pls.pls_y_predictor_all_lv
.. autofunction:: libscientific.pls.print_pls
.. autofunction:: libscientific.pls.write_pls
.. autofunction:: libscientific.pls.read_pls
