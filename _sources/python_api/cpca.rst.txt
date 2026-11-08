CPCA — Consensus (multi-block) PCA
==================================

.. automodule:: libscientific.cpca
   :no-members:

.. autoclass:: libscientific.cpca.CPCA
   :members:

Low-level API
-------------

Thin ctypes wrappers around the C functions, working on raw C
pointers. They are needed only to interoperate with other low-level code.

.. autoclass:: libscientific.cpca.CPCAMODEL
   :no-members:

.. autofunction:: libscientific.cpca.new_cpca_model
.. autofunction:: libscientific.cpca.del_cpca_model
.. autofunction:: libscientific.cpca.cpca_algorithm
.. autofunction:: libscientific.cpca.cpca_score_predictor
.. autofunction:: libscientific.cpca.print_cpca
.. autofunction:: libscientific.cpca.write_cpca
.. autofunction:: libscientific.cpca.read_cpca
