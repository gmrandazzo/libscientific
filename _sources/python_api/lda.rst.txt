LDA — Linear Discriminant Analysis
==================================

.. automodule:: libscientific.lda
   :no-members:

.. autoclass:: libscientific.lda.LDA
   :members:

Low-level API
-------------

Thin ctypes wrappers around the C functions, working on raw C
pointers. They are needed only to interoperate with other low-level code.

.. autoclass:: libscientific.lda.LDAMODEL
   :no-members:

.. autofunction:: libscientific.lda.new_lda_model
.. autofunction:: libscientific.lda.del_lda_model
.. autofunction:: libscientific.lda.lda_algorithm
.. autofunction:: libscientific.lda.lda_prediction
.. autofunction:: libscientific.lda.write_lda
.. autofunction:: libscientific.lda.read_lda
