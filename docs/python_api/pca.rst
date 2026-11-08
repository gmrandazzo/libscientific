PCA — Principal Component Analysis
==================================

.. automodule:: libscientific.pca
   :no-members:

.. autoclass:: libscientific.pca.PCA
   :members:

Low-level API
-------------

Thin ctypes wrappers around the C functions, working on raw C
pointers. They are needed only to interoperate with other low-level code.

.. autoclass:: libscientific.pca.PCAMODEL
   :no-members:

.. autofunction:: libscientific.pca.new_pca_model
.. autofunction:: libscientific.pca.del_pca_model
.. autofunction:: libscientific.pca.pca_algorithm
.. autofunction:: libscientific.pca.pca_score_predictor
.. autofunction:: libscientific.pca.pca_ind_var_predictor
.. autofunction:: libscientific.pca.print_pca
.. autofunction:: libscientific.pca.write_pca
.. autofunction:: libscientific.pca.read_pca
