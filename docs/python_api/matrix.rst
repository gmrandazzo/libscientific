Matrix
======

.. automodule:: libscientific.matrix
   :no-members:

.. autoclass:: libscientific.matrix.Matrix
   :members:
   :special-members: __getitem__, __setitem__, __len__

.. autodata:: libscientific.matrix.MatrixLike
   :annotation:

Functions
---------

.. autofunction:: libscientific.matrix.init_matrix
.. autofunction:: libscientific.matrix.new_matrix_from_list
.. autofunction:: libscientific.matrix.new_matrix_from_csv
.. autofunction:: libscientific.matrix.new_matrix
.. autofunction:: libscientific.matrix.as_c_matrix
.. autofunction:: libscientific.matrix.resize_matrix
.. autofunction:: libscientific.matrix.del_matrix
.. autofunction:: libscientific.matrix.matrix_check
.. autofunction:: libscientific.matrix.print_matrix
.. autofunction:: libscientific.matrix.val_in_matrix
.. autofunction:: libscientific.matrix.matrix_set
.. autofunction:: libscientific.matrix.matrix_copy
.. autofunction:: libscientific.matrix.set_matrix_value
.. autofunction:: libscientific.matrix.get_matrix_value
.. autofunction:: libscientific.matrix.get_matrix_row
.. autofunction:: libscientific.matrix.get_matrix_column
.. autofunction:: libscientific.matrix.matrix_to_list
.. autofunction:: libscientific.matrix.matrix_to_numpy
.. autofunction:: libscientific.matrix.matrix_from_numpy
.. autofunction:: libscientific.matrix.set_missing_matrix_value
.. autofunction:: libscientific.matrix.matrix_append_row
.. autofunction:: libscientific.matrix.matrix_append_col
.. autofunction:: libscientific.matrix.matrix_dvector_dot_product
.. autofunction:: libscientific.matrix.mt_matrix_dvector_dot_product
.. autofunction:: libscientific.matrix.dvector_matrix_dot_product
.. autofunction:: libscientific.matrix.mt_dvector_matrix_dot_product
.. autofunction:: libscientific.matrix.dvector_transposed_dvector_dot_product
.. autofunction:: libscientific.matrix.dvector_transposed_matrix_division
.. autofunction:: libscientific.matrix.matrix_dot_product
.. autofunction:: libscientific.matrix.row_col_outer_product
.. autofunction:: libscientific.matrix.matrix_transpose
.. autofunction:: libscientific.matrix.matrix_inversion
.. autofunction:: libscientific.matrix.evect_eval
.. autofunction:: libscientific.matrix.svd
.. autofunction:: libscientific.matrix.svd_lapack

Low-level API
-------------

Thin ctypes wrappers around the C functions, working on raw C
pointers. They are needed only to interoperate with other low-level code.

.. autoclass:: libscientific.matrix.MATRIX
   :no-members:
