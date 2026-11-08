Tensor
======

.. automodule:: libscientific.tensor
   :no-members:

.. autoclass:: libscientific.tensor.Tensor
   :members:
   :special-members: __getitem__, __setitem__, __len__

Functions
---------

.. autofunction:: libscientific.tensor.init_tensor
.. autofunction:: libscientific.tensor.new_tensor_matrix
.. autofunction:: libscientific.tensor.add_tensor_matrix
.. autofunction:: libscientific.tensor.new_tensor
.. autofunction:: libscientific.tensor.del_tensor
.. autofunction:: libscientific.tensor.print_tensor
.. autofunction:: libscientific.tensor.set_tensor_value
.. autofunction:: libscientific.tensor.tensor_append_column
.. autofunction:: libscientific.tensor.tensor_append_row
.. autofunction:: libscientific.tensor.tensor_set
.. autofunction:: libscientific.tensor.tensor_copy
.. autofunction:: libscientific.tensor.tensor_tolist

Low-level API
-------------

Thin ctypes wrappers around the C functions, working on raw C
pointers. They are needed only to interoperate with other low-level code.

.. autoclass:: libscientific.tensor.TENSOR
   :no-members:
