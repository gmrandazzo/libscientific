Vectors
=======

.. automodule:: libscientific.vector
   :no-members:

.. autoclass:: libscientific.vector.DVector
   :members:
   :special-members: __getitem__, __setitem__, __len__

.. autoclass:: libscientific.vector.UIVector
   :members:
   :special-members: __getitem__, __setitem__, __len__

Functions
---------

.. autofunction:: libscientific.vector.xvector_tolist
.. autofunction:: libscientific.vector.init_dvector
.. autofunction:: libscientific.vector.new_dvector
.. autofunction:: libscientific.vector.del_dvector
.. autofunction:: libscientific.vector.dvector_resize
.. autofunction:: libscientific.vector.print_dvector
.. autofunction:: libscientific.vector.dvector_has_value
.. autofunction:: libscientific.vector.dvector_set
.. autofunction:: libscientific.vector.set_dvector_value
.. autofunction:: libscientific.vector.get_dvector_value
.. autofunction:: libscientific.vector.dvector_append
.. autofunction:: libscientific.vector.dvector_tolist
.. autofunction:: libscientific.vector.dvector_remove_at
.. autofunction:: libscientific.vector.dvector_copy
.. autofunction:: libscientific.vector.init_uivector
.. autofunction:: libscientific.vector.new_uivector
.. autofunction:: libscientific.vector.del_uivector
.. autofunction:: libscientific.vector.uivector_resize
.. autofunction:: libscientific.vector.print_uivector
.. autofunction:: libscientific.vector.set_uivector_value
.. autofunction:: libscientific.vector.get_uivector_value
.. autofunction:: libscientific.vector.uivector_append
.. autofunction:: libscientific.vector.uivector_tolist
.. autofunction:: libscientific.vector.uivector_remove_at
.. autofunction:: libscientific.vector.uivector_copy
.. autofunction:: libscientific.vector.init_ivector
.. autofunction:: libscientific.vector.new_ivector
.. autofunction:: libscientific.vector.del_ivector
.. autofunction:: libscientific.vector.set_ivector_value
.. autofunction:: libscientific.vector.get_ivector_value
.. autofunction:: libscientific.vector.ivector_tolist
.. autofunction:: libscientific.vector.ivector_append
.. autofunction:: libscientific.vector.ivector_remove_at

Low-level API
-------------

Thin ctypes wrappers around the C functions, working on raw C
pointers. They are needed only to interoperate with other low-level code.

.. autoclass:: libscientific.vector.STRVECTOR
   :no-members:

.. autoclass:: libscientific.vector.DVECTOR
   :no-members:

.. autoclass:: libscientific.vector.UIVECTOR
   :no-members:

.. autoclass:: libscientific.vector.IVECTOR
   :no-members:
