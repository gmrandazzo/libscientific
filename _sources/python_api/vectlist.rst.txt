List of vectors
===============

.. automodule:: libscientific.vectlist
   :no-members:

.. autoclass:: libscientific.vectlist.DVectorList
   :members:
   :special-members: __getitem__, __setitem__, __len__

Functions
---------

.. autofunction:: libscientific.vectlist.init_dvector_list
.. autofunction:: libscientific.vectlist.new_dvector_list
.. autofunction:: libscientific.vectlist.del_dvector_list
.. autofunction:: libscientific.vectlist.dvector_list_append
.. autofunction:: libscientific.vectlist.dvector_list_tolist

Low-level API
-------------

Thin ctypes wrappers around the C functions, working on raw C
pointers. They are needed only to interoperate with other low-level code.

.. autoclass:: libscientific.vectlist.DVECTLIST
   :no-members:
