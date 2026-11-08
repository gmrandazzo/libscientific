Model files and bundles
=======================

.. automodule:: libscientific.model_io
   :no-members:

.. autoclass:: libscientific.model_io.ModelType
   :members:

.. autoclass:: libscientific.model_io.ModelBundle
   :members:

Functions
---------

.. autofunction:: libscientific.model_io.inspect_model_file
.. autofunction:: libscientific.model_io.is_bundle
.. autofunction:: libscientific.model_io.save_model
.. autofunction:: libscientific.model_io.import_model
.. autofunction:: libscientific.model_io.save_model_bundle
.. autofunction:: libscientific.model_io.import_model_bundle
.. autofunction:: libscientific.model_io.import_auto

Low-level API
-------------

Thin ctypes wrappers around the C functions, working on raw C
pointers. They are needed only to interoperate with other low-level code.

.. autoclass:: libscientific.model_io.MODEL_FILE_INFO
   :no-members:

.. autoclass:: libscientific.model_io.MODEL_ENTRY
   :no-members:

.. autoclass:: libscientific.model_io.MODEL_BUNDLE
   :no-members:
