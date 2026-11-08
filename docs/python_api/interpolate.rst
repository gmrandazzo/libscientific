Interpolation
=============

.. automodule:: libscientific.interpolate
   :no-members:

.. autoclass:: libscientific.interpolate.CubicSplineInterpolation
   :members:

Low-level API
-------------

Thin ctypes wrappers around the C functions, working on raw C
pointers. They are needed only to interoperate with other low-level code.

.. autofunction:: libscientific.interpolate.cubic_spline_interpolation
.. autofunction:: libscientific.interpolate.cubic_spline_predict
