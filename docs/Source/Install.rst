.. libscientific documentation master file, created by
   sphinx-quickstart on Tue Jul 12 10:35:29 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Install
=======

The installation works for Linux/OSX and Windows.


Requirements
------------

* A development environment (On windows msys/msys2 or visual studio)
* A c/fortran compiler
* Cmake


Installation process
--------------------

.. code-block::
   
   git clone https://github.com/gmrandazzo/libscientific.git
   cd libscientific
   mkdir build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=/usr/local/ ..
   make -j5
   sudo make install


Packages
--------
On OSX you can install libscientific via homebrew

.. code-block::
   
   brew install --HEAD libscientific





