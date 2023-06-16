.. libscientific documentation master file, created by
   sphinx-quickstart on Tue Jul 12 10:35:29 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Install
=======

The installation works for Linux/macOS and Windows.

Requirements
------------

* A development environment (On windows msys/msys2 or visual studio)
* A c/fortran compiler
* Cmake
* python3 (if you whant to use the library in python)


Installation process
--------------------

First you need to install the C library following these instructions:

.. code-block::
   
   git clone https://github.com/gmrandazzo/libscientific.git
   cd libscientific
   mkdir build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=/usr/local/ ..
   make -j5
   sudo make install

Then, if you want to use the library in python, you have also to install the python package.

.. code-block::
   
   cd ../src/python_bindings
   python3 setup.py install

   or 

   pip3 install libscientific


Packages
--------
On macOS you can install libscientific via homebrew

.. code-block::
   
   brew install --HEAD libscientific
   brew install --HEAD libscientific-python3




