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
* sqlite3 library
* Cmake
* python3 (if you whant to use the library in python)


Installation process
--------------------

First you need to install the C library following these instructions:

.. code-block::

   sudo apt-get -y install liblapack-dev cmake libsqlite3-dev
   git clone https://github.com/gmrandazzo/libscientific.git
   cd libscientific
   mkdir build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=/usr/ ..
   make -j5
   sudo make install

Then, if you want to use the library in python, you have also to install the python package.

.. code-block::
   
   cd ../src/python_bindings
   python3 setup.py install

   or 

   pip3 install libscientific


If you have installed the library in a non-standard path and you want to use the python bindings
please export the LIBSCIENTIFIC_LIB_DIR variable pointing to the libscientific.so installation path.

For example, if the library is installed in /home/<user>/mylocalenv, then you need to export the 
following variable

.. code-block::

    export LISCIENTIFIC_LIB_DIR=/home/<user>/mylocalenv


Packages
--------
On macOS you can install libscientific via homebrew

.. code-block::
   
   brew install --HEAD libscientific
   brew install --HEAD libscientific-python3




