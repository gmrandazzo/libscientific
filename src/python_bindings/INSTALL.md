Dependencies
============

The required dependencies to use libscientific are:

- fortran compiler
- c compiler (gcc or clang for osx)
- python
- cmake

Install
=======

Compile from source
-------------------
  cd libscientific
  mkdir build
  cd build
  cmake -DCMAKE_INSTALL_PREFIX=/usr/local/ ..
  make -j5
  sudo make install
  cd ../src/python_bindings
  python3 setup.py install --user 


Homebrew OSX
------------

brew install --HEAD libscientific
pip3 install libscientific


Development
===========
GIT
---

You can check the latest sources with the command::

  git clone https://github.com/gmrandazzo/libscientific.git


Contributing
------------

To contribute you can fork the project, or if you have already forked the project
update to the latest version of libscientific, make the changes and open a Pull Request.

However some recommendations before open a Pull Request:
  * Be sure that your code it's working.
  * Comment your code with Parameters, Attribute, Return, Notes and References.
  * A test example is necessary.

Probably your code will be integrated but some quality controls and goals
have to be respected.
