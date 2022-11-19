from setuptools import setup
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='libscientific',
      version='1.4.0',
      description='Libscientific python foreign function interface',
      long_description=long_description,
      long_description_content_type='text/markdown',
      url='http://github.com/gmrandazzo/libscientific',
      author='Giuseppe Marco Randazzo',
      author_email='gmrandazzo@gmail.com',
      license='GPLv3',
      packages=['libscientific'],
      python_requires=">=3.0",
      zip_safe=False)
