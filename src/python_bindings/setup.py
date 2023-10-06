from setuptools import setup
import platform

if platform.system() == 'Windows':
    library = 'libscientific.dll'
elif platform.system() == 'Darwin':
    library = 'libscientific.dylib'
else:
    library = 'libscientific.so'

arch = platform.machine()

setup(
    name="libscientific",
    version="1.6.0",
    packages=["libscientific"],
    package_data={"": [library]},
    platforms=[arch],
)

