from setuptools import setup
import platform

if platform.system() == 'Windows':
    library = '*.dll'
elif platform.system() == 'Darwin':
    library = '*.dylib'
else:
    library = '*.so*'

arch = platform.machine()

setup(
    name="libscientific",
    version="1.6.0",
    packages=["libscientific"],
    package_data={"": [library]},
    platforms=[arch],
)

