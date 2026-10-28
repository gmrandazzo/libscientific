import platform

from setuptools import setup

if platform.system() == "Windows":
    library = "*.dll"
elif platform.system() == "Darwin":
    library = "*.dylib"
else:
    library = "*.so*"

arch = platform.machine()

setup(
    name="libscientific",
    version="1.8.0",
    packages=["libscientific"],
    package_data={"": [library]},
    platforms=[arch],
)
