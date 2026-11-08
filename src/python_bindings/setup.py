import platform

from setuptools import Distribution, setup

try:  # setuptools >= 70.1 ships bdist_wheel
    from setuptools.command.bdist_wheel import bdist_wheel as _bdist_wheel
except ImportError:  # older setuptools: provided by the wheel package
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel

if platform.system() == "Windows":
    library = "*.dll"
elif platform.system() == "Darwin":
    library = "*.dylib"
else:
    library = "*.so*"


class bdist_wheel(_bdist_wheel):
    """
    The package bundles a native shared library loaded through ctypes:
    the wheel is platform specific but independent of the Python version
    and ABI, so tag it py3-none-<platform> (e.g. py3-none-win_amd64).
    """

    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False

    def get_tag(self):
        _, _, plat = super().get_tag()
        return "py3", "none", plat


class BinaryDistribution(Distribution):
    """The package ships a prebuilt shared library: install it as platlib."""

    def has_ext_modules(self):
        return True


setup(
    distclass=BinaryDistribution,
    packages=["libscientific"],
    package_data={"libscientific": [library]},
    cmdclass={"bdist_wheel": bdist_wheel},
)
