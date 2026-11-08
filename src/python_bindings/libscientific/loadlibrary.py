# loadlibrary.py method to load libscientific library
#
# Copyright (C) <2023>  Giuseppe Marco Randazzo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""Locate and load the libscientific C library.

Used internally at import time: the library bundled in the wheel is loaded
first; otherwise it is searched in ``LD_LIBRARY_PATH`` /
``DYLD_LIBRARY_PATH`` and the system library directories.
"""

import ctypes
import os
import pathlib
import platform


def get_platform_info() -> tuple[str, str]:
    """Return the operating system and the CPU architecture.

    Returns
    -------
    tuple[str, str]
        ``(system, machine)`` in lower case, e.g. ``("linux", "x86_64")`` or
        ``("darwin", "arm64")``; ``amd64`` is reported as ``x86_64`` and
        ``aarch64`` as ``arm64``.
    """
    system = platform.system().lower()
    machine = platform.machine().lower()
    arch_map = {"x86_64": "x86_64", "amd64": "x86_64", "arm64": "arm64", "aarch64": "arm64"}
    arch = arch_map.get(machine, machine)
    return system, arch


def get_platform_specific_paths() -> list[pathlib.Path]:
    """Return the standard library directories of the current platform.

    Returns
    -------
    list[pathlib.Path]
        Directories such as ``/usr/lib``, ``/usr/local/lib``, Homebrew prefixes
        on macOS, ...
    """
    system, arch = get_platform_info()
    paths = [
        pathlib.Path("/usr/lib"),
        pathlib.Path("/usr/local/lib"),
    ]
    if system == "linux":
        # Get machine architecture
        if arch == "x86_64":
            paths.extend(
                [
                    pathlib.Path("/usr/lib/x86_64-linux-gnu"),  # Debian/Ubuntu
                    pathlib.Path("/usr/lib64"),  # RedHat/CentOS
                    pathlib.Path("/lib64"),  # Some distributions
                ]
            )
        elif arch.startswith("arm") or arch.startswith("aarch"):
            paths.extend(
                [
                    pathlib.Path("/usr/lib/arm-linux-gnueabi"),
                    pathlib.Path("/usr/lib/aarch64-linux-gnu"),
                ]
            )
    elif system == "darwin":  # macOS
        paths.extend(
            [
                pathlib.Path("/opt/homebrew/lib"),
                pathlib.Path("/usr/local/lib"),
                pathlib.Path("/usr/local/opt/gcc/lib"),
                pathlib.Path("/opt/homebrew/opt/gcc/lib"),
                pathlib.Path("/opt/homebrew/opt/openblas/lib"),
                pathlib.Path("/opt/homebrew/opt/lapack/lib"),
                pathlib.Path("/usr/local/opt/opt/openblas/lib"),
                pathlib.Path("/usr/local/opt/opt/lapack/lib"),
            ]
        )
    return paths


def get_search_paths(local_path: pathlib.Path, verbose: bool = False) -> list[pathlib.Path]:
    """Return all the directories where the libraries are searched.

    Parameters
    ----------
    local_path : pathlib.Path
        The package directory (searched first, together with its ``lib``
        subdirectory).
    verbose : bool, default False
        Print the directories.

    Returns
    -------
    list[pathlib.Path]
        Existing directories: the package directory, the directories in
        ``LD_LIBRARY_PATH``, ``LIBRARY_PATH``, ``DYLD_LIBRARY_PATH`` and the
        platform directories.
    """
    search_paths: set[pathlib.Path] = set()

    # Add local paths
    search_paths.add(local_path)
    search_paths.add(local_path / "lib")

    # Add paths from environment variables
    env_vars = ["LD_LIBRARY_PATH", "LIBRARY_PATH", "DYLD_LIBRARY_PATH"]
    for env_var in env_vars:
        paths = os.environ.get(env_var, "")
        if paths:
            for p in paths.split(":"):
                if p:
                    path = pathlib.Path(p)
                    search_paths.add(path)
                    search_paths.add(path / "lib")

    # Add platform-specific paths
    search_paths.update(get_platform_specific_paths())

    # Filter non-existent paths and convert to list
    valid_paths = [p for p in search_paths if p.is_dir()]

    if verbose:
        print("Searching in the following paths:")
        for path in valid_paths:
            print(f"  - {path}")

    return valid_paths


def get_lib_extensions() -> list[str]:
    """Return the shared library file extensions of the current platform.

    Returns
    -------
    list[str]
        ``[".so", ".so.*"]`` on Linux, ``[".dylib", ".dylib.*"]`` on macOS,
        ``[".dll"]`` on Windows (all of them on other systems).
    """
    system, _ = get_platform_info()
    return {"linux": [".so", ".so.*"], "darwin": [".dylib", ".dylib.*"], "windows": [".dll"]}.get(
        system, [".so", ".dylib", ".dll"]
    )


def find_versioned_library(base_path: pathlib.Path, lib_name: str) -> pathlib.Path | None:
    """Find a versioned shared library in a directory (e.g. ``libfoo.so.3``).

    Parameters
    ----------
    base_path : pathlib.Path
        Directory to search.
    lib_name : str
        Library file name without version, e.g. ``"libgfortran.so"``.

    Returns
    -------
    pathlib.Path or None
        The matching file with the highest version, or ``None``.
    """
    system, _ = get_platform_info()
    extensions = get_lib_extensions()
    for ext in extensions:
        if lib_name.endswith(ext):
            base_name = lib_name[: -len(ext)]
            patterns = [f"{base_name}{ext}*"]
            if system == "darwin":
                version_patterns = [".*.dylib"]
                patterns.extend(f"{base_name}{pat}" for pat in version_patterns)
            for pattern in patterns:
                matches = list(base_path.glob(pattern))
                if matches:
                    matches.sort(reverse=True)
                    return matches[0]
    return None


def find_library(
    lib_name: str, local_path: pathlib.Path, verbose: bool = False
) -> pathlib.Path | None:
    """Find a shared library in the package directory and the system paths.

    Parameters
    ----------
    lib_name : str
        File name, e.g. ``"liblapack.so"`` (``.so`` is replaced by ``.dylib``
        on macOS). Versioned files (``liblapack.so.3``) are also accepted.
    local_path : pathlib.Path
        The package directory, searched first.
    verbose : bool, default False
        Print every path checked.

    Returns
    -------
    pathlib.Path or None
        Path of the library, or ``None`` if it is not found.

    Examples
    --------
    >>> import pathlib
    >>> find_library("liblapack.so", pathlib.Path("/usr/lib"))  # doctest: +SKIP
    PosixPath('/usr/lib/x86_64-linux-gnu/liblapack.so.3')
    """
    if not isinstance(local_path, pathlib.Path):
        local_path = pathlib.Path(local_path)

    system, _ = get_platform_info()
    search_paths = get_search_paths(local_path, verbose)
    if system == "darwin" and lib_name.endswith(".so"):
        lib_name = f"{lib_name[:-3]}.dylib"

    for path in search_paths:
        full_path = path / lib_name
        if verbose:
            print(f"Checking {full_path}")
        if full_path.exists():
            return full_path
        versioned_lib = find_versioned_library(path, lib_name)
        if versioned_lib:
            return versioned_lib

    if verbose:
        print(f"Could not find library: {lib_name}")
    return None


def load_dependencies_for_darwin(package_lib_path: pathlib.Path):
    """Preload the BLAS/LAPACK and OpenMP libraries on macOS, if found.

    Parameters
    ----------
    package_lib_path : pathlib.Path
        The package directory.
    """
    system, _ = get_platform_info()
    if system != "darwin":
        return

    # Try to load OpenMP runtime (required for some BLAS/LAPACK implementations)
    omp_names = ["libomp.dylib", "libgomp.dylib"]
    for name in omp_names:
        omp_path = find_library(name, package_lib_path)
        if omp_path:
            try:
                _ = ctypes.CDLL(str(omp_path))
                break
            except OSError:
                continue
    # Load BLAS and LAPACK
    for lib_name in ["libblas.dylib", "liblapack.dylib"]:
        lib_path = find_library(lib_name, package_lib_path)
        if lib_path:
            try:
                _ = ctypes.CDLL(str(lib_path))
            except OSError as e:
                print(f"Warning: Failed to load {lib_name}: {e}")


def load_gfortran_for_posix(package_lib_path: pathlib.Path):
    """Preload ``libgfortran`` (versions 3 to 5) if found (Linux).

    Parameters
    ----------
    package_lib_path : pathlib.Path
        The package directory.
    """
    for i in range(3, 6):
        gfortran_path = find_library(f"libgfortran.so.{i}", package_lib_path)
        if gfortran_path:
            _ = ctypes.CDLL(str(gfortran_path))
            break


def load_quadmath_for_posix(package_lib_path: pathlib.Path):
    """Preload ``libquadmath`` if found (Linux; not available on every architecture).

    Parameters
    ----------
    package_lib_path : pathlib.Path
        The package directory.
    """
    quadmath_names = ["libquadmath.so.0", "libquadmath.so.0.0.0"]
    for name in quadmath_names:
        quadmath_path = find_library(name, package_lib_path)
        if quadmath_path:
            try:
                _ = ctypes.CDLL(str(quadmath_path))
                break
            except OSError:
                continue


def load_library_for_posix():
    """Load libscientific on Linux and macOS.

    The library bundled in the wheel is loaded first (its dependencies are
    bundled too). Otherwise the dependencies (gfortran, quadmath, BLAS, LAPACK)
    are preloaded and ``libscientific.so``/``.dylib`` is searched in the system
    paths (see :func:`get_search_paths`).

    Returns
    -------
    ctypes.CDLL
        The loaded library.

    Raises
    ------
    RuntimeError
        If the library or one of its dependencies cannot be found or loaded; the
        message explains how to install it or set ``LD_LIBRARY_PATH`` /
        ``DYLD_LIBRARY_PATH``.
    """
    try:
        package_lib_path = pathlib.Path(__file__).parent.absolute()
        system, _ = get_platform_info()

        if system == "darwin":
            load_dependencies_for_darwin(package_lib_path)
            lib_path = find_library("libscientific.dylib", package_lib_path)
            if lib_path:
                return ctypes.CDLL(str(lib_path))
            raise OSError("Could not find libscientific.dylib")

        # Handle Linux
        # Wheels built by the release pipeline (auditwheel) bundle the
        # dependencies and set RUNPATH: load the packaged library directly.
        bundled = package_lib_path / "libscientific.so"
        if bundled.exists():
            try:
                return ctypes.CDLL(str(bundled))
            except OSError:
                pass  # dependencies not resolvable: fall back to preloading

        load_gfortran_for_posix(package_lib_path)
        load_quadmath_for_posix(package_lib_path)
        # Load BLAS and LAPACK
        for lib_name in ["libblas.so", "liblapack.so"]:
            lib_path = find_library(lib_name, package_lib_path)
            if lib_path:
                _ = ctypes.CDLL(str(lib_path))
            else:
                raise OSError(f"Could not find {lib_name}")
        # Load libscientific
        scientific_path = find_library("libscientific.so", package_lib_path)
        if scientific_path:
            return ctypes.CDLL(str(scientific_path))
        raise OSError("Could not find libscientific.so")
    except (OSError, RuntimeError) as err:
        msg = [
            f"Failed to load libscientific on {system}",
            "Please either:",
            "1. Install the library from source:",
            "   https://github.com/gmrandazzo/libscientific",
            "2. If already installed, specify location in appropriate environment variable:",
            "   Linux: export LD_LIBRARY_PATH=<path to libscientific>",
            "   macOS: export DYLD_LIBRARY_PATH=<path to libscientific>",
        ]
        raise RuntimeError("\n".join(msg)) from err


def load_library_for_nt():
    """Load ``libscientific.dll`` from the package directory (Windows).

    The runtime DLLs it needs (LAPACK, BLAS, gfortran, ...) are bundled in the
    same directory.

    Returns
    -------
    ctypes.WinDLL or None
        The loaded library, or ``None`` if it cannot be loaded.
    """
    try:
        lib_path = f"{pathlib.Path(__file__).parent}"
        win_dll = getattr(ctypes, "WinDLL", ctypes.CDLL)
        return win_dll(f"{lib_path}\\libscientific.dll")
    except (TypeError, OSError):
        return None


def load_libscientific_library():
    """Load the libscientific C library for the current platform.

    Every module of the package calls this function at import time; the result
    is the ``ctypes`` handle used to call the C functions.

    Returns
    -------
    ctypes.CDLL
        The loaded library.

    Raises
    ------
    RuntimeError
        If the library cannot be loaded (POSIX systems).
    """
    if os.name == "nt":
        return load_library_for_nt()
    return load_library_for_posix()
