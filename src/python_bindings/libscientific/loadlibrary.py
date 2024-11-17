"""loadlibrary.py method to load libscientific library

Copyright (C) <2023>  Giuseppe Marco Randazzo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
import pathlib
import ctypes
from typing import (
    Dict,
    Optional,
    Tuple,
    List,
    Set
)

import platform

def get_platform_info() -> Tuple[str, str]:
    """Get detailed platform and architecture information.
    
    Returns:
        Tuple of (system_name, architecture)
    """
    system = platform.system().lower()
    machine = platform.machine().lower()
    arch_map = {
        'x86_64': 'x86_64',
        'amd64': 'x86_64',
        'arm64': 'arm64',
        'aarch64': 'arm64'
    }
    arch = arch_map.get(machine, machine)
    return system, arch

def get_platform_specific_paths() -> List[pathlib.Path]:
    """Get platform-specific library paths based on the OS and architecture.
    
    Returns:
        List of platform-specific library paths
    """
    system, arch = get_platform_info()
    paths = [
        pathlib.Path('/usr/lib'),
        pathlib.Path('/usr/local/lib'),
    ]
    if system == 'linux':
        # Get machine architecture
        if arch == 'x86_64':
            paths.extend([
                pathlib.Path('/usr/lib/x86_64-linux-gnu'),  # Debian/Ubuntu
                pathlib.Path('/usr/lib64'),                 # RedHat/CentOS
                pathlib.Path('/lib64'),                     # Some distributions
            ])
        elif arch.startswith('arm') or arch.startswith('aarch'):
            paths.extend([
                pathlib.Path('/usr/lib/arm-linux-gnueabi'),
                pathlib.Path('/usr/lib/aarch64-linux-gnu'),
            ])
    elif system == 'darwin':  # macOS
        paths.extend([
            pathlib.Path('/opt/homebrew/lib'),
            pathlib.Path('/usr/local/lib'),
            pathlib.Path('/usr/local/opt/gcc/lib'),
            pathlib.Path('/opt/homebrew/opt/gcc/lib'),
            pathlib.Path('/opt/homebrew/opt/openblas/lib'),
            pathlib.Path('/opt/homebrew/opt/lapack/lib'),
            pathlib.Path('/usr/local/opt/opt/openblas/lib'),
            pathlib.Path('/usr/local/opt/opt/lapack/lib'),
        ])
    return paths

def get_search_paths(local_path: pathlib.Path, verbose: bool = False) -> List[pathlib.Path]:
    """Get all possible library search paths including environment variables.
    
    Args:
        local_path: Local directory to include in search
        verbose: Whether to print debug information
    
    Returns:
        List of unique, existing library paths
    """
    search_paths: Set[pathlib.Path] = set()

    # Add local paths
    search_paths.add(local_path)
    search_paths.add(local_path / 'lib')

    # Add paths from environment variables
    env_vars = ['LD_LIBRARY_PATH', 'LIBRARY_PATH', 'DYLD_LIBRARY_PATH']
    for env_var in env_vars:
        paths = os.environ.get(env_var, '')
        if paths:
            for p in paths.split(':'):
                if p:
                    path = pathlib.Path(p)
                    search_paths.add(path)
                    search_paths.add(path/'lib')

    # Add platform-specific paths
    search_paths.update(get_platform_specific_paths())

    # Filter non-existent paths and convert to list
    valid_paths = [p for p in search_paths if p.is_dir()]

    if verbose:
        print("Searching in the following paths:")
        for path in valid_paths:
            print(f"  - {path}")

    return valid_paths

def get_lib_extensions() -> Dict[str, List[str]]:
    """Get library extensions and patterns for different platforms.
    
    Returns:
        Dictionary of platform-specific library patterns
    """
    system, _ = get_platform_info()
    return {
        'linux': ['.so', '.so.*'],
        'darwin': ['.dylib', '.dylib.*'],
        'windows': ['.dll']
    }.get(system, ['.so', '.dylib', '.dll'])

def find_versioned_library(base_path: pathlib.Path, lib_name: str) -> Optional[pathlib.Path]:
    """Search for versioned library files.
    
    Args:
        base_path: Directory to search in
        lib_name: Library name to search for
    
    Returns:
        Path to the library if found, None otherwise
    """
    system, _ = get_platform_info()
    extensions = get_lib_extensions()
    for ext in extensions:
        if lib_name.endswith(ext):
            base_name = lib_name[:-len(ext)]
            patterns = [f'{base_name}{ext}*']
            if system == 'darwin':
                version_patterns = ['.*.dylib']
                patterns.extend(f'{base_name}{pat}' for pat in version_patterns)
            for pattern in patterns:
                matches = list(base_path.glob(pattern))
                if matches:
                    matches.sort(reverse=True)
                    return matches[0]
    return None

def find_library(
        lib_name: str,
        local_path: pathlib.Path,
        verbose: bool = False
    ) -> Optional[pathlib.Path]:
    """Try to find a library in system paths and local directory.
    
    Args:
        lib_name: Name of the library to find
        local_path: Local directory to check
        verbose: Whether to print debug information
    
    Returns:
        Path to library if found, None otherwise
    
    Example:
        >>> find_library('liblapack.so', Path('/usr/local'), verbose=True)
    """
    if not isinstance(local_path, pathlib.Path):
        local_path = pathlib.Path(local_path)

    system, _ = get_platform_info()
    search_paths = get_search_paths(local_path, verbose)
    if system == 'darwin' and lib_name.endswith('.so'):
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
    """Load required dependencies for macOS (both ARM and x86_64).
    
    Args:
        package_lib_path: Path to the package library directory
    """
    system, _ = get_platform_info()
    if system != 'darwin':
        return

    # Try to load OpenMP runtime (required for some BLAS/LAPACK implementations)
    omp_names = ['libomp.dylib', 'libgomp.dylib']
    for name in omp_names:
        omp_path = find_library(name, package_lib_path)
        if omp_path:
            try:
                _ = ctypes.CDLL(str(omp_path))
                break
            except OSError:
                continue
    # Load BLAS and LAPACK
    for lib_name in ['libblas.dylib', 'liblapack.dylib']:
        lib_path = find_library(lib_name, package_lib_path)
        if lib_path:
            try:
                _ = ctypes.CDLL(str(lib_path))
            except OSError as e:
                print(f"Warning: Failed to load {lib_name}: {e}")

def load_gfortran_for_posix(package_lib_path: pathlib.Path):
    """
    Load libgfortran
    """
    for i in range(3, 6):
        gfortran_path = find_library(f'libgfortran.so.{i}', package_lib_path)
        if gfortran_path:
            _ = ctypes.CDLL(str(gfortran_path))
            break

def load_quadmath_for_posix(package_lib_path: pathlib.Path):
    """
    Load libquadmath (optional for non-ARM systems)
    """
    quadmath_names = ['libquadmath.so.0', 'libquadmath.so.0.0.0']
    for name in quadmath_names:
        quadmath_path = find_library(name, package_lib_path)
        if quadmath_path:
            try:
                _ = ctypes.CDLL(str(quadmath_path))
                break
            except OSError:
                continue

def load_library_for_posix():
    """Load the libscientific library and its dependencies for POSIX systems.
    
    Returns:
        Loaded libscientific library
    
    Raises:
        RuntimeError: If libraries cannot be loaded
    """
    try:
        package_lib_path = pathlib.Path(__file__).parent.absolute()
        system, _ = get_platform_info()

        if system == 'darwin':
            load_dependencies_for_darwin(package_lib_path)
            lib_path = find_library('libscientific.dylib', package_lib_path)
            if lib_path:
                return ctypes.CDLL(str(lib_path))
            raise OSError("Could not find libscientific.dylib")

        # Handle Linux
        load_gfortran_for_posix(package_lib_path)
        load_quadmath_for_posix(package_lib_path)
        # Load BLAS and LAPACK
        for lib_name in ['libblas.so', 'liblapack.so']:
            lib_path = find_library(lib_name, package_lib_path)
            if lib_path:
                _ = ctypes.CDLL(str(lib_path))
            else:
                raise OSError(f"Could not find {lib_name}")
        # Load libscientific
        scientific_path = find_library('libscientific.so', package_lib_path)
        if scientific_path:
            return ctypes.CDLL(str(scientific_path))
        raise OSError("Could not find libscientific.so")
    except (OSError,RuntimeError) as err:
        msg = [
            f"Failed to load libscientific on {system}",
            "Please either:",
            "1. Install the library from source:",
            "   https://github.com/gmrandazzo/libscientific",
            "2. If already installed, specify location in appropriate environment variable:",
            "   Linux: export LD_LIBRARY_PATH=<path to libscientific>",
            "   macOS: export DYLD_LIBRARY_PATH=<path to libscientific>"
        ]
        raise RuntimeError('\n'.join(msg)) from err

def load_library_for_nt():
    """Load the libscientific library for NT systems
    """
    try:
        lib_path = f'{pathlib.Path(__file__).parent}'
        return ctypes.WinDLL(f'{lib_path}\\libscientific.dll')
    except TypeError:
        return None

def load_libscientific_library():
    """Load the libscientific library dynamically.

    Returns
    -------
    CDLL
        A ctypes dynamic link library object representing the loaded libscientific library.

    Raises
    ------
    RuntimeError
        If the library cannot be loaded or the platform is not supported.
    """
    if os.name == "nt":
        return load_library_for_nt()
    return load_library_for_posix()
