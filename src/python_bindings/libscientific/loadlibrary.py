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
    Optional,
    List,
    Set
)

import platform

def load_library_for_nt():
    """Load the libscientific library for NT systems
    """
    try:
        lib_path = f'{pathlib.Path(__file__).parent}'
        return ctypes.WinDLL(f'{lib_path}\\libscientific.dll')
    except TypeError:
        return None

def get_platform_specific_paths() -> List[pathlib.Path]:
    """Get platform-specific library paths based on the OS and architecture.
    
    Returns:
        List of platform-specific library paths
    """
    paths = [
        pathlib.Path('/usr/lib'),
        pathlib.Path('/usr/local/lib'),
    ]
    
    # Add platform-specific paths
    if platform.system() == 'Linux':
        # Get machine architecture
        machine = platform.machine()
        if machine == 'x86_64':
            paths.extend([
                pathlib.Path('/usr/lib/x86_64-linux-gnu'),  # Debian/Ubuntu
                pathlib.Path('/usr/lib64'),                 # RedHat/CentOS
                pathlib.Path('/lib64'),                     # Some distributions
            ])
        elif machine.startswith('arm') or machine.startswith('aarch'):
            paths.extend([
                pathlib.Path('/usr/lib/arm-linux-gnueabi'),
                pathlib.Path('/usr/lib/aarch64-linux-gnu'),
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
    for env_var in ['LD_LIBRARY_PATH', 'LIBRARY_PATH']:
        paths = os.environ.get(env_var, '')
        if paths:
            for p in paths.split(':'):
                if p:
                    path = pathlib.Path(p)
                    search_paths.add(path)
                    search_paths.add(path / 'lib')
                    
    # Add platform-specific paths
    search_paths.update(get_platform_specific_paths())
    
    # Filter non-existent paths and convert to list
    valid_paths = [p for p in search_paths if p.is_dir()]
    
    if verbose:
        print("Searching in the following paths:")
        for path in valid_paths:
            print(f"  - {path}")
            
    return valid_paths

def find_versioned_library(base_path: pathlib.Path, lib_name: str) -> Optional[pathlib.Path]:
    """Search for versioned library files.
    
    Args:
        base_path: Directory to search in
        lib_name: Library name to search for
    
    Returns:
        Path to the library if found, None otherwise
    """
    if lib_name.endswith('.so'):
        base_name = lib_name[:-3]  # Remove .so
        # Try common version patterns
        patterns = [
            f'{base_name}.so.*',    # Standard version pattern
            f'{lib_name}.*',        # Alternative version pattern
        ]
        for pattern in patterns:
            matches = list(base_path.glob(pattern))
            if matches:
                # Sort by version number to get the newest version
                matches.sort(reverse=True)
                return matches[0]
    return None

def find_library(lib_name: str, local_path: pathlib.Path, verbose: bool = False) -> Optional[pathlib.Path]:
    """Try to find a library in system paths and local directory.
    
    Args:
        lib_name: Name of the library to find
        local_path: Local directory to check
        verbose: Whether to print debug information
    
    Returns:
        Path to library if found, None otherwise
    
    Example:
        >>> find_library('libblas.so', Path('/usr/local'), verbose=True)
    """
    if not isinstance(local_path, pathlib.Path):
        local_path = pathlib.Path(local_path)
        
    search_paths = get_search_paths(local_path, verbose)
    
    for path in search_paths:
        # Try exact match
        full_path = path / lib_name
        if verbose:
            print(f"Checking {full_path}")
            
        if full_path.exists():
            return full_path
            
        # Try versioned library
        versioned_lib = find_versioned_library(path, lib_name)
        if versioned_lib:
            return versioned_lib
            
    if verbose:
        print(f"Could not find library: {lib_name}")
    return None

def load_library_for_posix():
    """Load the libscientific library and its dependencies for POSIX systems.
    
    Returns:
        Loaded libscientific library
    
    Raises:
        RuntimeError: If libraries cannot be loaded
    """
    try:
        package_lib_path = pathlib.Path(__file__).parent.absolute()
        # Handle MacOS
        if os.uname()[0] == "Darwin":
            return ctypes.CDLL(str(package_lib_path / 'libscientific.dylib'))
        else:
            # Handle Linux
            # Load libgfortran
            for i in range(3, 6):
                gfortran_path = find_library(f'libgfortran.so.{i}', package_lib_path)
                if gfortran_path:
                    _ = ctypes.CDLL(str(gfortran_path))
                    break
            
            # Load libquadmath (optional for non-ARM systems)
            quadmath_names = ['libquadmath.so.0', 'libquadmath.so.0.0.0']
            for name in quadmath_names:
                quadmath_path = find_library(name, package_lib_path)
                if quadmath_path:
                    try:
                        _ = ctypes.CDLL(str(quadmath_path))
                        break
                    except OSError:
                        continue
            
            # Load BLAS and LAPACK
            blas_path = find_library('libblas.so', package_lib_path)
            lapack_path = find_library('liblapack.so', package_lib_path)
            
            if blas_path:
                _ = ctypes.CDLL(str(blas_path))
            else:
                raise OSError("Could not find BLAS library")
                
            if lapack_path:
                _ = ctypes.CDLL(str(lapack_path))
            else:
                raise OSError("Could not find LAPACK library")
            
            # Finally load libscientific
            scientific_path = find_library('libscientific.so', package_lib_path)
    
            if scientific_path:
                return ctypes.CDLL(str(scientific_path))
            else:
                raise OSError("Could not find libscientific.so")
    except OSError as err:
        msg = [
            f"Failed to load libscientific on {os.name}",
            "Please either:",
            "1. Install the library from source:",
            "   https://github.com/gmrandazzo/libscientific",
            "2. If already installed, specify location in LD_LIBRARY_PATH:",
            "   export LD_LIBRARY_PATH=<path to libscientific>"
        ]
        raise RuntimeError('\n'.join(msg)) from err

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
