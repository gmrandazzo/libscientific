#!/usr/bin/env bash
# Build a manylinux wheel of the libscientific python bindings.
# Runs inside a quay.io/pypa/manylinux_2_28_<arch> container with the
# repository mounted on /io; the repaired wheel is written to /io/wheelhouse.
#
#   docker run --rm -v "$PWD":/io quay.io/pypa/manylinux_2_28_x86_64 \
#       /io/.github/scripts/build_linux_wheel.sh
set -euxo pipefail

PY=/opt/python/cp312-cp312/bin/python
BUILD_DIR=/tmp/libscientific-build
DIST_DIR=/tmp/libscientific-dist

dnf install -y lapack-devel blas-devel sqlite-devel

"$PY" -m pip install --upgrade pip build cmake

# Build only the shared library (tests are not needed for packaging).
"$(dirname "$PY")/cmake" -S /io -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local
"$(dirname "$PY")/cmake" --build "$BUILD_DIR" --target libscientific -j"$(nproc)"

PKG=/io/src/python_bindings/libscientific
rm -f "$PKG"/*.so*
cp -L "$BUILD_DIR/src/libscientific.so" "$PKG/libscientific.so"

rm -rf "$DIST_DIR" /io/src/python_bindings/build
"$PY" -m build --wheel --outdir "$DIST_DIR" /io/src/python_bindings

# auditwheel grafts liblapack, libblas, libgfortran, libquadmath, libsqlite3
# into the wheel, sets RUNPATH and retags it as manylinux_2_28_<arch>.
auditwheel repair "$DIST_DIR"/*.whl -w /io/wheelhouse
auditwheel show /io/wheelhouse/*.whl

rm -f "$PKG/libscientific.so"
