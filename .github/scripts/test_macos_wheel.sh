#!/usr/bin/env bash
# Build, check, install and test the libscientific macOS wheel locally,
# with the same steps as the macOS job in .github/workflows/release.yml.
#
# Usage:
#   ./test_macos_wheel.sh [path/to/libscientific]   (default: current git repo)
#
# Environment:
#   PYTHON=python3.12   python interpreter to use (>= 3.10, default python3)
#   KEEP=1              keep the built wheel and the logs (path printed at the end)
#
# Everything is built in a temporary directory; the only files written in the
# repository (the copied libscientific.dylib, python build/ and *.egg-info
# directories created by the build) are removed on exit, also on failure.
set -euo pipefail

if [ "$(uname -s)" != "Darwin" ]; then
  echo "error: this script must run on macOS" >&2
  exit 1
fi

REPO="${1:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"
REPO="$(cd "$REPO" && pwd)"
PKG_SRC="$REPO/src/python_bindings"
PKG_DIR="$PKG_SRC/libscientific"
if [ ! -f "$REPO/CMakeLists.txt" ] || [ ! -f "$PKG_SRC/setup.py" ]; then
  echo "error: $REPO is not a libscientific checkout" >&2
  exit 1
fi

PYTHON="${PYTHON:-python3}"
command -v cmake >/dev/null || { echo "error: cmake not found (brew install cmake)" >&2; exit 1; }
command -v "$PYTHON" >/dev/null || { echo "error: $PYTHON not found" >&2; exit 1; }
"$PYTHON" -c 'import sys; sys.exit(sys.version_info < (3, 10))' \
  || { echo "error: $PYTHON is older than 3.10 (set PYTHON=python3.12)" >&2; exit 1; }

# Same architecture / minimum macOS as the CI job
case "$(uname -m)" in
  arm64)  ARCH=arm64;  export MACOSX_DEPLOYMENT_TARGET=11.0 ;;
  x86_64) ARCH=x86_64; export MACOSX_DEPLOYMENT_TARGET=10.15 ;;
  *) echo "error: unsupported architecture $(uname -m)" >&2; exit 1 ;;
esac

WORK="$(mktemp -d -t libscientific-wheel)"
BUILD_DIR="$WORK/build"
DIST_DIR="$WORK/dist"
WHEELHOUSE="$WORK/wheelhouse"
VENV="$WORK/venv"

# Remember what already exists, so cleanup removes only what we create
HAD_PY_BUILD=0;  [ -e "$PKG_SRC/build" ] && HAD_PY_BUILD=1
HAD_EGG_INFO=0;  [ -e "$PKG_SRC/libscientific.egg-info" ] && HAD_EGG_INFO=1
if [ -e "$PKG_DIR/libscientific.dylib" ]; then
  echo "error: $PKG_DIR/libscientific.dylib already exists, remove it first" >&2
  exit 1
fi

cleanup() {
  status=$?
  rm -f "$PKG_DIR/libscientific.dylib"
  [ "$HAD_PY_BUILD" = 0 ] && rm -rf "$PKG_SRC/build"
  [ "$HAD_EGG_INFO" = 0 ] && rm -rf "$PKG_SRC/libscientific.egg-info"
  if [ "${KEEP:-0}" = 1 ]; then
    echo "kept: $WORK"
  else
    rm -rf "$WORK"
  fi
  if [ "$status" = 0 ]; then
    echo "== ALL OK"
  else
    echo "== FAILED (exit $status)" >&2
  fi
}
trap cleanup EXIT

step() { printf '\n== %s\n' "$*"; }

step "Build libscientific ($ARCH, macOS >= $MACOSX_DEPLOYMENT_TARGET, Release)"
cmake -S "$REPO" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_OSX_ARCHITECTURES="$ARCH"
cmake --build "$BUILD_DIR" --target libscientific -j"$(sysctl -n hw.ncpu)"
cp -L "$BUILD_DIR/src/libscientific.dylib" "$PKG_DIR/libscientific.dylib"
otool -L "$PKG_DIR/libscientific.dylib"

step "Create the virtual environment"
"$PYTHON" -m venv "$VENV"
# shellcheck disable=SC1091
source "$VENV/bin/activate"
python -m pip install -q --upgrade pip
python -m pip install -q build delocate pytest numpy

step "Build the wheel"
python -m build --wheel --outdir "$DIST_DIR" "$PKG_SRC"

step "Check the wheel with delocate"
delocate-wheel -v --require-archs "$ARCH" -w "$WHEELHOUSE" "$DIST_DIR"/*.whl
delocate-listdeps --all "$WHEELHOUSE"/*.whl
WHEEL="$(ls "$WHEELHOUSE"/*.whl)"
echo "wheel: $(basename "$WHEEL")"
case "$(basename "$WHEEL")" in
  *-py3-none-macosx_*_"$ARCH".whl) ;;
  *) echo "error: unexpected wheel tag (expected py3-none-macosx_*_$ARCH)" >&2; exit 1 ;;
esac
if delocate-listdeps --all "$WHEELHOUSE"/*.whl | grep -E '^\s*/(opt/homebrew|usr/local)/' ; then
  echo "note: the library links Homebrew libraries above (delocate bundled them into the wheel)"
fi
unzip -l "$WHEEL" | grep -E '\.dylib$'

step "Install the wheel"
python -m pip install -q --force-reinstall "$WHEEL"

step "Run the test suite against the installed wheel"
cp -r "$PKG_SRC/tests" "$WORK/tests"
cd "$WORK"
python -c "import libscientific; print(libscientific.__file__, libscientific.__version__)"
python -m pytest -q tests
