#!/usr/bin/env bash
# Build the extension with AddressSanitizer and UndefinedBehaviorSanitizer
# and run the test suite against it (Linux, GCC or Clang).
#
#     tools/sanitize.sh [pytest arguments]
#
# The extension is compiled by hand into a scratch directory together with a
# copy of the Python sources, so the working tree and the installed package
# stay untouched and the compiler flags are known exactly. Requires numpy,
# pytest and tifffile in the active Python environment.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$(mktemp -d)"
trap 'rm -rf "$OUT"' EXIT
CXX="${CXX:-g++}"
PYTHON="${PYTHON:-python}"

INCLUDES="$("$PYTHON" -c "import sysconfig, numpy; print('-I' + sysconfig.get_paths()['include'], '-I' + numpy.get_include())")"
SUFFIX="$("$PYTHON" -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")"

mkdir -p "$OUT/pystackreg"
cp -r "$REPO"/pystackreg/*.py "$REPO"/pystackreg/py.typed "$REPO"/pystackreg/util "$OUT/pystackreg/"

# shellcheck disable=SC2086
"$CXX" -std=c++11 -ffp-contract=off -O1 -g -fPIC -shared -fno-omit-frame-pointer \
    -fsanitize=address,undefined,float-cast-overflow \
    -fno-sanitize-recover=undefined,float-cast-overflow \
    -I"$REPO/inc" $INCLUDES \
    "$REPO"/src/pymain.cpp "$REPO"/src/TurboReg.cpp "$REPO"/src/TurboRegMask.cpp \
    "$REPO"/src/TurboRegImage.cpp "$REPO"/src/TurboRegTransform.cpp \
    "$REPO"/src/TurboRegPointHandler.cpp \
    -o "$OUT/pystackreg/turboreg$SUFFIX"

if ! nm -D "$OUT/pystackreg/turboreg$SUFFIX" | grep -q __asan; then
    echo "the extension was not instrumented" >&2
    exit 1
fi

# The sanitizer runtimes must be loaded before the (uninstrumented) Python
# executable. Leak detection is off: CPython itself does not free everything.
export LD_PRELOAD="$("$CXX" -print-file-name=libasan.so):$("$CXX" -print-file-name=libubsan.so)"
export ASAN_OPTIONS="detect_leaks=0:halt_on_error=1"
export UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1"
export PYTHONPATH="$OUT"

cd "$OUT"
"$PYTHON" -m pytest "$REPO/tests" --rootdir "$REPO" -c "$REPO/pyproject.toml" \
    -p no:cacheprovider "$@"
