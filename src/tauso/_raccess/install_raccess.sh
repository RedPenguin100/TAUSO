#!/usr/bin/env bash
set -euo pipefail

# --- Check for GCC availability ---
if ! command -v gcc >/dev/null 2>&1; then
    echo "ERROR: GCC compiler not found."
    echo "You must install GCC before running tauso install-raccess."
    echo ""
    echo "On Debian/Ubuntu:"
    echo "    sudo apt-get install build-essential"
    echo ""
    echo "On Fedora:"
    echo "    sudo dnf install gcc gcc-c++ make"
    echo ""
    echo "On macOS (Homebrew):"
    echo "    brew install gcc"
    echo ""
    echo "Or install a conda compiler toolchain:"
    echo "    conda install -c conda-forge gcc"
    exit 1
fi

UPSTREAM_URL="https://github.com/gterai/raccess.git"
TARGET_COMMIT="ad29d098b62606aa4051951469780b4e7be54536"

FORCE_CLONE=0

if [ "${1:-}" = "--force-clone" ] || [ "${1:-}" = "-f" ]; then
    FORCE_CLONE=1
    shift
fi

DEFAULT_DATA_DIR="${TAUSO_DATA_DIR:-$HOME/.local/share/tauso}"
DATA_DIR="${1:-$DEFAULT_DATA_DIR}"
DATA_DIR="$(readlink -f "$DATA_DIR")"

BUILD_DIR="$DATA_DIR/src/raccess"

if [ "$FORCE_CLONE" = "1" ] && [ -d "$BUILD_DIR" ]; then
    echo "Force flag detected. Removing existing clone..."
    rm -rf "$BUILD_DIR"
fi

if [ -d "$BUILD_DIR/.git" ]; then
    local_head=$(git -C "$BUILD_DIR" rev-parse HEAD)
    if [ "$local_head" = "$TARGET_COMMIT" ]; then
        if [ -x "$BUILD_DIR/src/raccess/run_raccess" ]; then
            echo "Already at desired commit $TARGET_COMMIT and binary exists — nothing to do."
            exit 0
        else
            echo "Already at desired commit, but binary missing. Proceeding to compilation..."
        fi
    else
        echo "Existing repo but wrong commit, checking out desired commit..."
        git -C "$BUILD_DIR" fetch --all
        git -C "$BUILD_DIR" checkout -f "$TARGET_COMMIT"
    fi
else
    echo "Cloning raccess..."
    git clone "$UPSTREAM_URL" "$BUILD_DIR"
    git -C "$BUILD_DIR" checkout -f "$TARGET_COMMIT"
fi

cd "$BUILD_DIR"

PATCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/patches"

for p in raccess_bugfix.patch \
         raccess_fix.patch \
         raccess_fix_link.patch \
         raccess_makefile.patch
do
    echo "Applying patch: $p"
    patch -N -p1 < "$PATCH_DIR/$p" || true
done

cd "src"

echo "Building raccess with matrix-math optimizations..."

export CXXFLAGS="-O3 -march=native -ffast-math -fno-rtti"
export CFLAGS="-O3 -march=native -ffast-math"
export LDFLAGS=""

make clean
make -j$(nproc)

echo
echo "raccess built at: $BUILD_DIR/src/raccess/run_raccess"

mkdir -p "$DATA_DIR/bin"
cp "$BUILD_DIR/src/raccess/run_raccess" "$DATA_DIR/bin/run_raccess"

echo "Successfully installed executable to: $DATA_DIR/bin/run_raccess"

PACKAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$PACKAGE_DIR/bin"
cp "$BUILD_DIR/src/raccess/run_raccess" "$PACKAGE_DIR/bin/"

echo "Successfully copied executable to: $PACKAGE_DIR/bin/run_raccess"
echo "It will now be auto-discovered by tauso."
echo
echo "Done."