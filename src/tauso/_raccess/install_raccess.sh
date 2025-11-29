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

# Now $1 is either the target dir OR empty
TARGET_DIR="${1:-$HOME/.local/share/tauso/raccess}"
mkdir -p "$(dirname "$TARGET_DIR")"
TARGET_DIR="$(readlink -f "$TARGET_DIR")"

echo "Installing raccess into: $TARGET_DIR"

if [ "$FORCE_CLONE" -eq 1 ]; then
    echo "Force re-clone requested, removing existing directory (if any)..."
    rm -rf "$TARGET_DIR"
fi

if [ -d "$TARGET_DIR/.git" ]; then
    # already cloned — check commit
    local_head=$(git -C "$TARGET_DIR" rev-parse HEAD)

    if [ "$local_head" = "$TARGET_COMMIT" ]; then
        echo "Already at desired commit $TARGET_COMMIT — nothing to do."
        exit 0
    fi

    echo "Existing repo but wrong commit, checking out desired commit..."
    git -C "$TARGET_DIR" fetch --all
    git -C "$TARGET_DIR" checkout -f "$TARGET_COMMIT"
else
    echo "Cloning raccess..."
    git clone "$UPSTREAM_URL" "$TARGET_DIR"
    git -C "$TARGET_DIR" checkout -f "$TARGET_COMMIT"
fi

cd "$TARGET_DIR"

PATCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/patches"

for p in raccess_bugfix.patch \
         raccess_fix.patch \
         raccess_fix_link.patch \
         raccess_makefile.patch
do
    echo "Applying patch: $p"
    patch -p1 < "$PATCH_DIR/$p"
done

cd "src"

echo "Building raccess..."
make

echo
echo "raccess built at: $TARGET_DIR/src/raccess/run_raccess"
echo "Add this to your shell (e.g. ~/.bashrc):"
echo "  export RACCESS_EXE=\"$TARGET_DIR/src/raccess/run_raccess\""
echo
echo "Done."
