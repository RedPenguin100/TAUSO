#!/usr/bin/env bash
set -euo pipefail

# Debug mode toggle: TAUSO_RACCESS_DEBUG=1 tauso-setup-raccess
if [[ "${TAUSO_RACCESS_DEBUG:-0}" != "0" ]]; then
    echo "DEBUG: install_raccess.sh starting"
    echo "DEBUG: PWD: $(pwd)"
    echo "DEBUG: ARGS: $*"
    set -x
fi

UPSTREAM_URL="https://github.com/gterai/raccess.git"
TARGET_COMMIT="ad29d098b62606aa4051951469780b4e7be54536"

FORCE_CLONE=0

# First arg may be a flag
if [ "${1:-}" = "--force-clone" ] || [ "${1:-}" = "-f" ]; then
    FORCE_CLONE=1
    shift
fi

TARGET_DIR="${1:-$HOME/.local/share/tauso/raccess}"

echo "INFO: Initial TARGET_DIR: $TARGET_DIR"
mkdir -p "$(dirname "$TARGET_DIR")"

# Make TARGET_DIR absolute (portable)
TARGET_DIR="$(cd "$TARGET_DIR" && pwd)"
echo "INFO: Absolute TARGET_DIR: $TARGET_DIR"
echo "INFO: FORCE_CLONE: $FORCE_CLONE"
echo "INFO: UPSTREAM_URL: $UPSTREAM_URL"
echo "INFO: TARGET_COMMIT: $TARGET_COMMIT"

# Show where this script lives & what files exist there
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "INFO: SCRIPT_DIR: $SCRIPT_DIR"
echo "INFO: Contents of SCRIPT_DIR:"
ls -al "$SCRIPT_DIR" || echo "WARN: cannot list SCRIPT_DIR"

PATCH_DIR="$SCRIPT_DIR/patches"
echo "INFO: PATCH_DIR: $PATCH_DIR"
ls -al "$PATCH_DIR" || echo "WARN: cannot list PATCH_DIR"

if [ "$FORCE_CLONE" -eq 1 ]; then
    echo "INFO: FORCE_CLONE=1, removing existing TARGET_DIR (if any): $TARGET_DIR"
    rm -rf "$TARGET_DIR"
fi

if [ -d "$TARGET_DIR/.git" ]; then
    echo "INFO: Existing git repo detected at $TARGET_DIR/.git"
    local_head=$(git -C "$TARGET_DIR" rev-parse HEAD || echo "UNKNOWN")
    echo "INFO: local HEAD: $local_head"

    if [ "$local_head" = "$TARGET_COMMIT" ] && [ "$FORCE_CLONE" -eq 0 ]; then
        echo "INFO: Already at desired commit $TARGET_COMMIT — nothing to do."
        exit 0
    fi

    echo "INFO: Repo exists but not at desired commit or FORCE_CLONE requested; updating..."
    git -C "$TARGET_DIR" fetch --all
    git -C "$TARGET_DIR" checkout -f "$TARGET_COMMIT"
else
    echo "INFO: No repo at $TARGET_DIR, cloning..."
    git clone "$UPSTREAM_URL" "$TARGET_DIR"
    git -C "$TARGET_DIR" checkout -f "$TARGET_COMMIT"
fi

echo "INFO: After clone/checkout, contents of TARGET_DIR:"
ls -al "$TARGET_DIR" || echo "WARN: cannot list TARGET_DIR"
echo "INFO: Git status:"
git -C "$TARGET_DIR" status || echo "WARN: git status failed"

cd "$TARGET_DIR"

echo "INFO: Applying patches from $PATCH_DIR"
for p in raccess_bugfix.patch \
         raccess_fix.patch \
         raccess_fix_link.patch \
         raccess_makefile.patch
do
    echo "INFO: Applying patch: $p"
    if [ -f "$PATCH_DIR/$p" ]; then
        patch -p1 < "$PATCH_DIR/$p"
    else
        echo "ERROR: Patch file not found: $PATCH_DIR/$p"
        exit 1
    fi
done

echo "INFO: Contents of TARGET_DIR/src before build:"
ls -al "$TARGET_DIR/src" || echo "WARN: cannot list $TARGET_DIR/src"

cd "$TARGET_DIR/src"

echo "INFO: Running make in $(pwd)"
make

echo
echo "INFO: raccess built at: $TARGET_DIR/src/raccess/run_raccess"
echo "INFO: To use raccess, add this to your shell (e.g. ~/.bashrc):"
echo "  export RACCESS_EXE=\"$TARGET_DIR/src/raccess/run_raccess\""
echo
echo "INFO: install_raccess.sh finis_
