#!/usr/bin/env bash
# Regenerate every figure in this folder -> figures/*.png.
# Run from this directory in the tauso_research env.
set -e
for s in *.py; do
  [ "$s" = "_feat.py" ] && continue
  echo "== $s =="
  python "$s"
done
