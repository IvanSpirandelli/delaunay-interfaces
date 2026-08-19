#!/usr/bin/env bash
# Profile a benchmark under a sampling profiler.
#
# Usage: scripts/profile.sh [--xctrace] [--no-open] [target] [-- <target args>]
#   target        benchmark executable name (default: benchmark_speed)
#   target args   passed through to it (default: --sizes 50000 --repeats 1)
#
# Records with samply (interactive flame graph in the Firefox Profiler UI;
# brew install samply) when installed, else falls back to Instruments' Time
# Profiler via xctrace; --xctrace forces the latter. Builds into
# build-profile/ with Release optimization plus debug info (-O3 -g) so
# inlined frames symbolicate; the normal build/ tree is untouched.
# Profiles land in benchmarks/profiles/ (gitignored).
set -euo pipefail
cd "$(dirname "$0")/.."

USE_XCTRACE=0
OPEN=1
TARGET=benchmark_speed
ARGS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --xctrace) USE_XCTRACE=1; shift ;;
        --no-open) OPEN=0; shift ;;
        --) shift; ARGS=("$@"); break ;;
        *) TARGET=$1; shift ;;
    esac
done
[ ${#ARGS[@]} -eq 0 ] && ARGS=(--sizes 50000 --repeats 1)

cmake -S . -B build-profile \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS_RELEASE="-O3 -g -DNDEBUG" \
    -DBUILD_BENCHMARKS=ON -DBUILD_TESTS=OFF -DBUILD_EXAMPLES=OFF \
    -DBUILD_PYTHON_BINDINGS=OFF -DBUILD_JULIA_BINDINGS=OFF > /dev/null
cmake --build build-profile --target "$TARGET" -j4 > /dev/null

BIN="build-profile/benchmarks/$TARGET"
STAMP=$(date +%Y%m%d-%H%M%S)
OUTDIR=benchmarks/profiles
mkdir -p "$OUTDIR"

if [ "$USE_XCTRACE" -eq 0 ] && command -v samply > /dev/null; then
    OUT="$OUTDIR/${TARGET}_${STAMP}.json.gz"
    if [ "$OPEN" -eq 0 ]; then
        samply record --save-only --output "$OUT" -- "$BIN" "${ARGS[@]}"
    else
        samply record --output "$OUT" -- "$BIN" "${ARGS[@]}"
    fi
    echo "profile saved: $OUT   (open later with: samply load $OUT)"
else
    OUT="$OUTDIR/${TARGET}_${STAMP}.trace"
    xcrun xctrace record --template 'Time Profiler' --output "$OUT" --launch -- "$BIN" "${ARGS[@]}"
    echo "profile saved: $OUT"
    [ "$OPEN" -eq 1 ] && open "$OUT"
fi
