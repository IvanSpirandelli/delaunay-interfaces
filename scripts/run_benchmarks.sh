#!/usr/bin/env bash
# Run the performance suite and save results stamped with branch + commit:
#   benchmarks/results/<branch>_<sha>.json  - full record incl. invariants
#
# Usage: scripts/run_benchmarks.sh "one-liner describing the change"
# Env overrides: SIZES=1000,5000 REPEATS=3 JULIA_SIZES=1000,10000 SKIP_JULIA=1
set -euo pipefail
cd "$(dirname "$0")/.."

NOTE=${1:-}

SIZES=${SIZES:-1000,5000,20000,50000}
REPEATS=${REPEATS:-3}
JULIA_SIZES=${JULIA_SIZES:-1000,10000,30000}

cmake -S . -B build -DBUILD_BENCHMARKS=ON > /dev/null
cmake --build build --target benchmark_speed -j4 > /dev/null
if [ -z "${SKIP_JULIA:-}" ]; then
    cmake --build build --target delaunay_interfaces_jl -j4 > /dev/null
fi

branch=$(git rev-parse --abbrev-ref HEAD)
sha=$(git rev-parse --short HEAD)
if git diff --quiet && git diff --cached --quiet; then worktree=clean; else worktree=dirty; fi
[ "$worktree" = dirty ] && echo "WARNING: working tree is dirty; results tagged 'dirty'"

outdir=benchmarks/results
mkdir -p "$outdir"
core_json=$(mktemp)
julia_json=$(mktemp)

./build/benchmarks/benchmark_speed --sizes "$SIZES" --repeats "$REPEATS" --out "$core_json"

if [ -z "${SKIP_JULIA:-}" ]; then
    julia --project=julia scripts/benchmark_julia_boundary.jl "$julia_json" "$JULIA_SIZES" "$REPEATS"
else
    echo "[]" > "$julia_json"
fi

out="$outdir/${branch}_${sha}.json"
CORE_JSON="$core_json" JULIA_JSON="$julia_json" OUT="$out" BRANCH="$branch" SHA="$sha" \
WORKTREE="$worktree" SIZES="$SIZES" REPEATS="$REPEATS" NOTE="$NOTE" \
python3 - <<'EOF'
import datetime, json, os, platform

core = json.load(open(os.environ["CORE_JSON"]))
julia = json.load(open(os.environ["JULIA_JSON"]))

doc = {
    "branch": os.environ["BRANCH"],
    "commit": os.environ["SHA"],
    "worktree": os.environ["WORKTREE"],
    "date": datetime.datetime.now().isoformat(timespec="seconds"),
    "machine": f"{platform.system()} {platform.machine()}",
    "sizes": os.environ["SIZES"],
    "repeats": int(os.environ["REPEATS"]),
    "note": os.environ["NOTE"],
    "core": core,
    "julia_boundary": julia,
}
with open(os.environ["OUT"], "w") as f:
    json.dump(doc, f, indent=1)
print("saved", os.environ["OUT"])
EOF

rm -f "$core_json" "$julia_json"
