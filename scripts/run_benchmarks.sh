#!/usr/bin/env bash
# Run the performance suite and save results stamped with branch + commit:
#   benchmarks/results/<branch>_<sha>.json  - full record incl. invariants
#   benchmarks/history.json                 - one entry per commit with note
#     and per-test [collect_ms, total_ms], for tracking across iterations
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

./build/tests/benchmark_speed --sizes "$SIZES" --repeats "$REPEATS" --out "$core_json"

if [ -z "${SKIP_JULIA:-}" ]; then
    julia --project=julia scripts/benchmark_julia_boundary.jl "$julia_json" "$JULIA_SIZES" "$REPEATS"
else
    echo "[]" > "$julia_json"
fi

out="$outdir/${branch}_${sha}.json"
CORE_JSON="$core_json" JULIA_JSON="$julia_json" OUT="$out" BRANCH="$branch" SHA="$sha" \
WORKTREE="$worktree" SIZES="$SIZES" REPEATS="$REPEATS" NOTE="$NOTE" \
HISTORY="benchmarks/history.json" python3 - <<'EOF'
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

# History: one entry per commit (re-running a commit replaces its entry),
# timings only - the per-run file above keeps the invariants.
entry = {
    "commit_id": os.environ["SHA"],
    "branch": os.environ["BRANCH"],
    "date": doc["date"],
    "note": os.environ["NOTE"],
    "tests": {f"{r['n']}/{r['colors']}/{r['scenario']}": [r["collect_ms"], r["total_ms"]]
              for r in core},
    "julia_boundary": {str(r["n"]): [r["cxx_ms"], r["wrapper_ms"]] for r in julia},
}
path = os.environ["HISTORY"]
history = json.load(open(path)) if os.path.exists(path) else []
history = [h for h in history if h["commit_id"] != entry["commit_id"]] + [entry]
with open(path, "w") as f:
    json.dump(history, f, indent=1)
print("appended", path)
EOF

rm -f "$core_json" "$julia_json"
