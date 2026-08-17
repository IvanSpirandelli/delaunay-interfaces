# Performance optimization plan

Working through the optimization findings from the 2026-08-17 review, one
iteration per commit, tracked in `history.json`. Run via:

```bash
scripts/run_benchmarks.sh "one-liner describing the change"
```

## Protocol (every iteration)

1. Implement one optimization batch.
2. `cmake --build build -j4` and `ctest` — all 5 suites must pass.
3. Invariant check against the baseline per-run file
   (`results/restructure-examples_f2babcf.json`): counts and
   `filtration_sum` must match. Bitwise-identical unless the iteration's
   notes below say otherwise; document any expected deviation in the
   commit message.
4. Commit + push (one commit per iteration, message explains the change).
5. `scripts/run_benchmarks.sh "<note>"` on the clean tree; the result file
   and the `history.json` entry are committed afterwards with the timing
   comparison in the commit message.
6. If invariants break unexpectedly or timings regress: revert the
   iteration, record why in this file, move on.

## Queue

- [x] **Iter 1 — subdivision phase** (committed): `filtration_set_` →
  vector + sort/unique; `get_or_create_vertex` reuses `faces[i].atoms` as
  key; part barycenters computed once per new vertex (value + position),
  `face_barycenters` removed. Invariants bitwise identical (verified).
- [x] **Iter 2 — result-assembly moves** (committed, f2f34af): take_*
  accessors on both subdivision classes; Python binding wrappers take the
  surface by rvalue and move members. Invariants bitwise identical.
  1.3-1.6x on copy-heavy pipelines; cumulative 2.0-2.4x on 50k
  full-Delaunay scenarios, Julia-side cxx call 15.1s -> 6.3s at 30k.
- [x] **Iter 3 — Julia FFI boundary** (committed, 08e3573): wrapper
  constructor uses bulk get_all_vertices plus new flat filtration
  accessors (get_all_simplex_vertices_flat / get_all_simplex_values),
  decoded in one pass; per-element accessors remain bound. Wrapper output
  verified element-identical. wrapper_ms 7.9x faster at 1k, ~2x at
  10k-30k (remaining cost is per-simplex Julia allocations, inherent to
  the public Vector-of-Vectors field types).
- [x] **Iter 4 — spatial-sorted range insert** (committed, 76c8454) for
  plain and weighted Delaunay paths in `interface_generation.cpp` (alpha
  paths already do this). Counts identical to baseline; filtration_sum
  identical except one last-ulp deviation (20000/two/
  delaunay_weighted_half), the expected summation-order effect of
  permuted equal-value ties. collect_ms 1.5-2.5x on delaunay scenarios
  (50k distinct: 319 -> 130 ms); untouched alpha scenarios within the
  run's ±15% noise band. Benchmarked with two unrelated processes
  pegging 2 of 12 cores (stuck Julia LSP indexer + busy IJulia kernel);
  a first, discarded run with 3 pegged cores was uniformly ~2x slow.
- [ ] **Iter 5 (optional) — vertex_map_ key packing**: `std::map` with
  `vector<int>` keys → packed fixed-size key (≤ 4 atoms) in a hash map.
  Vertex ids are assigned in creation order (map order unused), so output
  is unchanged. Only worth it if iter 1 left the map visible in profiles.

Deliberately skipped: pybind `def_readonly` re-conversion caching
(semantic change to attribute access), visualization-side performance
(cold path).

## Baseline reference (f2babcf, Darwin arm64, median of 3)

Worst-case pipelines: 50k/two delaunay total 13.9s (collect 0.4s — the
subdivision dominates); 50k/distinct alpha_uniform_half collect 2.1s.
Julia boundary at 30k: cxx 15.1s + wrapper 3.9s.
