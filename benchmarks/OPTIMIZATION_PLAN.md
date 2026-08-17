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
- [ ] **Iter 2 — result-assembly moves**: move-out accessors in
  `simplified_subdivision.cpp` / `barycentric_subdivision.cpp` result
  assembly (currently copies from dying locals); Python bindings take
  `InterfaceSurface&&` / `SimplifiedSurface&&` and `std::move` members;
  `std::move(filtration)` in `compute_barycentric_..._py` returns.
  Expect wins mainly in total_ms at large N; invariants must stay bitwise.
- [ ] **Iter 3 — Julia FFI boundary**: use the existing bulk
  `get_all_vertices` in the wrapper constructor (`reshape`); add bulk flat
  accessors for the filtration (follow `get_vertex_atom_indices_flat`
  pattern) and decode in one pass; drop per-element `get_vertex` /
  `get_simplex_vertices` / `get_simplex_value` loops. Measured by the
  `julia_boundary` section (wrapper_ms was ~26% of cxx_ms at 30k).
  Invariants bitwise.
- [ ] **Iter 4 — spatial-sorted range insert** for plain and weighted
  Delaunay paths in `interface_generation.cpp` (alpha paths already do
  this). NOTE: traversal order changes relabel vertex ids, so the
  filtration's tie order may permute — counts and filtration_sum stay
  identical (output is value-sorted; ties share equal values), but a raw
  diff of simplex ids against old output is not expected to match.
  Affects collect_ms; baseline says this is a modest win (collect is
  ~3-30% of total depending on scenario).
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
