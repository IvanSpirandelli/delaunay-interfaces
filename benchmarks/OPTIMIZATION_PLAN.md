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
- [x] **Iter 5 — vertex_map_ key packing** (committed, 3ac5c56):
  `std::map` with `vector<int>` keys → packed `std::array<int,4>` key
  (pad -1) in an `unordered_map` with a splitmix64-mixed hash;
  process_simplex now throws on > 4 vertices, making the documented
  edge/triangle/tetrahedron contract explicit. Output bitwise identical
  (verified at 1k/5k pre-commit and in the full run). Timing: raw run
  was noisy (machine load drifted mid-suite); normalizing subdivision
  time by the collect_ms drift of the same scenario (collect is
  untouched by this change) gives ~1.0-1.15x on the subdivision phase.
  Small real win, kept.

Queue below extended 2026-08-18 from a code-reading investigation
(no profiler; estimates anchored on the per-run invariant counts:
50k/two/delaunay pushes ~8.3M raw filtration entries of which ~2.3M
(27.6%) are per-tet vertex-singleton duplicates deduping to 989k, and
each entry is a separately heap-allocated vector; ~25-68M transient
allocations per run in the per-tet enumeration).

- [ ] **Iter 6 — singleton emission + alpha predicate order + reserves**
  (S, expect ~15% on 50k delaunay totals, bitwise identical):
  (a) stop pushing per-tet vertex singletons
  (barycentric_subdivision.cpp process_simplex); emit one `({id}, value)`
  per `vertex_map_` entry in finalize_filtration instead — the sets are
  provably equal, and the sort's strict total order makes output
  depend only on the set. Guard so repeated get_filtration calls don't
  re-append. (b) interface_generation.cpp free-facet/free-edge loops:
  test cheap `is_multicolored` before the expensive `classify` —
  classify(Edge) in GENERAL mode is a map lookup per edge. (c) exact
  `reserve()` for filtration_/barycenters_/vertex_map_: entry and
  vertex-per-tet counts are exact functions of the partition shape, so
  a cheap pre-pass over mc_simplices gives capacities.
- [ ] **Iter 7 — POD filtration accumulator + dim buckets** (M, expect
  25-40% on 50k delaunay, bitwise identical): accumulate
  `{double value; int32_t v[3]; int32_t dim}` (every emitted simplex has
  <= 3 vertices: one chain level per distinct atom-set size in {2,3,4});
  sort+unique three per-dimension buckets (emit sites know dim);
  materialize the public Filtration once at the end with one allocation
  per surviving entry. Concatenating dim1|dim2|dim3 reproduces the
  global order exactly.
- [ ] **Iter 8 — bitmask faces + per-shape tables** (M-L, expect another
  20-30%, bitwise identical only if tables are generated by running the
  EXISTING enumerators on canonical shape representatives): faces are in
  bijection with multicolored atom subsets, so represent them as 4-bit
  masks; is_subset becomes mask AND; only 7 reachable partition shapes
  exist, so precompute face masks / inclusion pairs / chain triples per
  shape at static init. Allocation-free compute_chromatic_partition on
  fixed arrays reproducing the current (size desc, color asc) order.
- [ ] **Iter 9 — alpha classification** (M, alpha scenarios only,
  expected value-identical but NOT guaranteed bitwise — threshold
  predicates replace constructed radii; validate counts/filtration_sum
  and test_alpha_construction carefully): `Alpha_shape_3` GENERAL mode
  builds four maps + a spectrum we never use; `CGAL::Fixed_alpha_shape_3`
  classifies cells/facets/vertices as O(1) field reads. ~1.5s of the
  50k alpha collect is classification overhead.
- [ ] **Iter 10 (optional) — deterministic parallel subdivision** (L):
  chunk tets, thread-local vertex tables, merge in chunk order to
  reproduce global creation order exactly; only after iters 6-8.
- [ ] **Anytime, orthogonal — Julia wrapper allocations** (S):
  `flat[idx:idx+k-1] .+ Int32(1)` allocates twice per simplex (~13M
  allocations at 30k); replace with preallocated vector + explicit
  loop; same at vertex_atom_indices; sizehint the push! loop. Expect
  0.8-1.3s of the 2.6s wrapper_ms at 30k.

Investigated and rejected: cross-tet structural dedup (residual dim-2/3
duplication is only 1.2x after iter 6's singleton removal, not worth a
collection-contract change); barycenter math (cache-miss-only, ~0.1s);
SimplifiedSubdivision (not on a benchmarked path);
collect_multicolored_tetrahedra (0.13s at 50k, nothing left).

Deliberately skipped: pybind `def_readonly` re-conversion caching
(semantic change to attribute access), visualization-side performance
(cold path).

## Baseline reference (f2babcf, Darwin arm64, median of 3)

Worst-case pipelines: 50k/two delaunay total 13.9s (collect 0.4s — the
subdivision dominates); 50k/distinct alpha_uniform_half collect 2.1s.
Julia boundary at 30k: cxx 15.1s + wrapper 3.9s.
