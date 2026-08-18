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

- [x] **Iter 6 — singleton emission + alpha predicate order + reserves**
  (committed, fdd23cd; invariants bitwise identical; delaunay totals
  1.20-1.29x, alpha ~1.1x, Julia cxx at 30k 1.28x, vs ~1.05-1.1x
  collect drift on untouched paths — most of the gain is real; per-shape
  entry counts corrected: tet 2-1-1 pushes 39 raw entries, not 38):
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
- [x] **Iter 7 — POD filtration accumulator + dim buckets** (committed,
  d33b49f; invariants bitwise identical; subdivision median 1.51x vs
  ~1.1x untouched-path drift, delaunay totals 1.4-1.6x, Julia cxx at
  30k 1.39x and wrapper 1.7x as a side effect of the cached
  materialized filtration; reserve() now takes per-dim counts): accumulate
  `{double value; int32_t v[3]; int32_t dim}` (every emitted simplex has
  <= 3 vertices: one chain level per distinct atom-set size in {2,3,4});
  sort+unique three per-dimension buckets (emit sites know dim);
  materialize the public Filtration once at the end with one allocation
  per surviving entry. Concatenating dim1|dim2|dim3 reproduces the
  global order exactly.
- [x] **Iter 8 — per-shape lookup tables** (committed 294de84, recorded;
  invariants bitwise identical; subdivision median 3.98x vs 1.1x
  untouched-path drift — far beyond the 20-30% estimate; 50k/two
  delaunay total 3.19s -> 1.05s, Julia cxx at 30k 3.9s -> 0.99s): faces are in
  bijection with multicolored atom subsets, so represent them as 4-bit
  masks; is_subset becomes mask AND; only 7 reachable partition shapes
  exist, so precompute face masks / inclusion pairs / chain triples per
  shape at static init. Allocation-free compute_chromatic_partition on
  fixed arrays reproducing the current (size desc, color asc) order.
- [x] **Iter 9 — alpha classification** (committed 5cd80d8, recorded;
  invariants exactly identical at ALL sizes incl. 50k — the feared
  threshold flips did not materialize on these inputs; alpha collect
  1.6-2.4x, 50k alpha totals 1.3-1.65s -> 0.70-0.86s; untouched
  delaunay drift 0.999 — cleanest run of the series): `Alpha_shape_3` GENERAL mode
  builds four maps + a spectrum we never use; `CGAL::Fixed_alpha_shape_3`
  classifies cells/facets/vertices as O(1) field reads. ~1.5s of the
  50k alpha collect is classification overhead.
- [-] **Iter 10 — DROPPED as specified** (2026-08-18 round-2
  investigation): its target, the process_simplex loop, shrank to
  ~170-265ms at 50k after iters 6-8; realistic net win ~150ms for
  effort L, superseded by iters 11/14/16 below at far lower cost.

Round-2 queue (2026-08-18, from a MEASURED profile: phase
instrumentation + sample call graphs + a radix prototype at <= 20k,
extrapolated by simplex count; instrumented totals matched the e5355a1
record). Post-iter-9 profile at 50k: subdivision is 57-61% per-bucket
std::sort; alpha collect is 61% CGAL initialize_alpha (Gabriel
predicate on ~96% of facets before the cheap radius test, plus a
global std::map over all edges); the public vector-of-vectors
Filtration is built on binding paths only to be re-flattened.

- [ ] **Iter 11 — radix-sort the finalization buckets** (S/M, bitwise
  identical, save ~650-700ms distinct / ~370ms two at 50k): stable LSD
  radix on the value's IEEE bit pattern (values >= 0, so bit order =
  numeric order; assert non-negativity), then order equal-value runs by
  ids. Prototype: 7.3M dim2-like entries, std::sort 535ms -> 112ms,
  output verified identical.
- [ ] **Iter 12 — custom fixed-alpha classification** (M/L, expected
  bitwise identical - same exact-predicate set, reordered cheap-first;
  weighted variant must replicate the weighted predicates exactly;
  save ~300-400ms on every alpha path): build the triangulation with
  the Fixed_alpha bases but skip the ctor's initialize_alpha; classify
  cells as CGAL does, facets radius-test-first (Gabriel only when
  needed), edges without the global map - status computed only for
  multicolored edges, Gabriel only when all incident facets are
  EXTERIOR.
- [ ] **Iter 13 — flat internal filtration, lazy public Filtration**
  (M, C++ struct field becomes accessor or additive API; save
  ~90-155ms C++ and ~170-250ms off the 30k Julia boundary): bindings
  consume the flat form directly; vector-of-vectors materialized only
  for C++ callers on demand.
- [ ] **Iter 14 — parallel finalization** (S/M, deterministic by
  construction - independent buckets on 3 threads; after iter 11; save
  ~60-100ms).
- [ ] **Iter 15 — vertex_atom_indices exact reserve** (S, bitwise
  identical, ~45ms distinct).
- [ ] **Iter 16 (optional) — open-addressing vertex map** (M, bitwise
  identical, ~80-120ms distinct incl. the ~50ms unordered_map
  teardown).

Round-2 rejected: iter 10 as written (above); CGAL Parallel_tag
triangulation (~85ms, permutes vertex ids); further Julia decode
micro-opts (measured at the one-allocation-per-simplex floor of the
public field type - only a breaking CSR field change helps, ~45-55% of
wrapper_ms is inherent); capacity pre-pass / uniques / delaunay
collect / alpha tet+facet loops (all <= ~35ms). Realistic endpoint at
50k: two-delaunay ~0.55s, distinct ~0.85s, alpha ~0.4-0.5s.
- [x] **Anytime, orthogonal — Julia wrapper allocations** (committed
  e5355a1, recorded): preallocated decode loops replace the
  slice+broadcast pattern; wrapper_ms at 30k 1.14s -> 0.64s (1.78x),
  10k 1.49x; the 1k wrapper number is GC/JIT-jitter-dominated and not
  meaningful. Element/type-identical output verified in-process.
  Observation from verification: the C++ alpha pipeline's RAW simplex
  ordering differs across processes (address-dependent CGAL internals);
  all invariants and the value-sorted filtration stay identical, but
  cross-process raw diffs are not a usable verification tool.

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
