# The filtration bucket sort, by example

`src/filtration_sort.hpp` (`detail::sort_bucket`, used by
`BarycentricSubdivision::finalize_filtration`) sorts each per-dimension
filtration bucket by `(value, ids lexicographic)`. The buckets are
`FlatFiltration`'s dim0/dim1/dim2, named by geometric dimension: dim0
holds the subdivision's vertex singletons, dim1 its edges (inclusion
pairs), dim2 its triangles (maximal chains) — the top cells of the
2-dimensional interface surface. That total order makes the
output independent of emission order (determinism) and makes exact duplicates
from shared tetrahedron faces adjacent, so `std::unique` can drop them.

Buckets of fewer than 10,000 entries just use `std::sort`. Larger buckets use
a stable LSD radix sort on the value's IEEE-754 bit pattern, followed by an
id-order fixup of equal-value runs. `benchmarks/benchmark_radix_sort.cpp`
compares the two on real pipeline buckets:

```bash
cmake -S . -B build -DBUILD_BENCHMARKS=ON && cmake --build build --target benchmark_radix_sort -j4
./build/benchmarks/benchmark_radix_sort            # defaults: --sizes 20000,50000 --repeats 5
```

Representative results (Darwin arm64, median of 7): radix is 2.8–3.0x faster
than `std::sort` on dim1 (the largest bucket; 139 ms vs 412 ms on 5.4M
entries), 2.4–2.6x on dim2, and 5.8–6.4x on dim0, whose 16-byte entries move
cheapest. The benchmark verifies both algorithms produce byte-identical
output on every repeat.

## Worked example

Take a tiny dim1 bucket (edge entries `{value, v[2]}`) in emission order,
with one exact duplicate from a face shared by two adjacent tetrahedra:

```
#   value    ids
0   0.25    (7, 9)
1   0.50    (2, 3)
2   0.25    (4, 8)
3   0.125   (5, 6)
4   0.50    (2, 3)    <- duplicate of #1 (shared face)
5   0.25    (1, 2)
```

(A bucket this small would take the `std::sort` fallback; pretend it goes
down the radix path.)

### Step 1: values become integers

Each `double` is reinterpreted as its raw 64-bit pattern. Filtration values
are mean pairwise distances — never negative, never NaN — and for
non-negative IEEE-754 doubles, larger value means larger bit pattern read as
an unsigned integer. That is what makes integer sorting legal on floats:

```
0.125 = 0x3FC0'0000'0000'0000
0.25  = 0x3FD0'0000'0000'0000
0.50  = 0x3FE0'0000'0000'0000
```

`0x3FC < 0x3FD < 0x3FE` — the integer order *is* the numeric order.

### Step 2: LSD radix — eight counting-sort passes, one per byte

Each pass is a stable counting sort on one byte, least significant byte
first. Every pass starts by building the byte's histogram: one linear read
counting how many entries carry each of the 256 possible byte values.

That histogram is needed for the scatter anyway, but it also enables a
skip: **if a single byte value's count equals n, every entry shares that
byte**. A stable counting sort with all n entries in one 256-bucket then
emits them in their current order — the pass is the identity permutation —
so the expensive part (scattering n entries to a scratch buffer) is skipped
and only the cheap histogram was paid.

In this example that skip does almost all the work: bytes 0–5 are `00` for
every entry and byte 7 is `3F` for every entry, so seven of the eight passes
are skipped. Only the byte-6 pass (`C0`/`D0`/`E0`) scatters:

1. *Count:* `C0` appears 1x, `D0` 3x, `E0` 2x.
2. *Prefix sums → start offsets:* `C0` at slot 0, `D0` at slot 1, `E0` at slot 4.
3. *Stable scatter* (walk input in order, drop each entry at its byte's next
   free slot):

```
0.125  (5,6)
0.25   (7,9)   <- the three 0.25s are still in EMISSION order:
0.25   (4,8)      the radix sorted by value only
0.25   (1,2)
0.50   (2,3)
0.50   (2,3)
```

Real data behaves the same way, just less extremely: distances in one point
cloud cluster in a narrow exponent range, so the high bytes (sign + exponent
+ top mantissa bits) are often shared and their passes get skipped, while
the noisy low mantissa bytes do the scattering.

### Step 3: fix the ties

One linear scan finds runs of equal values (equal means bit-identical —
duplicates come from identical arithmetic) and re-sorts each run by ids with
a small `std::sort`. Runs are typically a handful of entries, so this is
cheap:

```
0.125  (5,6)
0.25   (1,2)
0.25   (4,8)
0.25   (7,9)
0.50   (2,3)
0.50   (2,3)
```

The result now matches `std::sort` with the full `(value, ids)` comparator
exactly — proven at scale by the byte-identical A/B dumps in the
optimization series (iter 11) and re-verified on every benchmark run.

### Step 4: dedup (back in `finalize_filtration`)

The two `0.50 (2,3)` entries are adjacent, so `std::unique` collapses them:

```
0.125  (5,6)
0.25   (1,2)
0.25   (4,8)
0.25   (7,9)
0.50   (2,3)
```

This whole pipeline runs once per dimension bucket (dim0/dim1/dim2); above
50,000 total entries the three buckets are finalized concurrently
(`std::async`, iter 14) — each thread runs the sequential algorithm on its
own disjoint bucket, so the parallel output is bitwise identical to the
sequential one. That measured 1.65–1.68x on the finalization phase (~13–17%
of end-to-end pipeline time at 50k).
