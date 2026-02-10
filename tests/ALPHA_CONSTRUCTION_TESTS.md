# Alpha Construction Tests

## Regular Tetrahedron (edge length 1)

Base centered at origin, apex on +y.

| Vertex | Label | x | y | z |
|--------|-------|-----|-----------|-----------|
| v0 | base front | 0.0 | 0.0 | 1/√3 ≈ 0.57735 |
| v1 | base right | 0.5 | 0.0 | -1/(2√3) ≈ -0.28868 |
| v2 | base left | -0.5 | 0.0 | -1/(2√3) ≈ -0.28868 |
| v3 | apex | 0.0 | √(2/3) ≈ 0.81650 | 0.0 |

## Alpha Filtration Thresholds

| Simplex | Count | α² | α |
|-------------|-------|--------------|---------|
| Vertices | 4 | 0 | 0 |
| Edges | 6 | 1/4 = 0.25 | 0.5 |
| Triangles | 4 | 1/3 ≈ 0.3333 | ≈ 0.5774 |
| Tetrahedron | 1 | 3/8 = 0.375 | ≈ 0.6124 |

Setting all radii to `r` simulates an unweighted alpha complex at threshold `r`: a simplex with circumradius² = C enters the weighted alpha complex when `C - r² ≤ 0`, i.e. when `r ≥ √C`.

---

## 3-1 Coloring

Colors: `{1, 1, 1, 2}` — v0, v1, v2 share one color, v3 has a different color.

Multicolored simplices are those containing at least one vertex from each color group. The monocolored triangle `{v0, v1, v2}` and its 3 edges are **not** multicolored.

### r = 0 — Empty

Nothing in the alpha complex.

→ **0 barycenters, 0 simplices**

### r = 0.51 — Edges Only

6 edges enter the alpha complex (0.25 < 0.51² = 0.2601). Only 3 are multicolored: `{v0,v3}`, `{v1,v3}`, `{v2,v3}`. The 3 monocolored edges `{v0,v1}`, `{v0,v2}`, `{v1,v2}` are skipped.

Each free bicolored edge produces **1 barycenter** at the midpoint.

→ **3 barycenters, 3v = 3 simplices**

### r = 0.58 — Edges + Triangles

6 edges + 4 triangles enter the alpha complex (1/3 < 0.58² = 0.3364). No tetrahedron.

Multicolored triangles (3 of 4):
- `{v0, v1, v3}` — partition `{{v0,v1}, {v3}}` (2-1)
- `{v0, v2, v3}` — partition `{{v0,v2}, {v3}}` (2-1)
- `{v1, v2, v3}` — partition `{{v1,v2}, {v3}}` (2-1)

The triangle `{v0, v1, v2}` is monocolored → skipped.

Each free 2-1 triangle (e.g. `{v0, v1, v3}` with partition `{{v0,v1}, {v3}}`) produces:
- **2 level-2 combinations**: `{{v0},{v3}}`, `{{v1},{v3}}` (edge midpoints)
- **1 level-3 combination**: `{{v0,v1},{v3}}` (triangle center)

Edge midpoints are shared across triangles (e.g. mid(v0,v3) appears in both `{v0,v1,v3}` and `{v0,v2,v3}`):
- Unique barycenters: 3 edge midpoints + 3 triangle centers = **6**

Filtration edges from subset inclusion (level-2 ⊂ level-3) within each triangle:
- 2 edges per triangle × 3 triangles = **6 edges**
- All distinct (each triangle center is unique)

No filtration triangles: only 2 levels per free triangle.

→ **6 barycenters, 6v + 6e = 12 simplices**

### r = 0.62 — Full Tetrahedron

The tetrahedron enters the alpha complex (3/8 < 0.62² = 0.3844). Processed as 3-1 partition `{{v0,v1,v2}, {v3}}`.

Combinations (7 total):
- **Level 2** (3): `{{v0},{v3}}`, `{{v1},{v3}}`, `{{v2},{v3}}`
- **Level 3** (3): `{{v0,v1},{v3}}`, `{{v0,v2},{v3}}`, `{{v1,v2},{v3}}`
- **Level 4** (1): `{{v0,v1,v2},{v3}}`

Edges from subset inclusion (12 total):
- Level 2 → Level 3: each level-3 contains 2 level-2 combos → 3 × 2 = **6**
- Level 2 → Level 4: all 3 level-2 combos ⊂ level-4 → **3**
- Level 3 → Level 4: all 3 level-3 combos ⊂ level-4 → **3**

Maximal chains (level 2 → 3 → 4): 3 level-2 × 2 paths each = **6 triangles**

Equivalent to the Delaunay (no alpha filtering) result.

→ **7 barycenters, 7v + 12e + 6t = 25 simplices**

---

## 2-1-1 Coloring

Colors: `{1, 1, 2, 3}` — v0, v1 are color 1; v2 is color 2; v3 is color 3.

Only the edge `{v0, v1}` is monocolored. All other simplices are multicolored.

### r = 0.52 — Edges Only

6 edges in the alpha complex (0.25 < 0.52² = 0.2704). 5 are multicolored: `{v0,v2}`, `{v0,v3}`, `{v1,v2}`, `{v1,v3}`, `{v2,v3}`. The monocolored edge `{v0,v1}` is skipped.

Each free bicolored edge produces **1 barycenter** at the midpoint.

→ **5 barycenters, 5v = 5 simplices**

### r = 0.58 — Edges + Triangles

6 edges + 4 triangles in the alpha complex (1/3 < 0.58² = 0.3364). No tetrahedron. All 4 triangles are multicolored, but with different partition types:

| Triangle | Partition | Type |
|----------|-----------|------|
| `{v0, v1, v2}` | `{{v0,v1}, {v2}}` | 2-1 |
| `{v0, v1, v3}` | `{{v0,v1}, {v3}}` | 2-1 |
| `{v0, v2, v3}` | `{{v0}, {v2}, {v3}}` | 1-1-1 |
| `{v1, v2, v3}` | `{{v1}, {v2}, {v3}}` | 1-1-1 |

Each 2-1 triangle produces 2 level-2 + 1 level-3 combinations.
Each 1-1-1 triangle produces 3 level-2 + 1 level-3 combinations.

Edge midpoints are shared across triangles (each multicolored edge appears in 2 of the 4 triangles):
- Unique barycenters: 5 edge midpoints + 4 triangle centers = **9**

Filtration edges from subset inclusion (level-2 ⊂ level-3):
- 2 per 2-1 triangle × 2 + 3 per 1-1-1 triangle × 2 = **10 edges**

No filtration triangles: only 2 levels per free triangle.

→ **9 barycenters, 9v + 10e = 19 simplices**

### r = 0.62 — Full Tetrahedron

The tetrahedron enters the alpha complex. Processed as 2-1-1 partition `{{v0,v1}, {v2}, {v3}}`.

Combinations (10 total):
- **Level 2** (5): `{{v0},{v2}}`, `{{v1},{v2}}`, `{{v0},{v3}}`, `{{v1},{v3}}`, `{{v2},{v3}}`
- **Level 3** (4): `{{v0,v1},{v2}}`, `{{v0,v1},{v3}}`, `{{v0},{v2},{v3}}`, `{{v1},{v2},{v3}}`
- **Level 4** (1): `{{v0,v1},{v2},{v3}}`

Edges from subset inclusion (19 total):
- Level 2 → Level 3: each level-2 is contained in 2 level-3 combos → 5 × 2 = **10**
- Level 2 → Level 4: all 5 level-2 combos ⊂ level-4 → **5**
- Level 3 → Level 4: all 4 level-3 combos ⊂ level-4 → **4**

Maximal chains (level 2 → 3 → 4): 5 level-2 × 2 paths each = **10 triangles**

Equivalent to the Delaunay (no alpha filtering) result.

→ **10 barycenters, 10v + 19e + 10t = 39 simplices**

---

## 2-2 Coloring

Colors: `{1, 1, 2, 2}` — v0, v1 are color 1; v2, v3 are color 2.

Multicolored simplices must contain vertices of both colors. The monocolored edges `{v0,v1}` and `{v2,v3}` are skipped.

### r = 0.51 — Edges Only

6 edges in the alpha complex. 4 are multicolored: `{v0,v2}`, `{v0,v3}`, `{v1,v2}`, `{v1,v3}`. The 2 monocolored edges `{v0,v1}`, `{v2,v3}` are skipped.

Each free bicolored edge produces **1 barycenter** at the midpoint.

→ **4 barycenters, 4v = 4 simplices**

### r = 0.58 — Edges + Triangles

6 edges + 4 triangles in the alpha complex. All 4 triangles are multicolored (each is a 2-1 partition):

| Triangle | Partition |
|----------|-----------|
| `{v0, v1, v2}` | `{{v0,v1}, {v2}}` |
| `{v0, v1, v3}` | `{{v0,v1}, {v3}}` |
| `{v0, v2, v3}` | `{{v0}, {v2,v3}}` |
| `{v1, v2, v3}` | `{{v1}, {v2,v3}}` |

Each free 2-1 triangle produces 2 level-2 (edge midpoints) + 1 level-3 (triangle center) combinations.

Edge midpoints are shared across triangles (each multicolored edge appears in 2 of the 4 triangles):
- Unique barycenters: 4 edge midpoints + 4 triangle centers = **8**

Filtration edges from subset inclusion (level-2 ⊂ level-3) within each triangle:
- 2 edges per triangle × 4 triangles = **8 edges**

No filtration triangles: only 2 levels per free triangle.

→ **8 barycenters, 8v + 8e = 16 simplices**

### r = 0.62 — Full Tetrahedron

The tetrahedron enters the alpha complex. Processed as 2-2 partition `{{v0,v1}, {v2,v3}}`.

Combinations (9 total) — cross-product of non-empty subsets of each part:
- **Level 2** (4): `{{v0},{v2}}`, `{{v0},{v3}}`, `{{v1},{v2}}`, `{{v1},{v3}}`
- **Level 3** (4): `{{v0},{v2,v3}}`, `{{v1},{v2,v3}}`, `{{v0,v1},{v2}}`, `{{v0,v1},{v3}}`
- **Level 4** (1): `{{v0,v1},{v2,v3}}`

Edges from subset inclusion (16 total):
- Level 2 → Level 3: each level-2 combo is contained in 2 level-3 combos → 4 × 2 = **8**
- Level 2 → Level 4: all 4 level-2 combos ⊂ level-4 → **4**
- Level 3 → Level 4: all 4 level-3 combos ⊂ level-4 → **4**

Maximal chains (level 2 → 3 → 4): 4 level-2 × 2 paths each = **8 triangles**

Equivalent to the Delaunay (no alpha filtering) result.

→ **9 barycenters, 9v + 16e + 8t = 33 simplices**

---

## 1-1-1-1 Coloring

Each vertex gets a distinct color, so **every simplex is multicolored**. The barycentric subdivision of the interface depends on which simplices are present in the alpha complex.

### r = 0 — Empty

Nothing in the alpha complex.

→ **0 barycenters, 0 simplices**

### r = 0.51 — Edges Only

6 edges enter the alpha complex (0.25 < 0.51² = 0.2601). No triangles or tetrahedra.

Each free bicolored edge `{vi, vj}` produces **1 barycenter** at the midpoint. The 6 edges are independent — no subset inclusion relationships between them.

→ **6 barycenters, 6v = 6 simplices**

### r = 0.58 — Edges + Triangles

6 edges + 4 triangles enter the alpha complex (1/3 < 0.58² = 0.3364). No tetrahedron. The 4 triangles are **free** (not faces of any alpha tetrahedron). The 6 edges are covered by the free triangles.

Each free 1-1-1 triangle (e.g. `{v0, v1, v2}`) produces:
- **3 level-2 combinations** (edge midpoints): `{v0,v1}`, `{v0,v2}`, `{v1,v2}`
- **1 level-3 combination** (triangle center): `{v0,v1,v2}`

Edge midpoints are shared across triangles (each edge appears in 2 of the 4 triangles), so:
- Unique barycenters: 6 edge midpoints + 4 triangle centers = **10**

Filtration edges come from subset inclusion (level-2 ⊂ level-3) within each triangle:
- 3 edges per triangle × 4 triangles = **12 edges**
- All distinct (each triangle center is unique)

No filtration triangles: only 2 levels per free triangle, so maximal chains have length 2 (= edges, not triangles).

→ **10 barycenters, 10v + 12e = 22 simplices**

### r = 0.62 — Full Tetrahedron

The tetrahedron enters the alpha complex (3/8 < 0.62² = 0.3844). Single multicolored tet processed as 1-1-1-1 partition.

Combinations (11 total):
- **Level 2** (6): all C(4,2) = 6 edge midpoints
- **Level 3** (4): all C(4,3) = 4 triangle centers
- **Level 4** (1): tetrahedron center

Edges from subset inclusion (22 total):
- Level 2 → Level 3: each of 4 level-3 combos contains 3 level-2 combos → **12**
- Level 2 → Level 4: all 6 level-2 combos ⊂ level-4 → **6**
- Level 3 → Level 4: all 4 level-3 combos ⊂ level-4 → **4**

Maximal chains (level 2 → 3 → 4): 4 × 3 = **12 triangles**

Equivalent to the Delaunay (no alpha filtering) result.

→ **11 barycenters, 11v + 22e + 12t = 45 simplices**
