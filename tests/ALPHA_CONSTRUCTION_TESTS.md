# Alpha Construction Tests — 1-1-1-1 Split

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
