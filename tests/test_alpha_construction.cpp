// Alpha Construction Tests on a regular tetrahedron
// See ALPHA_CONSTRUCTION_TESTS.md for detailed derivation and expected values.

#include <iostream>
#include <cassert>
#include <cmath>
#include <delaunay_interfaces/interface_generation.hpp>
#include <delaunay_interfaces/barycentric_subdivision.hpp>

using namespace delaunay_interfaces;

static const Points tet_points = {
    { 0.0,  0.0,       1.0 / std::sqrt(3.0)},          // v0: base front
    { 0.5,  0.0,      -1.0 / (2.0 * std::sqrt(3.0))},  // v1: base right-back
    {-0.5,  0.0,      -1.0 / (2.0 * std::sqrt(3.0))},  // v2: base left-back
    { 0.0,  std::sqrt(2.0 / 3.0), 0.0}                  // v3: apex
};

static const ColorLabels colors_31 = {1, 1, 1, 2};
static const ColorLabels colors_211 = {1, 1, 2, 3};
static const ColorLabels colors_22 = {1, 1, 2, 2};
static const ColorLabels colors_1111 = {1, 2, 3, 4};

struct FiltrationCounts {
    size_t vertices;
    size_t edges;
    size_t triangles;
    size_t total;
};

static FiltrationCounts count_filtration(const Filtration& f) {
    FiltrationCounts c = {0, 0, 0, 0};
    for (const auto& [simplex, val] : f) {
        if (simplex.size() == 1) c.vertices++;
        else if (simplex.size() == 2) c.edges++;
        else if (simplex.size() == 3) c.triangles++;
    }
    c.total = f.size();
    return c;
}

// ========================== 3-1 Coloring ==========================

void test_31_alpha_0() {
    std::cout << "Test: 3-1 alpha=0 (empty)\n";

    Radii radii(4, 0.0);

    auto [verts, filtration] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_31, radii, true, true);

    assert(verts.size() == 0);
    assert(filtration.size() == 0);

    std::cout << "  0 vertices, 0 simplices - PASS\n";
}

void test_31_alpha_051() {
    std::cout << "Test: 3-1 alpha=0.51 (edges only)\n";

    Radii radii(4, 0.51);

    auto [verts, filtration] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_31, radii, true, true);

    auto c = count_filtration(filtration);

    std::cout << "  " << verts.size() << " vertices, "
              << c.vertices << "v " << c.edges << "e " << c.triangles << "t = "
              << c.total << " simplices\n";

    assert(verts.size() == 3);
    assert(c.vertices == 3);
    assert(c.edges == 0);
    assert(c.triangles == 0);
    assert(c.total == 3);

    std::cout << "  PASS\n";
}

void test_31_alpha_058() {
    std::cout << "Test: 3-1 alpha=0.58 (edges + triangles)\n";

    Radii radii(4, 0.58);

    auto [verts, filtration] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_31, radii, true, true);

    auto c = count_filtration(filtration);

    std::cout << "  " << verts.size() << " vertices, "
              << c.vertices << "v " << c.edges << "e " << c.triangles << "t = "
              << c.total << " simplices\n";

    assert(verts.size() == 6);
    assert(c.vertices == 6);
    assert(c.edges == 6);
    assert(c.triangles == 0);
    assert(c.total == 12);

    std::cout << "  PASS\n";
}

void test_31_alpha_062() {
    std::cout << "Test: 3-1 alpha=0.62 (full tet)\n";

    Radii radii(4, 0.62);

    auto [verts_alpha, filt_alpha] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_31, radii, true, true);

    auto ca = count_filtration(filt_alpha);

    std::cout << "  Alpha:   " << verts_alpha.size() << " vertices, "
              << ca.vertices << "v " << ca.edges << "e " << ca.triangles << "t = "
              << ca.total << " simplices\n";

    auto [verts_del, filt_del] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_31, radii, true, false);

    auto cd = count_filtration(filt_del);

    std::cout << "  Delaunay: " << verts_del.size() << " vertices, "
              << cd.vertices << "v " << cd.edges << "e " << cd.triangles << "t = "
              << cd.total << " simplices\n";

    assert(verts_alpha.size() == 7);
    assert(ca.vertices == 7);
    assert(ca.edges == 12);
    assert(ca.triangles == 6);
    assert(ca.total == 25);

    assert(verts_alpha.size() == verts_del.size());
    assert(ca.total == cd.total);

    std::cout << "  PASS\n";
}

// ========================= 2-1-1 Coloring =========================

void test_211_alpha_052() {
    std::cout << "Test: 2-1-1 alpha=0.52 (edges only)\n";

    Radii radii(4, 0.52);

    auto [verts, filtration] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_211, radii, true, true);

    auto c = count_filtration(filtration);

    std::cout << "  " << verts.size() << " vertices, "
              << c.vertices << "v " << c.edges << "e " << c.triangles << "t = "
              << c.total << " simplices\n";

    assert(verts.size() == 5);
    assert(c.vertices == 5);
    assert(c.edges == 0);
    assert(c.triangles == 0);
    assert(c.total == 5);

    std::cout << "  PASS\n";
}

void test_211_alpha_058() {
    std::cout << "Test: 2-1-1 alpha=0.58 (edges + triangles)\n";

    Radii radii(4, 0.58);

    auto [verts, filtration] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_211, radii, true, true);

    auto c = count_filtration(filtration);

    std::cout << "  " << verts.size() << " vertices, "
              << c.vertices << "v " << c.edges << "e " << c.triangles << "t = "
              << c.total << " simplices\n";

    assert(verts.size() == 9);
    assert(c.vertices == 9);
    assert(c.edges == 10);
    assert(c.triangles == 0);
    assert(c.total == 19);

    std::cout << "  PASS\n";
}

void test_211_alpha_062() {
    std::cout << "Test: 2-1-1 alpha=0.62 (full tet)\n";

    Radii radii(4, 0.62);

    auto [verts_alpha, filt_alpha] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_211, radii, true, true);

    auto ca = count_filtration(filt_alpha);

    std::cout << "  Alpha:   " << verts_alpha.size() << " vertices, "
              << ca.vertices << "v " << ca.edges << "e " << ca.triangles << "t = "
              << ca.total << " simplices\n";

    auto [verts_del, filt_del] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_211, radii, true, false);

    auto cd = count_filtration(filt_del);

    std::cout << "  Delaunay: " << verts_del.size() << " vertices, "
              << cd.vertices << "v " << cd.edges << "e " << cd.triangles << "t = "
              << cd.total << " simplices\n";

    assert(verts_alpha.size() == 10);
    assert(ca.vertices == 10);
    assert(ca.edges == 19);
    assert(ca.triangles == 10);
    assert(ca.total == 39);

    assert(verts_alpha.size() == verts_del.size());
    assert(ca.total == cd.total);

    std::cout << "  PASS\n";
}

// ========================== 2-2 Coloring ==========================

void test_22_alpha_051() {
    std::cout << "Test: 2-2 alpha=0.51 (edges only)\n";

    Radii radii(4, 0.51);

    auto [verts, filtration] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_22, radii, true, true);

    auto c = count_filtration(filtration);

    std::cout << "  " << verts.size() << " vertices, "
              << c.vertices << "v " << c.edges << "e " << c.triangles << "t = "
              << c.total << " simplices\n";

    assert(verts.size() == 4);
    assert(c.vertices == 4);
    assert(c.edges == 0);
    assert(c.triangles == 0);
    assert(c.total == 4);

    std::cout << "  PASS\n";
}

void test_22_alpha_058() {
    std::cout << "Test: 2-2 alpha=0.58 (edges + triangles)\n";

    Radii radii(4, 0.58);

    auto [verts, filtration] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_22, radii, true, true);

    auto c = count_filtration(filtration);

    std::cout << "  " << verts.size() << " vertices, "
              << c.vertices << "v " << c.edges << "e " << c.triangles << "t = "
              << c.total << " simplices\n";

    assert(verts.size() == 8);
    assert(c.vertices == 8);
    assert(c.edges == 8);
    assert(c.triangles == 0);
    assert(c.total == 16);

    std::cout << "  PASS\n";
}

void test_22_alpha_062() {
    std::cout << "Test: 2-2 alpha=0.62 (full tet)\n";

    Radii radii(4, 0.62);

    auto [verts_alpha, filt_alpha] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_22, radii, true, true);

    auto ca = count_filtration(filt_alpha);

    std::cout << "  Alpha:   " << verts_alpha.size() << " vertices, "
              << ca.vertices << "v " << ca.edges << "e " << ca.triangles << "t = "
              << ca.total << " simplices\n";

    auto [verts_del, filt_del] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_22, radii, true, false);

    auto cd = count_filtration(filt_del);

    std::cout << "  Delaunay: " << verts_del.size() << " vertices, "
              << cd.vertices << "v " << cd.edges << "e " << cd.triangles << "t = "
              << cd.total << " simplices\n";

    assert(verts_alpha.size() == 9);
    assert(ca.vertices == 9);
    assert(ca.edges == 16);
    assert(ca.triangles == 8);
    assert(ca.total == 33);

    assert(verts_alpha.size() == verts_del.size());
    assert(ca.total == cd.total);

    std::cout << "  PASS\n";
}

// ======================== 1-1-1-1 Coloring ========================

// alpha = 0: nothing in the alpha complex
void test_1111_alpha_0() {
    std::cout << "Test: 1-1-1-1 alpha=0 (empty)\n";

    // r=0 → weight=0, all circumradii² > 0 → nothing passes
    Radii radii(4, 0.0);

    auto [verts, filtration] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_1111, radii, true, true);

    assert(verts.size() == 0);
    assert(filtration.size() == 0);

    std::cout << "  0 vertices, 0 simplices - PASS\n";
}

// alpha = 0.51: only edges in alpha complex (α²=0.25 < 0.51²=0.2601)
// 6 free bicolored edges → 6 barycenters, 6 vertex simplices
void test_1111_alpha_051() {
    std::cout << "Test: 1-1-1-1 alpha=0.51 (edges only)\n";

    Radii radii(4, 0.51);

    auto [verts, filtration] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_1111, radii, true, true);

    auto c = count_filtration(filtration);

    std::cout << "  " << verts.size() << " vertices, "
              << c.vertices << "v " << c.edges << "e " << c.triangles << "t = "
              << c.total << " simplices\n";

    assert(verts.size() == 6);
    assert(c.vertices == 6);
    assert(c.edges == 0);
    assert(c.triangles == 0);
    assert(c.total == 6);

    std::cout << "  PASS\n";
}

// alpha = 0.58: edges + triangles in alpha complex (α²=1/3 < 0.58²=0.3364)
// 4 free tricolored triangles → 10 barycenters (6 shared edge + 4 tri centers)
// 12 edges (3 per triangle × 4), no triangles in filtration
void test_1111_alpha_058() {
    std::cout << "Test: 1-1-1-1 alpha=0.58 (edges + triangles)\n";

    Radii radii(4, 0.58);

    auto [verts, filtration] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_1111, radii, true, true);

    auto c = count_filtration(filtration);

    std::cout << "  " << verts.size() << " vertices, "
              << c.vertices << "v " << c.edges << "e " << c.triangles << "t = "
              << c.total << " simplices\n";

    assert(verts.size() == 10);
    assert(c.vertices == 10);
    assert(c.edges == 12);
    assert(c.triangles == 0);
    assert(c.total == 22);

    std::cout << "  PASS\n";
}

// alpha = 0.62: full tet in alpha complex (α²=3/8 < 0.62²=0.3844)
// Same as Delaunay: 11 barycenters, 22 edges, 12 triangles
void test_1111_alpha_062() {
    std::cout << "Test: 1-1-1-1 alpha=0.62 (full tet)\n";

    Radii radii(4, 0.62);

    auto [verts_alpha, filt_alpha] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_1111, radii, true, true);

    auto ca = count_filtration(filt_alpha);

    std::cout << "  Alpha:   " << verts_alpha.size() << " vertices, "
              << ca.vertices << "v " << ca.edges << "e " << ca.triangles << "t = "
              << ca.total << " simplices\n";

    // Compare with Delaunay (no alpha filtering)
    auto [verts_del, filt_del] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_1111, radii, true, false);

    auto cd = count_filtration(filt_del);

    std::cout << "  Delaunay: " << verts_del.size() << " vertices, "
              << cd.vertices << "v " << cd.edges << "e " << cd.triangles << "t = "
              << cd.total << " simplices\n";

    assert(verts_alpha.size() == 11);
    assert(ca.vertices == 11);
    assert(ca.edges == 22);
    assert(ca.triangles == 12);
    assert(ca.total == 45);

    // Must match Delaunay exactly
    assert(verts_alpha.size() == verts_del.size());
    assert(ca.total == cd.total);

    std::cout << "  PASS\n";
}

// =================== Bipyramid — 3-1 × 2 Doubling ===================
// v4 at -apex, colors {1,1,1,2,2}: two disjoint 3-1 tets → double counts

static const Points bipyramid_points = {
    { 0.0,  0.0,       1.0 / std::sqrt(3.0)},          // v0: base front
    { 0.5,  0.0,      -1.0 / (2.0 * std::sqrt(3.0))},  // v1: base right-back
    {-0.5,  0.0,      -1.0 / (2.0 * std::sqrt(3.0))},  // v2: base left-back
    { 0.0,  std::sqrt(2.0 / 3.0), 0.0},                 // v3: apex (top)
    { 0.0, -std::sqrt(2.0 / 3.0), 0.0}                  // v4: -apex (bottom)
};

static const ColorLabels colors_bipyramid = {1, 1, 1, 2, 2};

void test_bipyramid_alpha_051() {
    std::cout << "Test: bipyramid 3-1x2 alpha=0.51 (edges only)\n";

    Radii radii(5, 0.51);

    auto [verts, filtration] = get_barycentric_subdivision_and_filtration(
        bipyramid_points, colors_bipyramid, radii, true, true);

    auto c = count_filtration(filtration);

    std::cout << "  " << verts.size() << " vertices, "
              << c.vertices << "v " << c.edges << "e " << c.triangles << "t = "
              << c.total << " simplices\n";

    assert(verts.size() == 6);
    assert(c.vertices == 6);
    assert(c.edges == 0);
    assert(c.triangles == 0);
    assert(c.total == 6);

    std::cout << "  PASS\n";
}

void test_bipyramid_alpha_058() {
    std::cout << "Test: bipyramid 3-1x2 alpha=0.58 (edges + triangles)\n";

    Radii radii(5, 0.58);

    auto [verts, filtration] = get_barycentric_subdivision_and_filtration(
        bipyramid_points, colors_bipyramid, radii, true, true);

    auto c = count_filtration(filtration);

    std::cout << "  " << verts.size() << " vertices, "
              << c.vertices << "v " << c.edges << "e " << c.triangles << "t = "
              << c.total << " simplices\n";

    assert(verts.size() == 12);
    assert(c.vertices == 12);
    assert(c.edges == 12);
    assert(c.triangles == 0);
    assert(c.total == 24);

    std::cout << "  PASS\n";
}

void test_bipyramid_alpha_062() {
    std::cout << "Test: bipyramid 3-1x2 alpha=0.62 (full tets)\n";

    Radii radii(5, 0.62);

    auto [verts_alpha, filt_alpha] = get_barycentric_subdivision_and_filtration(
        bipyramid_points, colors_bipyramid, radii, true, true);

    auto ca = count_filtration(filt_alpha);

    std::cout << "  Alpha:   " << verts_alpha.size() << " vertices, "
              << ca.vertices << "v " << ca.edges << "e " << ca.triangles << "t = "
              << ca.total << " simplices\n";

    auto [verts_del, filt_del] = get_barycentric_subdivision_and_filtration(
        bipyramid_points, colors_bipyramid, radii, true, false);

    auto cd = count_filtration(filt_del);

    std::cout << "  Delaunay: " << verts_del.size() << " vertices, "
              << cd.vertices << "v " << cd.edges << "e " << cd.triangles << "t = "
              << cd.total << " simplices\n";

    assert(verts_alpha.size() == 14);
    assert(ca.vertices == 14);
    assert(ca.edges == 24);
    assert(ca.triangles == 12);
    assert(ca.total == 50);

    assert(verts_alpha.size() == verts_del.size());
    assert(ca.total == cd.total);

    std::cout << "  PASS\n";
}

int main() {
    std::cout << "Alpha Construction Tests\n";
    std::cout << "========================\n\n";

    try {
        std::cout << "--- 3-1 Coloring ---\n";
        test_31_alpha_0();
        test_31_alpha_051();
        test_31_alpha_058();
        test_31_alpha_062();

        std::cout << "\n--- 2-1-1 Coloring ---\n";
        test_211_alpha_052();
        test_211_alpha_058();
        test_211_alpha_062();

        std::cout << "\n--- 2-2 Coloring ---\n";
        test_22_alpha_051();
        test_22_alpha_058();
        test_22_alpha_062();

        std::cout << "\n--- 1-1-1-1 Coloring ---\n";
        test_1111_alpha_0();
        test_1111_alpha_051();
        test_1111_alpha_058();
        test_1111_alpha_062();

        std::cout << "\n--- Bipyramid 3-1x2 ---\n";
        test_bipyramid_alpha_051();
        test_bipyramid_alpha_058();
        test_bipyramid_alpha_062();

        std::cout << "\nAll alpha construction tests passed!\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nTest failed with exception: " << e.what() << "\n";
        return 1;
    }
}
