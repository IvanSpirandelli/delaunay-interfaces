// Alpha Construction Tests on a regular tetrahedron
// See ALPHA_CONSTRUCTION_TESTS.md for detailed derivation and expected values.

#include <iostream>
#include <cassert>
#include <cmath>
#include <random>
#include <set>
#include <map>
#include <fstream>
#include <filesystem>
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

    auto [verts, filtration, mc_simplices, vai] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_31, radii, true, true);

    assert(verts.size() == 0);
    assert(filtration.size() == 0);

    std::cout << "  0 vertices, 0 simplices - PASS\n";
}

void test_31_alpha_051() {
    std::cout << "Test: 3-1 alpha=0.51 (edges only)\n";

    Radii radii(4, 0.51);

    auto [verts, filtration, mc_simplices, vai] = get_barycentric_subdivision_and_filtration(
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

    auto [verts, filtration, mc_simplices, vai] = get_barycentric_subdivision_and_filtration(
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

    auto [verts_alpha, filt_alpha, mc_alpha, vai_alpha] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_31, radii, true, true);

    auto ca = count_filtration(filt_alpha);

    std::cout << "  Alpha:   " << verts_alpha.size() << " vertices, "
              << ca.vertices << "v " << ca.edges << "e " << ca.triangles << "t = "
              << ca.total << " simplices\n";

    auto [verts_del, filt_del, mc_del, vai_del] = get_barycentric_subdivision_and_filtration(
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

    auto [verts, filtration, mc_simplices, vai] = get_barycentric_subdivision_and_filtration(
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

    auto [verts, filtration, mc_simplices, vai] = get_barycentric_subdivision_and_filtration(
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

    auto [verts_alpha, filt_alpha, mc_alpha, vai_alpha] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_211, radii, true, true);

    auto ca = count_filtration(filt_alpha);

    std::cout << "  Alpha:   " << verts_alpha.size() << " vertices, "
              << ca.vertices << "v " << ca.edges << "e " << ca.triangles << "t = "
              << ca.total << " simplices\n";

    auto [verts_del, filt_del, mc_del, vai_del] = get_barycentric_subdivision_and_filtration(
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

    auto [verts, filtration, mc_simplices, vai] = get_barycentric_subdivision_and_filtration(
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

    auto [verts, filtration, mc_simplices, vai] = get_barycentric_subdivision_and_filtration(
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

    auto [verts_alpha, filt_alpha, mc_alpha, vai_alpha] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_22, radii, true, true);

    auto ca = count_filtration(filt_alpha);

    std::cout << "  Alpha:   " << verts_alpha.size() << " vertices, "
              << ca.vertices << "v " << ca.edges << "e " << ca.triangles << "t = "
              << ca.total << " simplices\n";

    auto [verts_del, filt_del, mc_del, vai_del] = get_barycentric_subdivision_and_filtration(
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

    auto [verts, filtration, mc_simplices, vai] = get_barycentric_subdivision_and_filtration(
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

    auto [verts, filtration, mc_simplices, vai] = get_barycentric_subdivision_and_filtration(
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

    auto [verts, filtration, mc_simplices, vai] = get_barycentric_subdivision_and_filtration(
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

    auto [verts_alpha, filt_alpha, mc_alpha, vai_alpha] = get_barycentric_subdivision_and_filtration(
        tet_points, colors_1111, radii, true, true);

    auto ca = count_filtration(filt_alpha);

    std::cout << "  Alpha:   " << verts_alpha.size() << " vertices, "
              << ca.vertices << "v " << ca.edges << "e " << ca.triangles << "t = "
              << ca.total << " simplices\n";

    // Compare with Delaunay (no alpha filtering)
    auto [verts_del, filt_del, mc_del, vai_del] = get_barycentric_subdivision_and_filtration(
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

    auto [verts, filtration, mc_simplices, vai] = get_barycentric_subdivision_and_filtration(
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

    auto [verts, filtration, mc_simplices, vai] = get_barycentric_subdivision_and_filtration(
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

    auto [verts_alpha, filt_alpha, mc_alpha, vai_alpha] = get_barycentric_subdivision_and_filtration(
        bipyramid_points, colors_bipyramid, radii, true, true);

    auto ca = count_filtration(filt_alpha);

    std::cout << "  Alpha:   " << verts_alpha.size() << " vertices, "
              << ca.vertices << "v " << ca.edges << "e " << ca.triangles << "t = "
              << ca.total << " simplices\n";

    auto [verts_del, filt_del, mc_del, vai_del] = get_barycentric_subdivision_and_filtration(
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

// ================= Bipyramid — 2-2 × 2 Shared Face ==================
// Same bipyramid, colors {1,2,2,1,1}: two 2-2 tets sharing multicolored face
// → double 2-2 counts minus 3 shared vertices and 2 shared edges

static const ColorLabels colors_bipyramid_22 = {1, 2, 2, 1, 1};

void test_bipyramid22_alpha_051() {
    std::cout << "Test: bipyramid 2-2x2 alpha=0.51 (edges only)\n";

    Radii radii(5, 0.51);

    auto [verts, filtration, mc_simplices, vai] = get_barycentric_subdivision_and_filtration(
        bipyramid_points, colors_bipyramid_22, radii, true, true);

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

void test_bipyramid22_alpha_058() {
    std::cout << "Test: bipyramid 2-2x2 alpha=0.58 (edges + triangles)\n";

    Radii radii(5, 0.58);

    auto [verts, filtration, mc_simplices, vai] = get_barycentric_subdivision_and_filtration(
        bipyramid_points, colors_bipyramid_22, radii, true, true);

    auto c = count_filtration(filtration);

    std::cout << "  " << verts.size() << " vertices, "
              << c.vertices << "v " << c.edges << "e " << c.triangles << "t = "
              << c.total << " simplices\n";

    assert(verts.size() == 13);
    assert(c.vertices == 13);
    assert(c.edges == 14);
    assert(c.triangles == 0);
    assert(c.total == 27);

    std::cout << "  PASS\n";
}

void test_bipyramid22_alpha_062() {
    std::cout << "Test: bipyramid 2-2x2 alpha=0.62 (full tets)\n";

    Radii radii(5, 0.62);

    auto [verts_alpha, filt_alpha, mc_alpha, vai_alpha] = get_barycentric_subdivision_and_filtration(
        bipyramid_points, colors_bipyramid_22, radii, true, true);

    auto ca = count_filtration(filt_alpha);

    std::cout << "  Alpha:   " << verts_alpha.size() << " vertices, "
              << ca.vertices << "v " << ca.edges << "e " << ca.triangles << "t = "
              << ca.total << " simplices\n";

    auto [verts_del, filt_del, mc_del, vai_del] = get_barycentric_subdivision_and_filtration(
        bipyramid_points, colors_bipyramid_22, radii, true, false);

    auto cd = count_filtration(filt_del);

    std::cout << "  Delaunay: " << verts_del.size() << " vertices, "
              << cd.vertices << "v " << cd.edges << "e " << cd.triangles << "t = "
              << cd.total << " simplices\n";

    assert(verts_alpha.size() == 15);
    assert(ca.vertices == 15);
    assert(ca.edges == 30);
    assert(ca.triangles == 16);
    assert(ca.total == 61);

    assert(verts_alpha.size() == verts_del.size());
    assert(ca.total == cd.total);

    std::cout << "  PASS\n";
}

// ================ Bipyramid — 1-1-1-1-1 Shared 1-1-1 Face ================
// All distinct colors: double 1-1-1-1 minus 4 shared vertices and 3 shared edges

static const ColorLabels colors_bipyramid_11111 = {1, 2, 3, 4, 5};

void test_bipyramid11111_alpha_051() {
    std::cout << "Test: bipyramid 1-1-1-1-1 alpha=0.51 (edges only)\n";

    Radii radii(5, 0.51);

    auto [verts, filtration, mc_simplices, vai] = get_barycentric_subdivision_and_filtration(
        bipyramid_points, colors_bipyramid_11111, radii, true, true);

    auto c = count_filtration(filtration);

    std::cout << "  " << verts.size() << " vertices, "
              << c.vertices << "v " << c.edges << "e " << c.triangles << "t = "
              << c.total << " simplices\n";

    assert(verts.size() == 9);
    assert(c.vertices == 9);
    assert(c.edges == 0);
    assert(c.triangles == 0);
    assert(c.total == 9);

    std::cout << "  PASS\n";
}

void test_bipyramid11111_alpha_058() {
    std::cout << "Test: bipyramid 1-1-1-1-1 alpha=0.58 (edges + triangles)\n";

    Radii radii(5, 0.58);

    auto [verts, filtration, mc_simplices, vai] = get_barycentric_subdivision_and_filtration(
        bipyramid_points, colors_bipyramid_11111, radii, true, true);

    auto c = count_filtration(filtration);

    std::cout << "  " << verts.size() << " vertices, "
              << c.vertices << "v " << c.edges << "e " << c.triangles << "t = "
              << c.total << " simplices\n";

    assert(verts.size() == 16);
    assert(c.vertices == 16);
    assert(c.edges == 21);
    assert(c.triangles == 0);
    assert(c.total == 37);

    std::cout << "  PASS\n";
}

void test_bipyramid11111_alpha_062() {
    std::cout << "Test: bipyramid 1-1-1-1-1 alpha=0.62 (full tets)\n";

    Radii radii(5, 0.62);

    auto [verts_alpha, filt_alpha, mc_alpha, vai_alpha] = get_barycentric_subdivision_and_filtration(
        bipyramid_points, colors_bipyramid_11111, radii, true, true);

    auto ca = count_filtration(filt_alpha);

    std::cout << "  Alpha:   " << verts_alpha.size() << " vertices, "
              << ca.vertices << "v " << ca.edges << "e " << ca.triangles << "t = "
              << ca.total << " simplices\n";

    auto [verts_del, filt_del, mc_del, vai_del] = get_barycentric_subdivision_and_filtration(
        bipyramid_points, colors_bipyramid_11111, radii, true, false);

    auto cd = count_filtration(filt_del);

    std::cout << "  Delaunay: " << verts_del.size() << " vertices, "
              << cd.vertices << "v " << cd.edges << "e " << cd.triangles << "t = "
              << cd.total << " simplices\n";

    assert(verts_alpha.size() == 18);
    assert(ca.vertices == 18);
    assert(ca.edges == 41);
    assert(ca.triangles == 24);
    assert(ca.total == 83);

    assert(verts_alpha.size() == verts_del.size());
    assert(ca.total == cd.total);

    std::cout << "  PASS\n";
}

// ============= Random Point Cloud Consistency Invariants =============

static std::string classify_tet_split(const Tetrahedron& tet, const ColorLabels& colors) {
    std::map<int, int> color_count;
    for (int v : tet) color_count[colors[v]]++;
    std::vector<int> sizes;
    for (const auto& [c, n] : color_count) sizes.push_back(n);
    std::sort(sizes.rbegin(), sizes.rend());
    if (sizes == std::vector<int>{3, 1}) return "3-1";
    if (sizes == std::vector<int>{2, 2}) return "2-2";
    if (sizes == std::vector<int>{2, 1, 1}) return "2-1-1";
    if (sizes == std::vector<int>{1, 1, 1, 1}) return "1-1-1-1";
    return "unknown";
}

static std::string classify_triangle_split(const std::vector<int>& tri, const ColorLabels& colors) {
    std::map<int, int> color_count;
    for (int v : tri) color_count[colors[v]]++;
    std::vector<int> sizes;
    for (const auto& [c, n] : color_count) sizes.push_back(n);
    std::sort(sizes.rbegin(), sizes.rend());
    if (sizes == std::vector<int>{2, 1}) return "2-1";
    if (sizes == std::vector<int>{1, 1, 1}) return "1-1-1";
    return "unknown";
}

struct FreeAnalysis {
    size_t free_edges;
    size_t isolated_vertices;
};

static FreeAnalysis analyze_free_simplices(const Filtration& f) {
    std::set<std::pair<int32_t, int32_t>> edges_in_triangles;
    std::set<std::pair<int32_t, int32_t>> all_edges;
    std::set<int32_t> all_vertices;
    std::set<int32_t> vertices_in_edges;

    for (const auto& [simplex, val] : f) {
        if (simplex.size() == 1) {
            all_vertices.insert(simplex[0]);
        } else if (simplex.size() == 2) {
            auto e = std::minmax(simplex[0], simplex[1]);
            all_edges.insert({e.first, e.second});
            vertices_in_edges.insert(simplex[0]);
            vertices_in_edges.insert(simplex[1]);
        } else if (simplex.size() == 3) {
            for (int i = 0; i < 3; i++) {
                for (int j = i + 1; j < 3; j++) {
                    auto e = std::minmax(simplex[i], simplex[j]);
                    edges_in_triangles.insert({e.first, e.second});
                }
                vertices_in_edges.insert(simplex[i]);
            }
        }
    }

    size_t free_edges = 0;
    for (const auto& e : all_edges) {
        if (edges_in_triangles.find(e) == edges_in_triangles.end())
            free_edges++;
    }

    size_t isolated = 0;
    for (int32_t v : all_vertices) {
        if (vertices_in_edges.find(v) == vertices_in_edges.end())
            isolated++;
    }

    return {free_edges, isolated};
}

static void dump_pointcloud(const std::string& dir, unsigned int seed,
                            const Points& points, const Radii& radii,
                            const ColorLabels& colors, const std::string& reason) {
    std::filesystem::create_directories(dir);
    std::string path = dir + "/seed_" + std::to_string(seed) + ".txt";
    std::ofstream out(path);
    out << "# Should this ever produce an erroneous example.\n";
    out << "# Please fix it :) ... or contact the developer.\n\n";
    out << "seed: " << seed << "\n";
    out << "failure: " << reason << "\n";
    out << "n_points: " << points.size() << "\n\n";
    out << "# index  x  y  z  radius  color\n";
    for (size_t i = 0; i < points.size(); i++) {
        out << i << "  "
            << points[i].x() << "  " << points[i].y() << "  " << points[i].z() << "  "
            << radii[i] << "  " << colors[i] << "\n";
    }
    out.close();
    std::cerr << "  Dumped failing point cloud to " << path << "\n";
}

static const std::string fail_dir =
    std::filesystem::path(__FILE__).parent_path().string() + "/failed-random-generated-tests";

bool test_random_pointcloud_consistency(unsigned int seed = 42) {
    std::cout << "Test: random 3-color point cloud consistency (seed=" << seed << ")\n";

    const int n_points = 25;
    const int n_colors = 3;
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> pos_dist(0.0, 1.0);
    std::uniform_real_distribution<double> rad_dist(0.25, 0.35);
    std::uniform_int_distribution<int> color_dist(1, n_colors);

    Points points(n_points);
    Radii radii(n_points);
    ColorLabels colors(n_points);

    for (int i = 0; i < n_points; i++) {
        points[i] = Point3D(pos_dist(rng), pos_dist(rng), pos_dist(rng));
        radii[i] = rad_dist(rng);
        colors[i] = color_dist(rng);
    }

    auto [verts, filtration, mc_simplices, vai] = get_barycentric_subdivision_and_filtration(
        points, colors, radii, true, true);

    auto c = count_filtration(filtration);

    std::cout << "  " << verts.size() << " vertices, "
              << c.vertices << "v " << c.edges << "e " << c.triangles << "t = "
              << c.total << " simplices\n";
    std::cout << "  generating tets: " << mc_simplices.generating_tetrahedra.size()
              << ", free tris: " << mc_simplices.generating_free_triangles.size()
              << ", free edges: " << mc_simplices.generating_free_edges.size() << "\n";

    // Invariant 1: Triangle count from generating tetrahedra
    size_t n_31 = 0, n_22 = 0, n_211 = 0, n_1111 = 0;
    for (const auto& tet : mc_simplices.generating_tetrahedra) {
        auto split = classify_tet_split(tet, colors);
        if (split == "3-1") n_31++;
        else if (split == "2-2") n_22++;
        else if (split == "2-1-1") n_211++;
        else if (split == "1-1-1-1") n_1111++;
        else {
            dump_pointcloud(fail_dir, seed, points, radii, colors,
                "unexpected tet split type: " + split);
            return false;
        }
    }

    size_t expected_triangles = n_31 * 6 + n_22 * 8 + n_211 * 10 + n_1111 * 12;
    std::cout << "  tet splits: " << n_31 << " x 3-1, " << n_22 << " x 2-2, "
              << n_211 << " x 2-1-1, " << n_1111 << " x 1-1-1-1\n";
    std::cout << "  expected triangles: " << expected_triangles
              << ", actual: " << c.triangles << "\n";
    if (c.triangles != expected_triangles) {
        dump_pointcloud(fail_dir, seed, points, radii, colors,
            "triangle count mismatch: expected " + std::to_string(expected_triangles)
            + ", got " + std::to_string(c.triangles));
        return false;
    }

    // Invariant 2: Free filtration edges from free generating triangles
    auto fa = analyze_free_simplices(filtration);

    size_t expected_free_edges = 0;
    for (const auto& tri : mc_simplices.generating_free_triangles) {
        auto split = classify_triangle_split(tri, colors);
        if (split == "2-1") expected_free_edges += 2;
        else if (split == "1-1-1") expected_free_edges += 3;
        else {
            dump_pointcloud(fail_dir, seed, points, radii, colors,
                "unexpected triangle split type: " + split);
            return false;
        }
    }

    std::cout << "  expected free edges: " << expected_free_edges
              << ", actual: " << fa.free_edges << "\n";
    if (fa.free_edges != expected_free_edges) {
        dump_pointcloud(fail_dir, seed, points, radii, colors,
            "free edge count mismatch: expected " + std::to_string(expected_free_edges)
            + ", got " + std::to_string(fa.free_edges));
        return false;
    }

    // Invariant 3: Isolated vertices from free generating edges
    size_t expected_isolated = mc_simplices.generating_free_edges.size();
    std::cout << "  expected isolated vertices: " << expected_isolated
              << ", actual: " << fa.isolated_vertices << "\n";
    if (fa.isolated_vertices != expected_isolated) {
        dump_pointcloud(fail_dir, seed, points, radii, colors,
            "isolated vertex count mismatch: expected " + std::to_string(expected_isolated)
            + ", got " + std::to_string(fa.isolated_vertices));
        return false;
    }

    std::cout << "  PASS\n";
    return true;
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

        std::cout << "\n--- Bipyramid 2-2x2 ---\n";
        test_bipyramid22_alpha_051();
        test_bipyramid22_alpha_058();
        test_bipyramid22_alpha_062();

        std::cout << "\n--- Bipyramid 1-1-1-1-1 ---\n";
        test_bipyramid11111_alpha_051();
        test_bipyramid11111_alpha_058();
        test_bipyramid11111_alpha_062();

        std::cout << "\n--- Random Point Cloud (100 random seeds) ---\n";
        std::random_device rd;
        std::mt19937 seed_gen(rd());
        int n_failed = 0;
        for (int i = 0; i < 100; i++) {
            unsigned int seed = seed_gen();
            if (!test_random_pointcloud_consistency(seed))
                n_failed++;
        }
        if (n_failed > 0) {
            std::cerr << "\n" << n_failed << " random seeds failed! "
                      << "See " << fail_dir << "/\n";
            return 1;
        }

        std::cout << "\nAll alpha construction tests passed!\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nTest failed with exception: " << e.what() << "\n";
        return 1;
    }
}
