// Validates the C++ implementation on large inputs against reference output
// from the original Julia implementation (DelaunayInterfaces.jl, via
// generate_ground_truth.jl):
//   - 4bmg protein dimer, weighted Delaunay (2259 atoms)
//   - 14 random point clouds: 7 unweighted Delaunay + 7 weighted Delaunay
// Expected values are reference output, not hand-derived; see
// ALPHA_CONSTRUCTION_TESTS.md for the hand-derived cases.

#include "ground_truth_utils.hpp"

int main() {
    bool all_passed = true;

    if (!ground_truth::run_single_case_file(
            "4bmg_dimer_delaunay",
            std::string(TEST_DATA_DIR) + "/ground_truth_4bmg_dimer_delaunay.json")) {
        all_passed = false;
    }

    if (!ground_truth::run_file(std::string(TEST_DATA_DIR) + "/ground_truth_random.json")) {
        all_passed = false;
    }

    if (!all_passed) {
        std::cout << "Some large example tests FAILED!\n";
        return 1;
    }
    std::cout << "All large example tests passed!\n";
    return 0;
}
