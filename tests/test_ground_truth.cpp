// Validates single-tetrahedron subdivisions (2-2, 3-1, 2-1-1, 1-1-1-1
// colorings, plus a weighted 2-2 case) against reference output from the
// original Julia implementation in tests/data/ground_truth.json.

#include "ground_truth_utils.hpp"

int main() {
    bool all_passed = ground_truth::run_file(std::string(TEST_DATA_DIR) + "/ground_truth.json");

    if (!all_passed) {
        std::cout << "Some ground truth tests FAILED!\n";
        return 1;
    }
    std::cout << "All ground truth tests passed!\n";
    return 0;
}
