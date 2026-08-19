// Internal header: bucket sort used by BarycentricSubdivision's filtration
// finalization, extracted so benchmarks/benchmark_radix_sort.cpp can compare
// it against std::sort on identical data. Not part of the public API.
#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

namespace delaunay_interfaces::detail {

// Sorts a filtration bucket by (value, ids lexicographic), matching
// std::sort with the same comparator exactly. Filtration values are mean
// pairwise distances: always >= 0 with a clear sign bit and never NaN, so
// the raw 64-bit IEEE-754 pattern ordered as an unsigned integer is exactly
// the numeric order, and a stable LSD radix over its 8 bytes sorts by value.
// The radix leaves equal-value entries in emission order, not id order, so
// each equal-value run is re-sorted by ids afterwards (equal values here
// means equal bit patterns — duplicates derive from identical arithmetic —
// and the runs are typically tiny). A pass whose byte is shared by every key
// is the identity permutation and is skipped; small buckets fall back to
// std::sort, which is faster there and keeps tiny cases fast.
template <typename Entry, typename FullLess, typename IdLess>
void sort_bucket(std::vector<Entry>& entries, FullLess full_less, IdLess id_less) {
    constexpr size_t kRadixMinSize = 10000;
    const size_t n = entries.size();
    if (n < kRadixMinSize) {
        std::sort(entries.begin(), entries.end(), full_less);
        return;
    }

#ifndef NDEBUG
    for (const Entry& e : entries) {
        assert(e.value >= 0.0 && !std::signbit(e.value));
    }
#endif

    auto key_bits = [](const Entry& e) {
        uint64_t bits;
        std::memcpy(&bits, &e.value, sizeof bits);
        return bits;
    };

    std::vector<Entry> scratch(n);
    Entry* src = entries.data();
    Entry* dst = scratch.data();
    for (int shift = 0; shift < 64; shift += 8) {
        size_t count[256] = {};
        for (size_t i = 0; i < n; ++i) {
            ++count[(key_bits(src[i]) >> shift) & 0xff];
        }
        if (std::find(std::begin(count), std::end(count), n) != std::end(count)) {
            continue;
        }
        size_t offset = 0;
        for (size_t& c : count) {
            const size_t bucket = c;
            c = offset;
            offset += bucket;
        }
        for (size_t i = 0; i < n; ++i) {
            dst[count[(key_bits(src[i]) >> shift) & 0xff]++] = src[i];
        }
        std::swap(src, dst);
    }
    if (src != entries.data()) {
        std::memcpy(entries.data(), src, n * sizeof(Entry));
    }

    size_t run_begin = 0;
    for (size_t i = 1; i <= n; ++i) {
        if (i == n || entries[i].value != entries[run_begin].value) {
            if (i - run_begin > 1) {
                std::sort(entries.begin() + run_begin, entries.begin() + i, id_less);
            }
            run_begin = i;
        }
    }
}

} // namespace delaunay_interfaces::detail
