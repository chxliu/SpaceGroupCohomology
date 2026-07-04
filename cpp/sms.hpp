#pragma once
// SMS sparse-triplet text format.
// Header: "nrows ncols mod"  (mod 2 = GF(2); mod 0 reserved for integers).
// Body:   1-based "i j v" triplets, duplicates accumulate; terminator "0 0 0".
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

struct Triplet { long i, j; long long v; };

struct SmsMatrix {
    long nrows = 0, ncols = 0;
    long mod   = 0;
    std::vector<Triplet> entries;
};

inline SmsMatrix read_sms(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "r");
    if (!f) throw std::runtime_error("cannot open " + path);
    SmsMatrix m;
    if (std::fscanf(f, "%ld %ld %ld", &m.nrows, &m.ncols, &m.mod) != 3) {
        std::fclose(f);
        throw std::runtime_error("bad SMS header in " + path);
    }
    if (m.nrows < 0 || m.ncols < 0)
        throw std::runtime_error("negative dimensions in " + path);
    long i, j;
    long long v;
    while (std::fscanf(f, "%ld %ld %lld", &i, &j, &v) == 3) {
        if (i == 0 && j == 0) break;
        if (i < 1 || i > m.nrows || j < 1 || j > m.ncols) {
            std::fclose(f);
            throw std::runtime_error("entry out of range in " + path);
        }
        m.entries.push_back({i, j, v});
    }
    std::fclose(f);
    return m;
}

inline void write_sms(const std::string& path, long nrows, long ncols, long mod,
                      const std::vector<Triplet>& entries) {
    FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + path);
    std::fprintf(f, "%ld %ld %ld\n", nrows, ncols, mod);
    for (const auto& t : entries)
        std::fprintf(f, "%ld %ld %lld\n", t.i, t.j, t.v);
    std::fprintf(f, "0 0 0\n");
    std::fclose(f);
}
