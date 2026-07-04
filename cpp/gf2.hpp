#pragma once
// Bit-packed GF(2) matrices and deterministic semi-echelon elimination.
// Row convention matches GAP: solve2 finds x with x*A = b; nullsp2 is the
// left nullspace {x : x*A = 0}. Results are deterministic but NOT required
// to match GAP's basis choices (relaxed-parity contract; see the plan).
#include "sms.hpp"
#include <cstdint>
#include <utility>
#include <vector>

struct GF2Row {
    std::vector<uint64_t> w;
    explicit GF2Row(long ncols = 0) : w((ncols + 63) / 64, 0) {}
    bool get(long j) const { return (w[j >> 6] >> (j & 63)) & 1; }
    void set(long j) { w[j >> 6] |= (uint64_t(1) << (j & 63)); }
    void flip(long j) { w[j >> 6] ^= (uint64_t(1) << (j & 63)); }
    void add(const GF2Row& o) {
        for (size_t k = 0; k < w.size(); ++k) w[k] ^= o.w[k];
    }
    bool zero() const {
        for (uint64_t x : w) if (x) return false;
        return true;
    }
    long first_set(long ncols) const {  // -1 if zero
        for (size_t k = 0; k < w.size(); ++k)
            if (w[k]) {
                long j = long(k) * 64 + __builtin_ctzll(w[k]);
                return j < ncols ? j : -1;
            }
        return -1;
    }
};

struct GF2Mat {
    long nrows = 0, ncols = 0;
    std::vector<GF2Row> rows;
    static GF2Mat from_sms(const SmsMatrix& m) {
        GF2Mat a;
        a.nrows = m.nrows;
        a.ncols = m.ncols;
        a.rows.assign(m.nrows, GF2Row(m.ncols));
        for (const auto& t : m.entries)
            if (t.v % 2 != 0) a.rows[t.i - 1].flip(t.j - 1);  // accumulate mod 2
        return a;
    }
};

// Semi-echelon data of A (rows processed in order, columns swept ascending):
//   vectors[k]  : semi-echelon basis rows of the row space
//   coeffs[k]   : coeffs[k] * A = vectors[k]
//   heads[j]    : 1-based index into vectors of the row with pivot j, or 0
//   relations   : basis of the left nullspace {x : x*A = 0}
struct SemiEch {
    std::vector<GF2Row> vectors, coeffs, relations;
    std::vector<long> heads;
};

inline SemiEch semi_echelon_transformation(const GF2Mat& A) {
    SemiEch se;
    se.heads.assign(A.ncols, 0);
    for (long i = 0; i < A.nrows; ++i) {
        GF2Row row = A.rows[i];
        GF2Row coeff(A.nrows);
        coeff.set(i);
        for (long j = 0; j < A.ncols; ++j)
            if (row.get(j) && se.heads[j] != 0) {
                row.add(se.vectors[se.heads[j] - 1]);
                coeff.add(se.coeffs[se.heads[j] - 1]);
            }
        long piv = row.first_set(A.ncols);
        if (piv >= 0) {
            se.vectors.push_back(row);
            se.coeffs.push_back(coeff);
            se.heads[piv] = long(se.vectors.size());
        } else {
            se.relations.push_back(coeff);
        }
    }
    return se;
}

// Solve x*A = b through precomputed semi-echelon data; {found, x}.
inline std::pair<bool, GF2Row> solution_mat(const SemiEch& se, long nrowsA, GF2Row b) {
    GF2Row x(nrowsA);
    for (long j = 0; j < long(se.heads.size()); ++j)
        if (b.get(j)) {
            if (se.heads[j] == 0) return {false, GF2Row(nrowsA)};
            b.add(se.vectors[se.heads[j] - 1]);
            x.add(se.coeffs[se.heads[j] - 1]);
        }
    return {true, x};
}

// Canonical RREF row basis, ordered by pivot column ascending.
inline std::vector<GF2Row> base_mat(const GF2Mat& A) {
    SemiEch se = semi_echelon_transformation(A);
    std::vector<GF2Row> out;
    std::vector<long> pivs;
    for (long j = 0; j < A.ncols; ++j)
        if (se.heads[j] != 0) {
            out.push_back(se.vectors[se.heads[j] - 1]);
            pivs.push_back(j);
        }
    for (size_t r = 0; r < out.size(); ++r)
        for (size_t s = 0; s < out.size(); ++s)
            if (s != r && out[s].get(pivs[r])) out[s].add(out[r]);
    return out;
}

// Deterministic complement basis of rowspace(I) inside rowspace(K):
// echelonize I, then sift rows of K in order through I's heads and through
// already-accepted complement rows; nonzero remainders join the complement.
// Same role as GAP BaseSteinitzVectors(K, I).factorspace (span-level).
inline std::vector<GF2Row> base_steinitz_factorspace(const GF2Mat& K, const GF2Mat& I) {
    long ncols = K.ncols;
    std::vector<GF2Row> ibasis, comp;
    std::vector<long> iheads(ncols, 0), cheads(ncols, 0);
    for (long i = 0; i < I.nrows; ++i) {
        GF2Row row = I.rows[i];
        for (long j = 0; j < ncols; ++j)
            if (row.get(j) && iheads[j] != 0) row.add(ibasis[iheads[j] - 1]);
        long piv = row.first_set(ncols);
        if (piv >= 0) {
            ibasis.push_back(row);
            iheads[piv] = long(ibasis.size());
        }
    }
    for (long i = 0; i < K.nrows; ++i) {
        GF2Row row = K.rows[i];
        for (long j = 0; j < ncols; ++j) {
            if (!row.get(j)) continue;
            if (iheads[j] != 0)      row.add(ibasis[iheads[j] - 1]);
            else if (cheads[j] != 0) row.add(comp[cheads[j] - 1]);
        }
        long piv = row.first_set(ncols);
        if (piv >= 0) {
            comp.push_back(row);
            cheads[piv] = long(comp.size());
        }
    }
    return comp;
}
