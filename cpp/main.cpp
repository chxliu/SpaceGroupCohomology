// sgclinalg — external sparse GF(2) linear algebra for SpaceGroupCohomology.
// Usage: sgclinalg <op> <in.sms> [in2.sms] <out>
// Ops:   rank2, rowbasis2, nullsp2, solve2 (A B out), semiech2, steinitz2 (K I out)
// All I/O in SMS sparse-triplet text format (see sms.hpp). Exit 0 on success.
#include "gf2.hpp"
#include "sms.hpp"
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

static std::vector<Triplet> rows_to_triplets(const std::vector<GF2Row>& rows, long ncols) {
    std::vector<Triplet> out;
    for (size_t r = 0; r < rows.size(); ++r)
        for (long j = 0; j < ncols; ++j)
            if (rows[r].get(j)) out.push_back({long(r) + 1, j + 1, 1});
    return out;
}

int main(int argc, char** argv) {
    try {
        if (argc < 4) {
            std::fprintf(stderr, "usage: sgclinalg <op> <in.sms> [in2.sms] <out>\n");
            return 2;
        }
        std::string op = argv[1];

        if (op == "rank2") {
            GF2Mat A = GF2Mat::from_sms(read_sms(argv[2]));
            SemiEch se = semi_echelon_transformation(A);
            FILE* f = std::fopen(argv[3], "w");
            if (!f) { std::fprintf(stderr, "cannot write %s\n", argv[3]); return 1; }
            std::fprintf(f, "rank %zu\n", se.vectors.size());
            std::fclose(f);

        } else if (op == "rowbasis2") {
            GF2Mat A = GF2Mat::from_sms(read_sms(argv[2]));
            auto basis = base_mat(A);
            write_sms(argv[3], long(basis.size()), A.ncols, 2,
                      rows_to_triplets(basis, A.ncols));

        } else if (op == "nullsp2") {
            GF2Mat A = GF2Mat::from_sms(read_sms(argv[2]));
            SemiEch se = semi_echelon_transformation(A);
            write_sms(argv[3], long(se.relations.size()), A.nrows, 2,
                      rows_to_triplets(se.relations, A.nrows));

        } else if (op == "solve2") {
            if (argc < 5) { std::fprintf(stderr, "solve2 needs A B out\n"); return 2; }
            GF2Mat A = GF2Mat::from_sms(read_sms(argv[2]));
            GF2Mat B = GF2Mat::from_sms(read_sms(argv[3]));
            if (B.ncols != A.ncols) {
                std::fprintf(stderr, "solve2: B ncols %ld != A ncols %ld\n", B.ncols, A.ncols);
                return 1;
            }
            SemiEch se = semi_echelon_transformation(A);
            FILE* f = std::fopen(argv[4], "w");
            if (!f) { std::fprintf(stderr, "cannot write %s\n", argv[4]); return 1; }
            for (long k = 0; k < B.nrows; ++k) {
                auto [ok, x] = solution_mat(se, A.nrows, B.rows[k]);
                if (!ok) {
                    std::fprintf(f, "fail\n");
                    continue;
                }
                std::fprintf(f, "sol ");
                for (long i = 0; i < A.nrows; ++i) std::fputc(x.get(i) ? '1' : '0', f);
                std::fputc('\n', f);
            }
            std::fclose(f);

        } else if (op == "semiech2") {
            GF2Mat A = GF2Mat::from_sms(read_sms(argv[2]));
            SemiEch se = semi_echelon_transformation(A);
            FILE* f = std::fopen(argv[3], "w");
            if (!f) { std::fprintf(stderr, "cannot write %s\n", argv[3]); return 1; }
            std::fprintf(f, "nvec %zu ncols %ld nrowsA %ld\n",
                         se.vectors.size(), A.ncols, A.nrows);
            for (const auto& v : se.vectors) {
                std::fprintf(f, "vec ");
                for (long j = 0; j < A.ncols; ++j) std::fputc(v.get(j) ? '1' : '0', f);
                std::fputc('\n', f);
            }
            for (const auto& c : se.coeffs) {
                std::fprintf(f, "coeff ");
                for (long i = 0; i < A.nrows; ++i) std::fputc(c.get(i) ? '1' : '0', f);
                std::fputc('\n', f);
            }
            std::fprintf(f, "heads");
            for (long j = 0; j < A.ncols; ++j) std::fprintf(f, " %ld", se.heads[j]);
            std::fputc('\n', f);
            std::fclose(f);

        } else if (op == "steinitz2") {
            if (argc < 5) { std::fprintf(stderr, "steinitz2 needs K I out\n"); return 2; }
            GF2Mat K = GF2Mat::from_sms(read_sms(argv[2]));
            GF2Mat I = GF2Mat::from_sms(read_sms(argv[3]));
            if (K.ncols != I.ncols) {
                std::fprintf(stderr, "steinitz2: ncols mismatch\n");
                return 1;
            }
            auto comp = base_steinitz_factorspace(K, I);
            write_sms(argv[4], long(comp.size()), K.ncols, 2,
                      rows_to_triplets(comp, K.ncols));

        } else {
            std::fprintf(stderr, "unknown op %s\n", op.c_str());
            return 2;
        }
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "sgclinalg: %s\n", e.what());
        return 1;
    }
}
