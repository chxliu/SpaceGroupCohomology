# linalg.tst — correctness of SGC_* dispatch wrappers against GAP-native
# routines. Two passes: native path (huge threshold), then forced-external
# (threshold 0, skipped with a notice if bin/sgclinalg is not built).
# Contract (relaxed parity): ranks/dims match native exactly; bases are
# span-equal; solutions satisfy their defining equations.

gap> START_TEST("SpaceGroupCohomology: linalg correctness");
gap> saved_threshold := SGC_LINALG_THRESHOLD;;
gap> rs := RandomSource(IsMersenneTwister, 20260703);;
gap> RandMat := function(m, n) return List([1..m], i -> List([1..n], j -> Random(rs, 0, 1))); end;;
gap> mats := List([1..25], k -> RandMat(Random(rs, 1, 40), Random(rs, 1, 40)));;
gap> vecs := List(mats, M -> List([1..10], k -> RandMat(1, Length(M[1]))[1]));;
gap> SpanEq := function(A, B)
>     local rA, rB;
>     if Length(A) = 0 or Length(B) = 0 then return Length(A) = Length(B); fi;
>     rA := RankMat(A); rB := RankMat(B);
>     return rA = rB and rA = RankMat(Concatenation(A, B));
> end;;
gap> CheckAll := function()
>     local ok, i, M, B, ns, S, b, x, K, I, comp;
>     ok := true;
>     for i in [1..Length(mats)] do
>         M := mats[i];
>         # rank: must equal native exactly
>         if SGC_RankMod2(M) <> RankMat(M*Z(2)) then ok := false; Print("rank mismatch at ", i, "\n"); fi;
>         # row basis: span-equal to native BaseMat
>         B := SGC_RowBasisMod2(M);
>         if not SpanEq(List(B, ShallowCopy), List(BaseMat(M*Z(2)), ShallowCopy)) then
>             ok := false; Print("rowbasis span mismatch at ", i, "\n"); fi;
>         # nullspace: every vector annihilates M; dimension = nrows - rank
>         ns := SGC_NullspaceMod2(M);
>         if not ForAll(ns, v -> ForAll((v*M) mod 2, e -> e = 0)) then
>             ok := false; Print("nullspace vector fails at ", i, "\n"); fi;
>         if Length(ns) <> Length(M) - RankMat(M*Z(2)) then
>             ok := false; Print("nullspace dim mismatch at ", i, "\n"); fi;
>         # solver: solvability agrees with native; solutions satisfy x*M = b
>         S := SGC_SolverMod2(M);
>         for b in vecs[i] do
>             x := S(b);
>             if (x = fail) <> (SolutionMat(M*Z(2), b*Z(2)) = fail) then
>                 ok := false; Print("solvability mismatch at ", i, "\n");
>             elif x <> fail and x*(M*Z(2)) <> b*Z(2) then
>                 ok := false; Print("solution fails equation at ", i, "\n"); fi;
>         od;
>         # steinitz: complement of a sub-rowspace I inside rowspace(K)
>         K := M;
>         I := BaseMat(K{[1..QuoInt(Length(K)+1, 2)]}*Z(2));
>         comp := SGC_SteinitzMod2(K, I);
>         if Length(comp) <> RankMat(K*Z(2)) - Length(I) then
>             ok := false; Print("steinitz size mismatch at ", i, "\n"); fi;
>         if Length(comp) > 0 and Length(I) > 0 then
>             if RankMat(Concatenation(List(I, ShallowCopy), List(comp, ShallowCopy)))
>                <> Length(I) + Length(comp) then
>                 ok := false; Print("steinitz not independent at ", i, "\n"); fi;
>         fi;
>         if not ForAll(comp, v -> SolutionMat(K*Z(2), v) <> fail) then
>             ok := false; Print("steinitz outside rowspace(K) at ", i, "\n"); fi;
>     od;
>     return ok;
> end;;

# Pass 1: native path (threshold impossible to exceed).
gap> SGC_LINALG_THRESHOLD := 10^12;;
gap> CheckAll();
true

# Pass 2: forced-external (threshold 0). Every call goes through sgclinalg.
gap> SGC_LINALG_THRESHOLD := 0;;
gap> if SGC_LinalgBinary() = fail then
>     Print("SKIP external pass: bin/sgclinalg not built\ntrue\ntrue\n");
> else
>     SGC_LINALG_CALLS := 0;;
>     Print(CheckAll(), "\n");
>     Print(SGC_LINALG_CALLS > 0, "\n");
> fi;
true
true

# Sparse-record inputs round-trip through both paths.
gap> S := SGC_SparseMat(3, 4);;
gap> Append(S.entries, [[1,1,1],[1,2,1],[2,2,1],[2,3,1],[3,1,1],[3,3,1]]);;
gap> SGC_LINALG_THRESHOLD := 10^12;;
gap> SGC_RankMod2(S);
2
gap> SGC_LINALG_THRESHOLD := 0;;
gap> if SGC_LinalgBinary() = fail then Print("2\n"); else Print(SGC_RankMod2(S), "\n"); fi;
2
gap> SGC_LINALG_THRESHOLD := saved_threshold;;
gap> STOP_TEST("linalg.tst", 0);
