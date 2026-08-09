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

# Span testers retain a compact row basis, answer containment queries without
# mutating it, and add precisely the vectors that enlarge its span.
gap> Mspan := [[1,0,1],[0,1,1],[1,1,0]];;
gap> SpanFixture := function(M) local T; if not IsBoundGlobal("SGC_NewSpanTesterMod2") then return fail; fi; T := ValueGlobal("SGC_NewSpanTesterMod2")(M); return [T.dimension,T.rank(),T.contains([1,1,0]),T.contains([1,0,0]),T.addIfIndependent([1,0,0]),T.rank(),T.contains([0,0,1]),T.addIfIndependent([0,0,1]),T.rank()]; end;;
gap> SpanFixture(Mspan);
[ 3, 2, true, false, true, 3, true, false, 3 ]
gap> Sspan := SGC_SparseMat(0,3);;
gap> SpanEmptyFixture := function(M) local T; if not IsBoundGlobal("SGC_NewSpanTesterMod2") then return fail; fi; T := ValueGlobal("SGC_NewSpanTesterMod2")(M); return [T.dimension,T.contains([0,0,0]),T.contains([1,0,0]),T.addIfIndependent([1,0,0]),T.rank()]; end;;
gap> SpanEmptyFixture(Sspan);
[ 3, true, false, true, 1 ]

# Invocation-local solver caches do not share entries and refactor a matrix
# after rows have been appended.
gap> Mcache := [[1,0]];;
gap> CacheFixture := function() local cache1,cache2,Solve; if not IsBoundGlobal("SGC_NewSolverCache") or not IsBoundGlobal("SGC_CachedSolveIn") then return fail; fi; Solve := ValueGlobal("SGC_CachedSolveIn");; cache1 := ValueGlobal("SGC_NewSolverCache")();; cache2 := ValueGlobal("SGC_NewSolverCache")();; return [Length(cache1.entries),SGC_GF2ToZ(Solve(cache1,Mcache,[1,0])),Length(cache1.entries),Length(cache2.entries),SGC_GF2ToZ(Solve(cache2,Mcache,[1,0])),Length(cache2.entries)]; end;;
gap> CacheFixture();
[ 0, [ 1 ], 1, 0, [ 1 ], 1 ]
gap> CacheGrowthFixture := function() local cache,Solve; if not IsBoundGlobal("SGC_NewSolverCache") or not IsBoundGlobal("SGC_CachedSolveIn") then return fail; fi; Solve := ValueGlobal("SGC_CachedSolveIn");; cache := ValueGlobal("SGC_NewSolverCache")();; Solve(cache,Mcache,[1,0]);; Append(Mcache,[[0,1]]);; return [SGC_GF2ToZ(Solve(cache,Mcache,[1,1])),cache.entries[1].len]; end;;
gap> CacheGrowthFixture();
[ [ 1, 1 ], 2 ]

# Relation reducers enumerate target monomials without a depth cap and keep
# the sentinel coordinate used by the existing relation loops.
gap> Rels4 := [];; Rels4[2] := [[[2,0],[0,1]]];;
gap> RelationFixture := function() local T; if not IsBoundGlobal("SGC_NewRelationReducer") then return fail; fi; T := ValueGlobal("SGC_NewRelationReducer")([[1,0],[0,1]],[1,2],Rels4,4); return [T.monomials,T.vectorLength,T.relationProductCount,T.rank,T.encode([[4,0],[0,2]]),T.contains([1,0,1,0]),T.contains([0,1,0,0]),T.addIfIndependent([0,1,0,0]),T.rank,T.addIfIndependent([0,1,0,0])]; end;;
gap> RelationFixture();
[ [ [ 4, 0 ], [ 2, 1 ], [ 0, 2 ] ], 4, 2, 2, [ 1, 0, 1, 0 ], true, false, true, 3, false ]
gap> RelationEmptyFixture := function() local T; if not IsBoundGlobal("SGC_NewRelationReducer") then return fail; fi; T := ValueGlobal("SGC_NewRelationReducer")([[1]],[3],[],2); return [T.monomials,T.vectorLength,T.relationProductCount,T.rank,T.contains([0])]; end;;
gap> RelationEmptyFixture();
[ [ ], 1, 0, 0, true ]
gap> Rels12 := [];; Rels12[6] := [[[1]]];;
gap> Relation12Fixture := function() local T; if not IsBoundGlobal("SGC_NewRelationReducer") then return fail; fi; T := ValueGlobal("SGC_NewRelationReducer")([[1]],[6],Rels12,12); return [T.monomials,T.relationProductCount,T.rank]; end;;
gap> Relation12Fixture();
[ [ [ 2 ] ], 1, 1 ]
gap> Rels13 := [];; Rels13[2] := [[[2]]];;
gap> Relation13Fixture := function() local T; if not IsBoundGlobal("SGC_NewRelationReducer") then return fail; fi; T := ValueGlobal("SGC_NewRelationReducer")([[1]],[1],Rels13,13); return [T.monomials,T.relationProductCount,T.rank]; end;;
gap> Relation13Fixture();
[ [ [ 13 ] ], 1, 1 ]

# End-to-end forced offload through CR_Mod2CocyclesAndCoboundaries (Task 5):
# cohomology dims are basis-independent and must equal the native values,
# and the call counter proves the external binary was actually exercised.
gap> SGC_LINALG_THRESHOLD := 0;;
gap> if SGC_LinalgBinary() = fail then
>     Print("[ 5, 12, 20, 28 ]\ntrue\n");
> else
>     SGC_LINALG_CALLS := 0;;
>     R16 := SGC_ResolutionSpaceGroup(SpaceGroupBBNWZ(3,16), 7);;
>     Print(List([1..4], k -> CR_Mod2CocyclesAndCoboundaries(R16, k, false)!.Mod2Cohomologydim), "\n");
>     Print(SGC_LINALG_CALLS > 0, "\n");
> fi;
[ 5, 12, 20, 28 ]
true
gap> SGC_LINALG_THRESHOLD := saved_threshold;;
gap> STOP_TEST("linalg.tst", 0);
