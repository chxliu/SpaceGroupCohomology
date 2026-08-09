# relations.tst — integration coverage for degree-generic relation reduction.

gap> START_TEST("SpaceGroupCohomology: relation reduction");
gap> SetRecursionTrapInterval(100000000);;

# Every relation-search degree must construct the generic reducer.  The
# observer is intentionally optional in production and records no output by
# itself; tests use it to characterize the integration path.
gap> relationReducerTrace := [];;
gap> relationReducerStats := [];;
gap> relationObserverSawMutableReducer := false;;
gap> relationObserverWasBound := IsBoundGlobal("SGC_RELATION_REDUCER_OBSERVER");;
gap> relationObserverWasReadOnly := relationObserverWasBound and IsReadOnlyGlobal("SGC_RELATION_REDUCER_OBSERVER");;
gap> if relationObserverWasBound then relationSavedObserver := ValueGlobal("SGC_RELATION_REDUCER_OBSERVER"); if relationObserverWasReadOnly then MakeReadWriteGlobal("SGC_RELATION_REDUCER_OBSERVER"); fi; fi;;
gap> SGC_RELATION_REDUCER_OBSERVER := function(reducer)
>     Add(relationReducerTrace,reducer.degree);
>     Add(relationReducerStats,
>         [reducer.degree,reducer.relationProductCount,reducer.rank]);
>     if IsBound(reducer.addIfIndependent) then
>         relationObserverSawMutableReducer := true;
>     fi;
> end;;

gap> degreeFamilyView := function(dims,degree)
>     return List(SGC_DegreeCupFamilies(dims,degree),
>         family->[family.factorDegree,family.sourceDegree,
>                  family.restrictWidth,family.candidateOrder]);
> end;;
gap> [degreeFamilyView([3,2,1,1,0,2],2),
>     degreeFamilyView([3,2,1,1,0,2],7),
>     degreeFamilyView([3,2,1,1,0,2],12)];
[ [ [ 1, 1, 0, "offDiagonalThenSquares" ] ],
  [ [ 1, 6, 0, "cartesian" ], [ 2, 5, 3, "cartesian" ],
    [ 3, 4, 5, "cartesian" ] ],
  [ [ 1, 11, 0, "cartesian" ], [ 2, 10, 3, "cartesian" ],
    [ 3, 9, 5, "cartesian" ], [ 4, 8, 6, "cartesian" ],
    [ 6, 6, 7, "unorderedPairs" ] ] ]

# IT219 exercises every relation degree through 12.  Coefficient solves in
# this invocation must not populate the compatibility process-global cache.
gap> savedSolverCacheList := SGC_SolverCacheList;;
gap> SGC_SolverCacheList := [];;
gap> D219relations := SGC_CohomologyData(219);;
gap> relationGlobalCacheGrowth := Length(SGC_SolverCacheList);;
gap> [relationReducerTrace,relationObserverSawMutableReducer];
[ [ 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ], false ]
gap> relationGlobalCacheGrowth;
0

# A supplied degree-6 generator list is authoritative for a named IT group;
# the compatibility fallback must not recompute Mod2RingGenerators(R,6).
gap> originalRingGenerators := Mod2RingGenerators;;
gap> ringGeneratorsWasReadOnly := IsReadOnlyGlobal("Mod2RingGenerators");;
gap> if ringGeneratorsWasReadOnly then MakeReadWriteGlobal("Mod2RingGenerators"); fi;;
gap> suppliedGen6FallbackCalls := 0;;
gap> Mod2RingGenerators := function(arg) suppliedGen6FallbackCalls := suppliedGen6FallbackCalls+1; return CallFuncList(originalRingGenerators,arg); end;;
gap> D219supplied := Mod2RingGensAndRels(219,3,D219relations.resolution,D219relations.ring.generators.classes,true,12);;
gap> Mod2RingGenerators := originalRingGenerators;;
gap> if ringGeneratorsWasReadOnly then MakeReadOnlyGlobal("Mod2RingGenerators"); fi;;
gap> [suppliedGen6FallbackCalls,D219supplied.generatorDegrees = D219relations.ring.generatorDegrees];
[ 0, true ]

# IT198 has no target monomials in degrees 3--5.  Empty reducer targets remain
# valid, and its only stored relation is the square of its degree-3 generator.
gap> relationReducerTrace := [];; relationReducerStats := [];;
gap> D198relations := SGC_CohomologyData(198);;
gap> [D198relations.ring.generatorDegrees,
>     Filtered(D198relations.ring.relations,x->x[1]=6),
>     relationReducerTrace];
[ [ 3 ], [ [ 6, [ [ [ 2 ] ] ] ] ], [ 3, 4, 5, 6 ] ]

# Preserve the raw degree-6 path used by a length-6 resolution (n=5).  Its
# accepted relation and its position in the degree sequence are pinned.
gap> relationReducerTrace := [];; relationReducerStats := [];;
gap> R221raw := SGC_ResolutionSpaceGroup(SpaceGroupIT(3,221),6);;
gap> G221raw := Mod2RingGenerators(R221raw,4,3);;
gap> D221raw := Mod2RingGensAndRels(221,3,R221raw,G221raw,true);;
gap> [D221raw.maxRelationDegree,D221raw.generatorDegrees,
>     Filtered(D221raw.relations,x->x[1]=6),relationReducerTrace];
[ 5, [ 1, 1, 1, 2, 2, 3, 3 ],
  [ [ 6, [ [ [ 0, 0, 0, 0, 0, 1, 1 ] ] ] ] ], [ 3, 4, 5, 6 ] ]

# IT136 is the canonical-row stress case from the performance audit.  Product
# counts exclude the legacy all-zero sentinel row; ranks pin the exact spans.
gap> relationReducerTrace := [];; relationReducerStats := [];;
gap> D136relations := SGC_CohomologyData(136);;
gap> relationReducerStats;
[ [ 3, 9, 6 ], [ 4, 45, 26 ], [ 5, 149, 76 ], [ 6, 399, 176 ],
  [ 7, 926, 352 ], [ 8, 1983, 656 ] ]

# An explicit cap may be any supported integer at least 2.  Supplied genuine
# generator slots remain visible through that cap, including empty high slots.
gap> relationReducerTrace := [];; relationReducerStats := [];;
gap> R1degree13 := SGC_ResolutionSpaceGroup(SpaceGroupIT(3,1),14);;
gap> G1degree13 := Mod2RingGenerators(R1degree13,4,3);;
gap> for r in [5..13] do G1degree13[r] := []; od;;
gap> D1degree13 := Mod2RingGensAndRels(fail,3,R1degree13,G1degree13,true,13,GENNAMES[1]);;
gap> [D1degree13.maxRelationDegree,
>     Length(D1degree13.generators.classes),
>     D1degree13.generatorDimensions,
>     relationReducerTrace=[3..13],
>     IsBound(D1degree13.relationsByDegree[13]),
>     D1degree13.relations,
>     D1degree13.bases];
[ 13, 13, [ 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], true, true,
  [ [ 2, [ [ [ 2, 0, 0 ] ], [ [ 0, 2, 0 ] ], [ [ 0, 0, 2 ] ] ] ] ],
  [ [ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ],
    [ [ 1, 1, 0 ], [ 1, 0, 1 ], [ 0, 1, 1 ] ],
    [ [ 1, 1, 1 ] ], [  ] ] ]

# The same untagged resolution may stop at any explicit cap supported by it.
# Degree 2 has no generic reducer, and only available basis states are returned.
gap> relationReducerTrace := [];; relationReducerStats := [];;
gap> D1degree2 := Mod2RingGensAndRels(fail,3,R1degree13,G1degree13,true,2,GENNAMES[1]);;
gap> [D1degree2.maxRelationDegree,
>     Length(D1degree2.generators.classes),
>     D1degree2.generatorDimensions,
>     relationReducerTrace=[],
>     IsBound(D1degree2.relationsByDegree[2]),
>     D1degree2.relations,
>     D1degree2.bases];
[ 2, 13, [ 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], true, true,
  [ [ 2, [ [ [ 2, 0, 0 ] ], [ [ 0, 2, 0 ] ], [ [ 0, 0, 2 ] ] ] ] ],
  [ [ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ],
    [ [ 1, 1, 0 ], [ 1, 0, 1 ], [ 0, 1, 1 ] ] ] ]

gap> SGC_SolverCacheList := savedSolverCacheList;;
gap> if relationObserverWasBound then SGC_RELATION_REDUCER_OBSERVER := relationSavedObserver; if relationObserverWasReadOnly then MakeReadOnlyGlobal("SGC_RELATION_REDUCER_OBSERVER"); fi; else UnbindGlobal("SGC_RELATION_REDUCER_OBSERVER"); fi;;
gap> STOP_TEST("relations.tst", 0);
