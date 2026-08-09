# smoke.tst — sanity check for the SpaceGroupCohomology package.

gap> START_TEST("SpaceGroupCohomology: smoke test");

# One assertion: all 230 data entries are bound after package load.
gap> Length(PGGens230) = 230 and Length(IWP) = 230 and Length(GENNAMES) = 230 and Length(funcs230) = 230;
true

# Test 1: Simple triclinic group.  Instrument the existing pinned call to
# guarantee that the compatibility wrapper constructs exactly one resolution,
# and retain that resolution for the focused public-presentation checks below.
gap> originalResolutionBuilder := SGC_ResolutionSpaceGroup;;
gap> resolutionCallCount := 0;;
gap> Rapi := fail;;
gap> resolutionBuilderWasReadOnly := IsReadOnlyGlobal("SGC_ResolutionSpaceGroup");;
gap> if resolutionBuilderWasReadOnly then MakeReadWriteGlobal("SGC_ResolutionSpaceGroup"); fi;;
gap> SGC_ResolutionSpaceGroup := function(G,n) resolutionCallCount := resolutionCallCount+1; Rapi := originalResolutionBuilder(G,n); return Rapi; end;;
gap> SpaceGroupCohomologyRingGapInterface(1);
===========================================
Mod-2 Cohomology Ring of Group No. 1:
Z2[Ax,Ay,Az]/<R2>
R2:  Ax^2  Ay^2  Az^2
===========================================
LSM:
1a Ax.Ay.Az
true
gap> SGC_ResolutionSpaceGroup := originalResolutionBuilder;;
gap> if resolutionBuilderWasReadOnly then MakeReadOnlyGlobal("SGC_ResolutionSpaceGroup"); fi;;
gap> resolutionCallCount;
1

# The batched cup helper reuses one left chain map and preserves right-factor
# order.  This is the degree-one IT1 fixture used throughout the ring pipeline.
gap> C1 := CR_Mod2CocyclesAndCoboundaries(Rapi,1,true);; C2 := CR_Mod2CocyclesAndCoboundaries(Rapi,2,true);;
gap> BatchedCupFixture := function() if not IsBoundGlobal("Mod2CupProducts") then return fail; fi; return ValueGlobal("Mod2CupProducts")(Rapi,[1,0,0],CohomologyBasis([1,1,1]),1,1,C1,C1,C2); end;;
gap> BatchedCupFixture();
[ [ 0, 0, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ]
gap> originalChainMapBuilder := CR_ChainMapFromCocycle;; chainMapBuildCount := 0;; chainMapEvaluationCount := 0;;
gap> chainMapBuilderWasReadOnly := IsReadOnlyGlobal("CR_ChainMapFromCocycle");;
gap> if chainMapBuilderWasReadOnly then MakeReadWriteGlobal("CR_ChainMapFromCocycle"); fi;;
gap> CR_ChainMapFromCocycle := function(arg) local f; chainMapBuildCount := chainMapBuildCount+1; f := CallFuncList(originalChainMapBuilder,arg); return function(x) chainMapEvaluationCount := chainMapEvaluationCount+1; return f(x); end; end;;
gap> batchedCounted := Mod2CupProducts(Rapi,[1,0,0],CohomologyBasis([1,1,1]),1,1,C1,C1,C2);;
gap> CR_ChainMapFromCocycle := originalChainMapBuilder;;
gap> if chainMapBuilderWasReadOnly then MakeReadOnlyGlobal("CR_ChainMapFromCocycle"); fi;;
gap> [batchedCounted,chainMapBuildCount,chainMapEvaluationCount];
[ [ [ 0, 0, 0 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ], 1, 3 ]

# Public 2.1 presentation APIs.  The IT-only forms reuse Rapi through a narrow
# substitution of the resolution factory; all cohomology and formatting work
# remains real.  The supplied-resolution forms use the same Rapi directly.
gap> originalResolutionForIT := SGC_ResolutionForIT;;
gap> resolutionForITWasReadOnly := IsReadOnlyGlobal("SGC_ResolutionForIT");;
gap> if resolutionForITWasReadOnly then MakeReadWriteGlobal("SGC_ResolutionForIT"); fi;;
gap> SGC_ResolutionForIT := function(IT) return Rapi; end;;
gap> GroupCohomologyMod2(1);;
Z2[Ax,Ay,Az]/<R2>
R2:  Ax^2  Ay^2  Az^2
gap> GroupCohomologyMod2(1,Rapi);;
Z2[Ax,Ay,Az]/<R2>
R2:  Ax^2  Ay^2  Az^2
gap> WPCohomologyTable(1);;
1a Ax.Ay.Az
gap> WPCohomologyTable(1,Rapi);;
1a Ax.Ay.Az
gap> WPCohomologyClass(1,["1a"]);;
Ax.Ay.Az
gap> WPCohomologyClass(1,Rapi,["1a"]);;
Ax.Ay.Az
gap> WPCohomologyClass(1,Rapi,["1a","1a"]);;
0

# Prepared Context data must avoid the unnecessary generic generator pass for
# IT1, and the public APIs must accept one caller-owned Context/WPData pair.
gap> originalMod2RingGenerators := Mod2RingGenerators;;
gap> mod2RingGeneratorDegrees := [];;
gap> mod2RingGeneratorsWasReadOnly := IsReadOnlyGlobal("Mod2RingGenerators");;
gap> if mod2RingGeneratorsWasReadOnly then MakeReadWriteGlobal("Mod2RingGenerators"); fi;;
gap> Mod2RingGenerators := function(arg) Add(mod2RingGeneratorDegrees,arg[2]); return CallFuncList(originalMod2RingGenerators,arg); end;;
gap> ContextApi := SGC_CohomologyData(1,Rapi);;
gap> Mod2RingGenerators := originalMod2RingGenerators;;
gap> if mod2RingGeneratorsWasReadOnly then MakeReadOnlyGlobal("Mod2RingGenerators"); fi;;
gap> mod2RingGeneratorDegrees;
[  ]

# Once the source IWP rows already span H3, WP construction must not evaluate
# any additional topological-invariant candidates.  The synthetic point-group
# matrix makes continued scanning observable without changing the IT1 IWP row.
gap> ContextEarlyExit := ShallowCopy(ContextApi);; originalTopoInvdeg3 := ContextApi.topoInvdeg3;; topoInvdeg3Calls := 0;; nonemptyIWPCount := Number(IWP[1],x->x[2]<>[]);;
gap> ContextEarlyExit.pointGroupMatrices := [[[-1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,1]]];;
gap> ContextEarlyExit.topoInvdeg3 := function(arg) topoInvdeg3Calls := topoInvdeg3Calls+1; if topoInvdeg3Calls <= nonemptyIWPCount then return CallFuncList(originalTopoInvdeg3,arg); fi; return List([1..Length(arg[2])],x->0); end;;
gap> WPDataEarlyExit := SGC_WPCohomologyData(ContextEarlyExit);;
gap> [topoInvdeg3Calls,Number(IWP[1],x->x[2]<>[])];
[ 1, 1 ]

# The two-argument topological-invariant form is exactly the identity transform;
# an explicitly supplied non-identity transform retains its old semantics.
gap> topoFixture := IWP[1][1][2];; topoLetters := ContextApi.base3Letters;;
gap> [ContextApi.topoInvdeg3(topoFixture,topoLetters)=ContextApi.topoInvdeg3(topoFixture,topoLetters,[[1]]),ContextApi.topoInvdeg3(topoFixture,topoLetters,[[0]])];
[ true, [ 0 ] ]
gap> WPDataApi := SGC_WPCohomologyData(ContextApi);;
gap> [ContextApi.isSGCCohomologyContext,WPDataApi.isSGCWPCohomologyData,IsIdenticalObj(WPDataApi.context,ContextApi)];
[ true, true, true ]
gap> originalCohomologyData := SGC_CohomologyData;; originalWPCohomologyData := SGC_WPCohomologyData;;
gap> cohomologyDataCalls := 0;; wpCohomologyDataCalls := 0;;
gap> cohomologyDataWasReadOnly := IsReadOnlyGlobal("SGC_CohomologyData");; wpCohomologyDataWasReadOnly := IsReadOnlyGlobal("SGC_WPCohomologyData");;
gap> if cohomologyDataWasReadOnly then MakeReadWriteGlobal("SGC_CohomologyData"); fi;; if wpCohomologyDataWasReadOnly then MakeReadWriteGlobal("SGC_WPCohomologyData"); fi;;
gap> SGC_CohomologyData := function(arg) cohomologyDataCalls := cohomologyDataCalls+1; return CallFuncList(originalCohomologyData,arg); end;;
gap> SGC_WPCohomologyData := function(Context) wpCohomologyDataCalls := wpCohomologyDataCalls+1; return originalWPCohomologyData(Context); end;;
gap> GroupCohomologyMod2(ContextApi);;
Z2[Ax,Ay,Az]/<R2>
R2:  Ax^2  Ay^2  Az^2
gap> WPCohomologyTable(ContextApi);;
1a Ax.Ay.Az
gap> WPCohomologyTable(ContextApi,WPDataApi);;
1a Ax.Ay.Az
gap> WPCohomologyClass(ContextApi,["1a"]);;
Ax.Ay.Az
gap> WPCohomologyClass(ContextApi,WPDataApi,["1a","1a"]);;
0
gap> SGC_CohomologyData := originalCohomologyData;; SGC_WPCohomologyData := originalWPCohomologyData;;
gap> if cohomologyDataWasReadOnly then MakeReadOnlyGlobal("SGC_CohomologyData"); fi;; if wpCohomologyDataWasReadOnly then MakeReadOnlyGlobal("SGC_WPCohomologyData"); fi;;
gap> [cohomologyDataCalls,wpCohomologyDataCalls];
[ 0, 2 ]
gap> SGC_ResolutionForIT := originalResolutionForIT;;
gap> if resolutionForITWasReadOnly then MakeReadOnlyGlobal("SGC_ResolutionForIT"); fi;;

# Multi-digit labels stay in source IWP order, aliases print their canonical
# labels and classes, and source-declared empty positions print zero.
gap> WPCohomologyTable(70);;
8a Acp^2.Ac+Acp.Ac^2+Cc
8b Cc
16c Ai^3+Ai.Bxyxzyz
16d Ai.Bxyxzyz
gap> WPCohomologyTable(143);;
1a Az.Bxy
1b Az.Bxy
1c Az.Bxy
gap> WPCohomologyTable(147);;
1a Ai^2.Az+Ai^3+Ai.Bxy+Az.Bxy
1b Ai^2.Az+Az.Bxy
2d 0
3e Ai.Bxy+Az.Bxy
3f Az.Bxy

# Test 2: Orthorhombic group
gap> SpaceGroupCohomologyRingGapInterface(16);
===========================================
Mod-2 Cohomology Ring of Group No. 16:
Z2[Ac,Acp,Ax,Ay,Az]/<R2>
R2:  Ax^2+Ac.Ax+Acp.Ax  Ay^2+Ac.Ay  Az^2+Acp.Az
===========================================
LSM:
1a Ac^2.Acp+Ac.Acp^2+Ac.Acp.Ax+Ac.Acp.Ay+Ac.Acp.Az+Ac.Ax.Az+Ac.Ay.Az+Ac^2.Az+A\
cp.Ax.Ay+Acp^2.Ay+Acp.Ay.Az+Ax.Ay.Az
1b Ac.Acp.Ax+Ac.Ax.Az+Acp.Ax.Ay+Ax.Ay.Az
1c Ac.Acp.Ay+Ac.Ay.Az+Acp.Ax.Ay+Acp^2.Ay+Acp.Ay.Az+Ax.Ay.Az
1d Ac.Acp.Az+Ac.Ax.Az+Ac.Ay.Az+Ac^2.Az+Acp.Ay.Az+Ax.Ay.Az
1e Acp.Ax.Ay+Ax.Ay.Az
1f Ac.Ax.Az+Ax.Ay.Az
1g Ac.Ay.Az+Acp.Ay.Az+Ax.Ay.Az
1h Ax.Ay.Az
true

# Test 3: Orthorhombic group with degree-3 generators
gap> SpaceGroupCohomologyRingGapInterface(22);
===========================================
Mod-2 Cohomology Ring of Group No. 22:
Z2[Ac,Acp,Axy,Axz,Cc,Cxyz]/<R2,R4,R6>
R2:  Acp.Axz+Ac.Axy  Axy^2+Acp.Axy  Axz^2+Ac.Axz
R4:  Ac.Cxyz+Ac.Cc+Axz.Cc  Acp.Cxyz+Acp.Cc+Axy.Cc  Axy.Cxyz  Axz.Cxyz
R6:  Cc^2+Ac^2.Acp.Cc+Ac.Acp^2.Cc  Cc.Cxyz+Ac^2.Acp.Cc+Ac^2.Axy.Cc+Ac.Acp^2.Cc\
+Ac.Acp.Axy.Cc  Cxyz^2+Ac^2.Acp.Cc+Ac^2.Axy.Cc+Ac.Acp^2.Cc+Ac.Acp.Axy.Cc
===========================================
LSM:
4a Ac^2.Acp+Ac.Acp^2+Ac.Acp.Axy+Ac.Acp.Axz+Cxyz
4b Cxyz
4c Cc+Cxyz
4d Ac.Acp.Axy+Ac.Acp.Axz+Cc+Cxyz
true

# Test 4: Tetragonal group with degree-4 generator
gap> mod2RingGeneratorDegrees := [];;
gap> if mod2RingGeneratorsWasReadOnly then MakeReadWriteGlobal("Mod2RingGenerators"); fi;;
gap> Mod2RingGenerators := function(arg) Add(mod2RingGeneratorDegrees,arg[2]); return CallFuncList(originalMod2RingGenerators,arg); end;;
gap> SpaceGroupCohomologyRingGapInterface(108);
===========================================
Mod-2 Cohomology Ring of Group No. 108:
Z2[Amp,Am,Axyz,Ba,Bzxy,Dd]/<R2,R3,R4,R5,R6,R8>
R2:  Amp.Am  Amp.Axyz  Axyz^2+Amp^2+Am^2
R3:  Am.Bzxy+Am^2.Axyz+Am.Axyz^2+Axyz.Ba  Axyz.Bzxy+Am^2.Axyz+Am.Axyz^2+Axyz.B\
a
R4:  Bzxy^2+Ba.Bzxy
R5:  Amp.Dd  Axyz.Dd+Am.Dd
R6:  Bzxy.Dd+Ba.Dd
R8:  Dd^2
===========================================
LSM:
4a Amp.Ba+Am.Ba+Axyz.Ba+Amp.Bzxy
4b Amp.Bzxy
true
gap> Mod2RingGenerators := originalMod2RingGenerators;;
gap> if mod2RingGeneratorsWasReadOnly then MakeReadOnlyGlobal("Mod2RingGenerators"); fi;;
gap> mod2RingGeneratorDegrees;
[ 4 ]

# Test 5: Cubic group with degree-6 generators
gap> mod2RingGeneratorDegrees := [];;
gap> if mod2RingGeneratorsWasReadOnly then MakeReadWriteGlobal("Mod2RingGenerators"); fi;;
gap> Mod2RingGenerators := function(arg) Add(mod2RingGeneratorDegrees,arg[2]); return CallFuncList(originalMod2RingGenerators,arg); end;;
gap> SpaceGroupCohomologyRingGapInterface(219);
===========================================
Mod-2 Cohomology Ring of Group No. 219:
Z2[Am,Ba,Bxyxzyz,Ca,Cb,Fd1,Fd2]/<R3,R4,R5,R6,R7,R8,R9,R12>
R3:  Am^3
R4:  Am^2.Ba  Am^2.Bxyxzyz  Am.Ca  Am.Cb  Bxyxzyz^2+Ba.Bxyxzyz
R5:  Ba.Cb+Bxyxzyz.Ca  Bxyxzyz.Cb+Bxyxzyz.Ca
R6:  Cb^2+Ca.Cb
R7:  Am.Fd1  Am.Fd2
R8:  Bxyxzyz.Fd1  Bxyxzyz.Fd2+Ba.Fd2
R9:  Cb.Fd1  Cb.Fd2+Ca.Fd2
R12:  Fd1^2+Ba^3.Ca^2+Ba^2.Bxyxzyz.Ca^2+Ca^4+Ca^3.Cb+Ca^2.Fd1  Fd1.Fd2  Fd2^2+\
Ba^2.Bxyxzyz.Ca^2+Ca^3.Cb+Ca^2.Fd2
===========================================
LSM:
8a Am.Ba+Am.Bxyxzyz+Ca+Cb
8b Am.Bxyxzyz+Cb
24c Am.Bxyxzyz
24d Am.Ba+Am.Bxyxzyz
true
gap> Mod2RingGenerators := originalMod2RingGenerators;;
gap> if mod2RingGeneratorsWasReadOnly then MakeReadOnlyGlobal("Mod2RingGenerators"); fi;;
gap> mod2RingGeneratorDegrees;
[ 6 ]

gap> STOP_TEST("smoke.tst", 0);
