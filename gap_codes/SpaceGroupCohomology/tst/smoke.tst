# smoke.tst — sanity check for the SpaceGroupCohomology package.

gap> START_TEST("SpaceGroupCohomology: smoke test");

# One assertion: all 230 data entries are bound after package load.
gap> Length(PGGens230) = 230 and Length(IWP) = 230 and Length(GENNAMES) = 230 and Length(funcs230) = 230;
true

# Test 1: Simple triclinic group
gap> SpaceGroupCohomologyRingGapInterface(1);
===========================================
Mod-2 Cohomology Ring of Group No. 1:
Z2[Ax,Ay,Az]/<R2>
R2:  Ax^2  Ay^2  Az^2  
===========================================
LSM:
1a Ax.Ay.Az
true

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

# Test 3: Tetragonal group
gap> SpaceGroupCohomologyRingGapInterface(89);
===========================================
Mod-2 Cohomology Ring of Group No. 89:
Z2[Acp,Acpp,Axy,Az,Ba,Bxy]/<R2,R3,R4>
R2:  Acp.Acpp  Acpp.Axy  Az^2+Acp.Az+Acpp.Az  
R3:  Acp.Bxy+Acp.Axy^2+Axy^3+Axy.Ba  Axy.Bxy+Acp.Axy^2+Axy^3+Axy.Ba  
R4:  Bxy^2+Ba.Bxy  
===========================================
LSM:
1a Acp^2.Axy+Acp.Axy.Az+Axy^2.Az+Axy^3+Acp.Ba+Acpp.Ba+Axy.Ba+Az.Ba+Acpp.Bxy+Az\
.Bxy
1b Acp.Axy.Az+Axy^2.Az+Az.Ba+Az.Bxy
1c Acp.Axy^2+Axy^3+Axy.Ba+Acpp.Bxy+Az.Bxy
1d Az.Bxy
2e Acp^2.Axy+Acp.Axy^2+Acp.Axy.Az+Axy^2.Az
2f Acp.Axy.Az+Axy^2.Az
true

gap> STOP_TEST("smoke.tst", 0);
