# smoke.tst — representative public API coverage.

gap> START_TEST("SpaceGroupCohomology: smoke test");
gap> SetRecursionTrapInterval(100000000);;

# Ordinary cohomology-ring presentation.
gap> GroupCohomologyMod2(16);;
Z2[Ac,Acp,Ax,Ay,Az]/<R2>
R2:  Ax^2+Ac.Ax+Acp.Ax  Ay^2+Ac.Ay  Az^2+Acp.Az

# Wyckoff table for a group with explicit degree-3 generators.
gap> WPCohomologyTable(22);;
4a Ac^2.Acp+Ac.Acp^2+Ac.Acp.Axy+Ac.Acp.Axz+Cxyz
4b Cxyz
4c Cc+Cxyz
4d Ac.Acp.Axy+Ac.Acp.Axz+Cc+Cxyz

# A nontrivial class in a group with a degree-4 generator.
gap> WPCohomologyClass(108,["4a"]);;
Amp.Ba+Am.Ba+Axyz.Ba+Amp.Bzxy

# End-to-end ring and LSM output for a group with degree-6 generators.
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

gap> STOP_TEST("smoke.tst", 0);
