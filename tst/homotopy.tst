gap> START_TEST("SpaceGroupCohomology: integral contracting homotopy");
gap> LoadPackage("SpaceGroupCohomology");;
gap> G := SpaceGroupBBNWZ(3,16);;
gap> P := SGC_CrystallographicComplexWithGeometry(G);;
gap> SGC_AttachCellularContraction(P);;
gap> ok := true;;
gap> for q in [0..3] do
>     for i in [1..P!.dimension(q)] do
>         for g in [1..Minimum(Length(P!.elts),4)] do
>             if not SGC_NonFreeHomotopyCheck(P, q, [i,g]) then
>                 ok := false; Print("P-level failure: ", [q,i,g], "\n");
>             fi;
>         od;
>     od;
> od;
gap> ok;    # d h + h d = id on the non-free complex, over Z
true
gap> R0 := FreeGResolution(P, 4);;
gap> R := SGC_FreeGResolutionWithHomotopy(P, 4);;
gap> EvaluateProperty(R, "characteristic");
0
gap> List([0..4], i -> R0!.dimension(i)) = List([0..4], i -> R!.dimension(i));
true
gap> ok := true;;
gap> for n in [1..4] do
>     for i in [1..R0!.dimension(n)] do
>         if R0!.boundary(n,i) <> R!.boundary(n,i) then
>             ok := false; Print("boundary mismatch: ", [n,i], "\n");
>         fi;
>     od;
> od;
gap> ok;    # fork reproduces HAP's FreeGResolution boundaries exactly
true
gap> ok := true;;
gap> for k in [0..2] do
>     for i in [1..R!.dimension(k)] do
>         for g in [1..Minimum(Length(R!.elts),4)] do
>             if not SGC_TotalHomotopyCheck(R, k, [i,g]) then
>                 ok := false; Print("R-level failure: ", [k,i,g], "\n");
>             fi;
>         od;
>     od;
> od;
gap> ok;    # d h + h d = id (resp. id - [1,1_G] at k=0) on the free resolution
true
gap> STOP_TEST("homotopy.tst", 0);
