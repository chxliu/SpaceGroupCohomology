# Code Review — `gap_codes/SpaceGroupCohomologyFunctions.gi`

**Reviewed file:** `gap_codes/SpaceGroupCohomologyFunctions.gi` (4,519 lines, GAP/HAP)
**Date:** 2026-05-21
**Method:** Static review of the full file. GAP language semantics that one finding depends on were verified empirically with an isolated snippet (`gap -q`). The code was **not** run end-to-end (that needs HAP plus the 14 MB `SpaceGroupCohomologyData.gi` / `Space_Group_Cocycles.gi`). No source files were modified.

---

## Conclusion (read this first)

The mathematical core (`CR_Mod2CocyclesAndCoboundaries`, `Mod2CupProduct`, the cocycle/coboundary and cup-product machinery) looks correct, and the **primary entry path** — calling `SpaceGroupCohomologyRingGapInterface(IT)`, which in turn calls `Mod2RingGensAndRels(IT,3,R,Gens)` with **four** arguments — avoids every correctness bug below.

However, there are:

- **2 real (latent) crash bugs** on documented/alternate call paths,
- several **fragility / coupling** problems that make the code correct only for one exact environment,
- a number of **performance hotspots** (recomputation inside hot loops, combinatorial loops), and
- large-scale **maintainability** debt (~2,900 lines of copy-pasted degree blocks, many dead/unused locals).

Severity legend: 🔴 blocking/correctness · 🟡 important · 🟢 minor/nit · ⚡ performance · 📚 note.

| # | Severity | Location | Issue |
|---|----------|----------|-------|
| B1 | 🔴 | `Mod2RingGensAndRels` L591 | Uses unassigned local `n` on the `Length(arg)<=2` path → crash |
| B2 | 🔴 | `TopoInvdeg3` L4063–4080 | `Length(gs)=4` with no matching group → unassigned `vallist` → crash (no `else`) |
| B3 | 🟡 | `CR_Mod2CocyclesAndCoboundaries` L162,267 | 3rd argument silently ignored; `toggle=false` branch is dead code |
| B4 | 🟡 | `ClassToCycle` L301, `TopoInvdeg3` L4017 | Loop variables `j` / `i` not declared local (closure-captured from parent) |
| D1 | 🟡 | header L558–562 | Doc says `arg[2] = n`; code uses `arg[2]` as `spacedim` — and the examples would crash (see B1) |
| D2 | 🟢 | two main functions | `arg[2]` means `n` in `Mod2RingGenerators` but `spacedim` in `Mod2RingGensAndRels` |
| F1 | 🟡 | L4225–4276, L2526 | Hardcoded cocycle vectors / H⁶ dimension table tied to one exact resolution |
| F2 | 🟡 | whole file | Depends on external globals `GENNAMES`, `PGGens230`, `funcs230`, `IWP` never loaded here |
| F3 | 🟡 | L4121–4167 | Point-group enumeration assumes the 2nd generator has order 2 (`o2 in [0..1]`) |
| F4 | 🟢 | L3916,3934,3945 | `Position(...) = fail` not guarded in `MatToPow`/`Invofg`/`Prodg1g2Pow` |
| P1 | ⚡ | r=4..8 prep loops | `Concatenation([Gen0],GensLett)` rebuilt on every nested iteration |
| P2 | ⚡ | every coboundary test | `Cohomology(TR,k)` + zero-vector rebuilt inside the inner loop |
| P3 | ⚡ | relation reduction | `Position(RelReduceLett,…)` linear scan in tight loops |
| P4 | ⚡ | r=2, r=3 | `Mod2CupProduct` recomputes the chain map per `v` (not hoisted as at r≥4) |
| P5 | ⚡ | r≥5 | `RankMatrix` recomputed over a growing matrix for every candidate |
| P6 | ⚡ | LSM loop L4425–4477 | `PGind × [-2..2]³ × PGind × [-2..2]³` combinatorial blow-up |
| P7 | ⚡ | L172–214, L1703–1715 | Boundary-matrix build rescans sparse `Boundary(...)` once per column — scatter directly instead |
| M1 | 🟢 | L716–3551 | ~2,900 lines of near-duplicated degree-2…8 blocks |
| M2 | 🟢 | many | Unused locals + large dead/commented blocks (e.g. `if false then …` L4297–4361) |
| N1 | 🟢 | `GF2ToZ` L51–58 | Verbose 8-line loop replaceable by idiomatic `List(v, IntFFE)` |

---

## 🔴 Correctness bugs

### B1 — `Mod2RingGensAndRels` references `n` before it is assigned (L591)

In the `Length(arg) <= 2` branch:

```gap
if Length(arg)<=2 then
    ...
    R := ResolutionAlmostCrystalGroup(GG, n+1);   # L591  — n is unbound here
    Gens := Mod2RingGenerators(R,4,spacedim);
else
    R := arg[3];
    Gens := arg[4];
fi;
...
n := Length(Size(R))-1;                            # L688  — n first assigned here
```

`n` is a local (declared L544) and is only assigned at **L688**, which is *after* L591. On the `<=2`-argument path GAP aborts with
`Error, Variable: 'n' must have an assigned value`.

I verified this is exactly how GAP treats an unassigned local with a minimal snippet:

```
calling f(89) ...
Error, Variable: 'n' must have an assigned value in
  R := n + 1; at t.g:4
```

**Impact:** the documented standalone calls `Mod2RingGensAndRels(89)` and `Mod2RingGensAndRels(89,3)` (examples at L559–560) crash immediately. The code only works through the 4-argument call made by `SpaceGroupCohomologyRingGapInterface` (L4373), where the `else` branch is taken and L591 is never reached.

**Fix direction:** decide the resolution length explicitly on this path (e.g. set `n` from the desired highest relation degree, or build `R` with a fixed length such as `n+1` after choosing `n`), or simply delete the dead `<=2` branch if the standalone entry point is not meant to be supported.

### B2 — `TopoInvdeg3` can return with `vallist` unassigned (L4063–4080)

The `Length(gs) = 4` branch only assigns `vallist` for space groups 19/198, 29, and 33:

```gap
elif Length(gs) = 4 then
    if  ... then       vallist := ...;     # No. 19 & 198
    elif ... then      vallist := ...;     # No. 29
    elif ... then      vallist := ...;     # No. 33
    fi;                                     # <-- no `else`
else
    Print("Wrong in checking ...");
fi;
return GF2ToZ((solrels*vallist)*Z(2));     # L4080 — vallist may be unbound
```

Unlike the `Length(gs)=3` branch (which has an `else` printing an error), the 4-element branch has **no fallback**. A 4-element entry in `IWP[IT]` for any other group reaches L4080 with `vallist` unbound → crash. The other branches at least print a diagnostic, but they *still* fall through to L4080 with `vallist` unset after printing, so a non-invariant input is also a latent crash rather than a clean error.

**Fix direction:** add an `else` that sets a defined fallback (or `return fail`/raises a clear `Error`) in both the `=3` and `=4` branches, and guard L4080 against an unassigned `vallist`.

---

## 🟡 Important issues

### B3 — Ignored 3rd argument and dead branch in `CR_Mod2CocyclesAndCoboundaries`

```gap
R := arg[1];
n := arg[2];
...
toggle := true;          # L162 — hardcoded; arg[3] is never read
...
if toggle=false then     # L267 — can never be true
   return rec( ... cocycleToClass:=fail, classToCocycle:=fail );
fi;
```

The function is called everywhere as `CR_Mod2CocyclesAndCoboundaries(R, p, true)` (e.g. L351, L448, L694, L4202), implying the 3rd argument is meant to select behavior. It is never read; `toggle` is hardwired to `true`, so the `toggle=false` block (L267–276) is unreachable. Behavior happens to be the intended one, so this is not a wrong result — but it is misleading dead code and a trap for the next maintainer.

**Fix direction:** either `if Length(arg) > 2 then toggle := arg[3]; fi;`, or drop the parameter and the dead branch and the `true` at every call site.

### B4 — Inner-function loop variables not declared `local` (closure capture)

- `ClassToCycle` (L293) declares `local v,w,i,temp;` but its inner loop uses `j` (L301): `for j in [1..cohdim] do`. `j` resolves to the **enclosing** `CR_Mod2CocyclesAndCoboundaries` local `j` (L155).
- `TopoInvdeg3` (L4009) declares `local gs,letters,solrels,vallist;` but uses `i` at L4017: `for i in [1..Length(letters)] do`. `i` resolves to the enclosing `SpaceGroupCohomologyRingGapInterface` local `i` (L3891).

These work today only because no *active* outer loop over the same variable wraps the call sites. They are latent bugs: any future caller that invokes these closures from inside an `i`/`j` loop will have its loop counter silently corrupted. GAP would normally warn about this; declaring the variables local removes the hazard.

Related (correctness-OK but wrong-looking): in `ClassToCycle` the assignment `w[i] := temp mod 2;` (L303) sits *inside* the `j` loop, so it is rewritten `cohdim` times per `i`. The final value is correct, but it should be one assignment after the loop.

### D1 — Header documentation contradicts the code (L558–562)

```gap
#Standard input: arg[1] = IT (# of space group), arg[2] = n (relations up to deg(n) is calculated)
#e.g.: Mod2RingGensAndRels(89,3);
```

In the code `arg[2]` is assigned to **`spacedim`** (L572), not `n`; the relation degree `n` is derived from the resolution length (`n := Length(Size(R))-1`, L688). So the comment is wrong, and (per B1) the example `Mod2RingGensAndRels(89,3)` would crash anyway. Update the doc to reflect `arg = [IT, spacedim, R, Gens]` and remove/repair the standalone examples.

### F1 — Hardcoded values pinned to one exact resolution

- **Generator vectors** for specific groups are hardcoded with the warning *"Must use standard resolution!!"*:
  - `Gen3` for IT 225/227/229 (L4228–4233)
  - `Gen4` for IT 108/120/140/142/230 (L4266–4274)
- **H⁶ dimension table** in the raw method: `[62,11,31,26,45,20,19,6,40,7][IT-220]` (L2526, L2529).

These literal cocycle coordinates and dimensions are only valid for the precise basis that this HAP version's resolution produces, and only for IT in 221–230 (for the table). If HAP changes its resolution, or the function is fed any other resolution, these silently become wrong. `IT-220` with `IT <= 220` would also index out of range. This is acceptable for a frozen reproducibility artifact but should be loudly documented as such; ideally the dimensions would be computed rather than tabulated.

### F2 — Undeclared external global dependencies

The file does `LoadPackage("HAP");` (L40) but the data `Read(...)` is commented out (L4099, L42). Yet it references the globals `GENNAMES`, `PGGens230`, `funcs230`, and `IWP` **56 times** (e.g. L3561, L3617, L3718, L4103–4104, L4409). These must come from `SpaceGroupCohomologyData.gi` / `Space_Group_Cocycles.gi`. If a user loads only this file, every such reference is an undefined-global error. There is no guard such as `if not IsBound(GENNAMES) then ... fi;`. Recommend documenting the required load order at the top, or adding an `IsBound` precondition check.

### F3 — Point-group enumeration assumes generator #2 has order 2

In `SpaceGroupCohomologyRingGapInterface`, the loops that enumerate point-group elements hardcode the second generator's exponent range as `for o2 in [0..1]` (L4123, L4131, L4141, L4153) with the comment *"the ordering of the group generators must strictly follow this condition"*. All other generators use `[0..Order(PGGen33[k])-1]`. If any group's 2nd point-group generator has order > 2 in the `PGGens230` data, `PGind`/`PGMat33` are built incompletely and downstream `Position(PGind, …)` lookups return `fail` (then crash, see F4) or produce wrong invariants — with no runtime check. Worth at least an assertion `Assert(0, Order(PGGen33[2]) = 2)`.

### F4 — Unguarded `Position(...)=fail` in matrix→power helpers

`MatToPow` (L3916 `i:=Position(PGMat33,mat33)`), `Invofg` (L3934), and `Prodg1g2Pow` (L3945) index `PGMatinv[Position(PGind, …)]` / `PGMatinv[i]` without checking for `fail`. For any element outside the enumerated set (e.g. if F3 is violated, or a translation outside the assumed range), `Position` returns `fail` and the array index crashes with an opaque message. A guarded error message would aid debugging.

---

## ⚡ Performance

### P1 — `Concatenation([Gen0],GensLett)` rebuilt inside nested loops
The relation-reduction preparation for r = 4,5,6,7,8 iterates with the list expression *inline* in each `for`:

```gap
for p in Concatenation([Gen0],GensLett) do
    for q in Concatenation([Gen0],GensLett) do
        for r in Concatenation([Gen0],GensLett) do
            for s in Concatenation([Gen0],GensLett) do ...
```
(e.g. L1055–1058, L1778–1822, L2613–2676, L3100–3186). The same list is reconstructed on every iteration of every level. For r=8 (up to 8 nested levels) this is a large, avoidable cost. Hoist once: `allGens := Concatenation([Gen0], GensLett);` and reuse. Likewise hoist `Sum(GenDim1to4)+1` out of loop bounds.

### P2 — `Cohomology(TR,k)` and the zero vector rebuilt per inner-loop iteration
The "is it a coboundary" test is written as, e.g.:
```gap
if cupped = List([1..Cohomology(TR,3)],x->0) then     # L909, and similarly L734,767,1130,1228,1445,1554,1882,1992,2130,2730,2842,2982,3240,3352,3492 ...
```
This recomputes `Cohomology(TR,k)` and allocates a fresh zero list **every time** the test runs (i.e. for every cup product in the double loops). Precompute once per degree: `zero_k := List([1..Cohomology(TR,k)], x->0);` and compare against `zero_k`. Depending on whether HAP caches `Cohomology`, this can be a major saving.

### P3 — Linear `Position(RelReduceLett, …)` lookups
Reducibility construction does `RelReduceVec[Position(RelReduceLett, u+w)] := 1;` repeatedly (L870, L1076, L1344, …). `Position` is O(length) and `RelReduceLett` grows large at high degree. A hashed lookup (e.g. build a record/`Dictionary` keyed by the monomial vector) would cut this from O(n) to O(1) per access.

### P4 — Low-degree cup products do not reuse the chain map
At r=2 (L732, L765) and r=3 (L887) the code calls `Mod2CupProduct(R,u,v,…)` inside the double loop, which recomputes `CR_ChainMapFromCocycle` for `u` on every `v`. The r≥4 blocks already optimize this with the explicit "part 1 / part 2" split (compute `uChainMap` once per `u`, then only the cheap inner sum per `v`, e.g. L1391–1417). Applying the same hoist at r=2/3 is a free win and would also remove an inconsistency.

### P5 — `RankMatrix` recomputed on a growing matrix
For r≥5 the reducibility test is `rk1 := RankMatrix(Concatenation(RelReduceMat*Z(2), [RelReduceVec]*Z(2)));` per candidate (L1482, L1591, L1919, L2284, …). Each call re-ranks the entire accumulated matrix. An incremental rank update (or a maintained row-echelon form / `SemiEchelonMat` you extend) would avoid the repeated full factorization.

### P6 — LSM search is combinatorial
The non-LSM topological-invariant search (L4425–4477) nests `for v2 in PGind` × `x2,y2,z2 in [-2..2]` and, in the mirror case, an inner `for v1 in PGind` × `x1,y1,z1 in [-2..2]`. Worst case ≈ `|PGind|·5³·|PGind|·5³`; for cubic groups (|PGind|=48) that is tens of millions of iterations, each running `TopoInvdeg3` (many `FuncVal` sums) and `SolutionMat`. This dominates runtime for high-symmetry groups. Consider early termination once `RankMatrix(Mat) = Length(Base3Lett)` (full rank reached → no further invariants possible), and narrowing the `[-2..2]` translation window to what is mathematically required.

Smaller items: `RankMatrix(Mat*Z(2))` is computed twice in the L4485 condition (store it); the identity `solrels` is rebuilt on every `TopoInvdeg3` call when the 3rd arg is omitted (L4016–4019); `Mod2CupProduct` recomputes `P/Q/N` cocycle bases when they are not passed (fine internally, costly if used standalone in a loop).

### P7 — Boundary matrices rescan the sparse boundary once per column
In `CR_Mod2CocyclesAndCoboundaries` the M1/M2 build (L172–214) is structured as a triple nest:

```gap
for i in [1..Dimension(n)] do
row:=[];
        for j in [1..Dimension(n-1)] do          # one pass per column j
        sum:=0;
                for x in Boundary(n,i) do          # ...rescans the whole sparse list each time
                if AbsoluteValue(x[1])=j then
                sum := sum + SignInt(x[1]);
                fi;
                od;
        row[j]:= RemInt(sum,2);
        od;
M1[i]:=row;
od;
```

`Boundary(n,i)` is a *sparse* list of `[target, …]` terms, but it is re-scanned for every column `j`, giving **O(Dimension(n-1) · |Boundary(n,i)|)** per row. Scatter each term straight into its target column to get **O(Dimension(n-1) + |Boundary(n,i)|)** per row (the `Dimension(n-1)` factor disappears):

```gap
for i in [1..Dimension(n)] do
    row := List([1..Dimension(n-1)], j -> 0);
    for x in Boundary(n,i) do
        row[AbsoluteValue(x[1])] := row[AbsoluteValue(x[1])] + SignInt(x[1]);
    od;
    M1[i] := List(row, s -> RemInt(s,2));
od;
```

Same dense, mod-2 result. The **M2** block (L193–205) has the identical pattern, and so does the "raw method" **M1** build at **L1703–1715** — all three sites share the fix. `CR_Mod2CocyclesAndCoboundaries` is called many times, so the saving compounds (though rank/solve steps still dominate overall runtime). *(Credit: raised by Gemini; distinct from P3, which concerns `Position` scans in `RelReduceLett`.)*

---

## 🟢 Maintainability / style

### M1 — Massive duplication of the degree blocks
The r = 2,3,4,5,6,7,8 sections of `Mod2RingGensAndRels` (≈ L716–3551, ~2,900 lines) are near-identical copy-paste: each repeats the "compute cup product → test coboundary → extend basis → build relation-reduction matrix → test reducibility" skeleton with only the degrees and source bases changed. This is the single largest risk factor for divergent bugs (a fix applied to one block but not the parallel ones). A degree-parameterized helper — e.g. `ProcessCupLevel(degU, degV, basesU, basesV, CB, …)` returning `(newBasis, newRelations)` — would collapse this to a few hundred lines and make the r=2/3 vs r≥4 inconsistencies (P4) impossible.

### M2 — Dead code and unused locals
- The whole generator-naming block is disabled by `if false then … fi;` (L4297–4361) and carries its only uses of `Rdim, lst1..lst6, lst1to4, GenName_standard`.
- Commented-out logic recurs in every degree block (`IO`, `NonNegativeVec`, `CupTemp`/`CupTempLett`, the alternate "raw" reinsertion). `CupTemp`/`CupTempLett` are still *assigned* `[]` but only *used* in comments.
- Unused declared locals include (non-exhaustive): in `CR_Mod2CocyclesAndCoboundaries` — `ColMat, InvColMat, RemoveRowsMat, InsertRowsList, Rels, Mod2Cohomologydim`; in `Letter2Monomial` — `v`; in `Mod2RingGenerators` — `sol, ln`; in `Mod2RingGensAndRels` — `CupBase2all, SmithRecord, cc, IO, Lett2`; in `SpaceGroupCohomologyRingGapInterface` — `o0, C2, C2p, M, P, C3`. GAP flags these; pruning them reduces noise.

### Smaller notes
- 🟢 **N1 —** `GF2ToZ` (L48–61) builds the 0/1 vector with an 8-line explicit loop (L51–58). It is correct and not wasteful (one pass, O(n)), but non-idiomatic; it reduces to `return List(v, IntFFE);` (`IntFFE(0*Z(2))=0`, `IntFFE(Z(2))=1`). Pure readability — no correctness or performance change. Both forms assume GF(2) field-element input, which every call site satisfies. (Aside: the current `v[k] = 0*Z(2)` test would map a *plain integer* `0` to `1`, since integer `0 ≠ 0*Z(2)` in GAP — never triggered in practice, but the `IntFFE` form sidesteps it.)
- 📚 `Letter2Monomial` hardcodes `["A","B","C","D","E","F"]` (L79–89): only 6 generator-degrees supported; degree-7+ would index out of range. Fine for space groups, but a silent cap.
- 📚 Comment/code mismatch at L461: comment says `QuoInt(j,2)` but code is `QuoInt(j+1,2)`.
- 📚 Mixed return types from `CR_Mod2CocyclesAndCoboundaries`: `false` (n<0), `[0]` (n=0), or a `rec` otherwise — callers must know which.
- 📚 `PrintMonomialString` returns `0` (L135); `PointGroupTranslationExtension`/`IrreducibleWyckoffPoints` mix `Print` side effects with return values. Consider HAP's `Info` levels instead of bare `Print` so verbosity is controllable.
- 📚 Inconsistent indentation in the M1/M2 construction (L172–214) and elsewhere hurts readability.
- 📚 Functions are defined by plain global assignment (`Name := function … end;`) rather than `BindGlobal`/`DeclareGlobalFunction`, so re-reading the file silently clobbers globals and there is no read-only protection. Acceptable for a script, worth noting for packaging.
- 📚 The file is one 4,519-line unit. Splitting into (cocycle/coboundary core) · (ring generators & relations) · (Wyckoff / LSM interface) would aid navigation and testing.

---

## Suggested priority order

1. **B1, B2** — fix the two latent crashes (or delete the unreachable standalone path and add the missing `else`/guard).
2. **B3, D1, B4** — remove the dead `toggle` argument, fix the misleading header doc, declare the captured loop variables local.
3. **P2, P1, P4, P7** — the cheap, high-value performance fixes (hoist `Cohomology`/zero-vector, hoist `Concatenation`, reuse the chain map at r=2/3, scatter the boundary matrices instead of rescanning).
4. **F2, F1, F3, F4** — document/guard the external-global and hardcoded-resolution assumptions.
5. **M1, M2** — refactor the duplicated degree blocks and prune dead code (largest long-term payoff, but only after the above are locked in by tests).

---

## What I checked vs. did not

- **Checked (static):** control flow and variable scoping of all top-level functions and their nested helpers; argument handling; data dependencies; the cocycle/coboundary/cup-product math at a structural level (the GF(2) nullspace/image construction in `CR_Mod2CocyclesAndCoboundaries` is sound, including the empty-/top-degree edge cases).
- **Verified empirically:** GAP's unassigned-local error semantics underlying B1 (isolated snippet).
- **Not done:** full end-to-end execution (requires HAP + the large data file); numerical confirmation of the hardcoded generator vectors / H⁶ table (F1); confirmation that `Size(R)` returns the resolution dimension list in your HAP version (used by `Length(Size(R))` at L688, L1699, L1836, L2555, L3639 — worth a one-line sanity check on your setup).
