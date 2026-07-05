# homotopy.gi — contracting homotopies for space-group resolutions.
# Read after linalg.gi and BEFORE functions.gi (see read.g).
#
# HAP's ResolutionSpaceGroup (Polymake-based) returns homotopy:=fail, but
# CR_ChainMapFromCocycle (hap/lib/Rings/cocycleChainMap.gi) needs R!.homotopy
# for every cup product.  The package-facing constructor here,
# SGC_ResolutionSpaceGroup, rebuilds the same resolution with a genuine
# integral contracting homotopy installed at construction time:
#   1. SGC_CrystallographicComplexWithGeometry — HAP's CrystallographicComplex
#      with the fundamental-domain geometry retained;
#   2. SGC_AttachCellularContraction — a Z-linear cellular contraction of the
#      Voronoi tessellation (distance-shelling discrete vector field);
#   3. SGC_FreeGResolutionWithHomotopy — HAP's FreeGResolution with the
#      total perturbation homotopy grafted in (FreeZGResolution-style).
# GF(2) does not enter the homotopy construction; it only appears later in
# the cohomology pipeline (functions.gi).
#
# The rest of this header documents the older GF(2) search-based fallback
# SGC_AttachHomotopyMod2 (kept for debugging via
# SGC_ResolutionSpaceGroupSearchHomotopy): h(letter) is found by recursively
# reducing to a boundary-filling problem d(w)=z and solving it over GF(2) on
# a growing finite patch of the orbit space.
#
# Adapted from attach_hap_homotopy.g (repo root) with these deliberate changes:
#  * All word arithmetic is over GF(2) (see the original rationale below).
#  * The group-translation convention is fixed to HAP's left-module word
#    convention (see the original rationale below).
#  * Performance (this revision): profiling SpaceGroupCohomologyRingGapInterface
#    on IT=89 showed the naive version ~50x slower than the old
#    ResolutionAlmostCrystalGroup path, almost entirely inside this file.
#    Three structural fixes, all local to this file:
#     - EnsureElt used Position(R!.elts, g), an O(n) scan repeated millions
#       of times against a list that keeps growing. Replaced with a
#       dictionary (hash) kept in sync with R!.elts incrementally -- catching
#       up on any elements OTHER code appended directly to R!.elts (HAP's own
#       CR_ChainMapFromCocycle does this with its own local Mult/Position) so
#       the cache can never go stale, at the cost of only rescanning the
#       delta since the last sync, not the whole list.
#     - UpLettersHittingLowLetter rescanned every (k+1)-generator's boundary
#       from scratch on every call. The incidence "which (k+1)-generators
#       have generator i in their boundary, and via which base element" is a
#       fixed combinatorial fact about the resolution; precomputed once per
#       degree (GetUpIncidence) and cached, each call becomes a table lookup.
#     - PartialHomotopy rebuilt a dense GF(2) linear system from scratch at
#       every radius-growth step and called SolutionMat on the whole thing.
#       Replaced with an incremental sparse XOR-basis (a standard GF(2)
#       linear-basis technique): each row is a mod-2 word (already this
#       file's native representation), keyed by its lexicographically-largest
#       letter as "pivot". New rows (from newly discovered up-letters) are
#       folded into the existing basis with InsertRow instead of re-deriving
#       everything; membership/fill-extraction (ReduceTarget) is a single
#       pass over the existing basis. This also only re-explores the newly
#       discovered frontier each step instead of the whole accumulated
#       low-basis (UpLettersHittingLowLetter is a pure function of a low
#       letter, so re-scanning previously-seen ones can never find anything
#       new).
#
# A mod-2 word here is a sorted list of letters [i,g], i>0, each of odd
# multiplicity — i.e. the support of the chain over GF(2). These are valid
# HAP words, so they can be returned to HAP code directly.

SGC_HOMOTOPY_MAXRADIUS := 24;

#####################################################################
SGC_AttachHomotopyMod2:=function(R)
local maxRadius, cache, idpos,
      eltDict, lastSynced, upIncidence,
      NormalizeWord, SymDiff, EnsureElt, Mul,
      BoundaryLetter, BoundaryWord, GetUpIncidence, UpLettersHittingLowLetter,
      InsertRow, ReduceTarget, PartialHomotopy, HomotopyWord,
      CacheLookup, CacheStore, HLetter;

maxRadius := SGC_HOMOTOPY_MAXRADIUS;
cache := [];
upIncidence := [];

#####################################################################
NormalizeWord:=function(word)     # HAP word -> mod-2 word (sorted set of [i,g])
local s, t, key;
s := [];
for t in word do
    key := [AbsInt(t[1]), t[2]];
    if key in s then
        RemoveSet(s, key);
    else
        AddSet(s, key);
    fi;
od;
return s;
end;
#####################################################################
SymDiff:=function(a, b)           # sum of two mod-2 words
return NormalizeWord(Concatenation(a, b));
end;
#####################################################################
# Dictionary mirror of R!.elts, kept in sync lazily: other code (e.g. HAP's
# own CR_ChainMapFromCocycle, which has its own local Position/Add on
# R!.elts) can append elements behind our back, so every lookup first catches
# the dictionary up on whatever was appended since the last sync -- only the
# delta is rescanned, never the whole list.
eltDict := NewDictionary(R!.elts[1], true);
lastSynced := 0;
EnsureElt:=function(g)
local pos, i;
if Length(R!.elts) > lastSynced then
    for i in [lastSynced+1..Length(R!.elts)] do
        if LookupDictionary(eltDict, R!.elts[i]) = fail then
            AddDictionary(eltDict, R!.elts[i], i);
        fi;
    od;
    lastSynced := Length(R!.elts);
fi;
pos := LookupDictionary(eltDict, g);
if pos <> fail then return pos; fi;
if IsBound(R!.appendToElts) then
    R!.appendToElts(g);
else
    Add(R!.elts, g);
fi;
pos := Length(R!.elts);
if R!.elts[pos] <> g then
    # appendToElts normalized to something else; fall back to a genuine scan
    # rather than risk a wrong index.
    pos := Position(R!.elts, g);
    if pos = fail then
        Error("SGC_AttachHomotopyMod2: could not add group element to R!.elts");
    fi;
fi;
AddDictionary(eltDict, g, pos);
lastSynced := Length(R!.elts);
return pos;
end;
#####################################################################
Mul:=function(gpos, hpos)         # position of Elts[gpos]*Elts[hpos]
return EnsureElt(R!.elts[gpos] * R!.elts[hpos]);
end;
#####################################################################

idpos := EnsureElt(One(R!.group));

#####################################################################
BoundaryLetter:=function(k, letter)   # boundary of translated letter, mod 2
local out, term;
if k <= 0 then return []; fi;
out := [];
for term in R!.boundary(k, AbsInt(letter[1])) do
    Add(out, [AbsInt(term[1]), Mul(letter[2], term[2])]);
od;
return NormalizeWord(out);
end;
#####################################################################
BoundaryWord:=function(k, word)
local out, t;
out := [];
for t in word do
    out := SymDiff(out, BoundaryLetter(k, t));
od;
return out;
end;
#####################################################################
# Precompute, once per degree and cached, which (k+1)-generators have a
# given k-generator i in their boundary and via which base (untranslated)
# element -- a fixed combinatorial fact about the resolution. Replaces a
# full rescan of R!.boundary(k+1, *) on every UpLettersHittingLowLetter call.
GetUpIncidence:=function(k)
local tbl, j, term, gen;
if IsBound(upIncidence[k+1]) then return upIncidence[k+1]; fi;
tbl := List([1..R!.dimension(k)], x -> []);
if R!.dimension(k+1) > 0 then
    for j in [1..R!.dimension(k+1)] do
        for term in R!.boundary(k+1, j) do
            gen := AbsInt(term[1]);
            Add(tbl[gen], [j, term[2]]);
        od;
    od;
fi;
upIncidence[k+1] := tbl;
return tbl;
end;
#####################################################################
UpLettersHittingLowLetter:=function(k, low)
# all (k+1)-letters whose boundary can contain the k-letter low
local ans, tbl, pr, upg;
ans := [];
if R!.dimension(k+1) = 0 then return ans; fi;
tbl := GetUpIncidence(k);
for pr in tbl[low[1]] do
    # up_g * term_g = low_g  (left convention)
    upg := R!.elts[low[2]] * R!.elts[pr[2]]^-1;
    AddSet(ans, [pr[1], EnsureElt(upg)]);
od;
return ans;
end;
#####################################################################
# Incremental GF(2) linear basis over mod-2 words, keyed by each row's
# lexicographically-largest letter as "pivot" (a standard sparse XOR-basis
# technique). InsertRow folds one new row (a boundary word, tagged with which
# single up-letter produced it) into the basis in place; ReduceTarget checks
# whether a word is in the current span and, if so, returns the combination
# of up-letters (as a mod-2 word) whose boundaries sum to it. Both take the
# basis (a dictionary) as an explicit argument so a fresh one can be started
# per PartialHomotopy call while growing it incrementally across radius steps
# within that call -- no full re-derivation when the search patch grows.
InsertRow:=function(pivots, rowWord, comboWord)
local piv, existing;
while rowWord <> [] do
    piv := rowWord[Length(rowWord)];
    if KnowsDictionary(pivots, piv) then
        existing := LookupDictionary(pivots, piv);
        rowWord := SymDiff(rowWord, existing[1]);
        comboWord := SymDiff(comboWord, existing[2]);
    else
        AddDictionary(pivots, piv, [rowWord, comboWord]);
        return;
    fi;
od;
# rowWord reduced to the zero word: it was already in the span: nothing to add.
end;
#####################################################################
ReduceTarget:=function(pivots, word)
local piv, existing, combo;
combo := [];
while word <> [] do
    piv := word[Length(word)];
    if not KnowsDictionary(pivots, piv) then
        return fail;
    fi;
    existing := LookupDictionary(pivots, piv);
    word := SymDiff(word, existing[1]);
    combo := SymDiff(combo, existing[2]);
od;
return combo;
end;
#####################################################################
PartialHomotopy:=function(k, z)
# given a mod-2 k-cycle z, find w in degree k+1 with d(w)=z (mod 2),
# growing a finite patch until the GF(2) basis spans z.
local zword, pivots, seenUp, seenLow, frontier, newFrontier, step,
      low, low2, u, bdy, upword, res;

zword := NormalizeWord(z);
if zword = [] then return []; fi;
if R!.dimension(k+1) = 0 then
    Error("SGC_AttachHomotopyMod2: no term R_", k+1,
          " available; build the resolution one degree higher");
fi;

pivots := NewDictionary([1,1], true);
seenUp := NewDictionary([1,1], true);
seenLow := NewDictionary([1,1], true);
for low in zword do AddDictionary(seenLow, low, true); od;
frontier := ShallowCopy(zword);

res := ReduceTarget(pivots, zword);
if res <> fail then return res; fi;

for step in [1..maxRadius] do
    newFrontier := [];
    for low in frontier do
        for u in UpLettersHittingLowLetter(k, low) do
            if not KnowsDictionary(seenUp, u) then
                AddDictionary(seenUp, u, true);
                bdy := BoundaryLetter(k+1, u);
                upword := NormalizeWord([[u[1], u[2]]]);
                InsertRow(pivots, bdy, upword);
                for low2 in bdy do
                    if not KnowsDictionary(seenLow, low2) then
                        AddDictionary(seenLow, low2, true);
                        Add(newFrontier, low2);
                    fi;
                od;
            fi;
        od;
    od;

    res := ReduceTarget(pivots, zword);
    if res <> fail then return res; fi;

    if newFrontier = [] then
        break;   # patch saturated (no new letters reachable): won't improve by looping further
    fi;
    frontier := newFrontier;
od;

Error("SGC_AttachHomotopyMod2: no filling found after maxRadius = ", maxRadius,
      "; raise SGC_HOMOTOPY_MAXRADIUS and retry");
end;
#####################################################################
HomotopyWord:=function(k, word)
local out, t;
out := [];
for t in word do
    out := SymDiff(out, HLetter(k, t));
od;
return out;
end;
#####################################################################
CacheLookup:=function(k, letter)
local slot, pr;
if not IsBound(cache[k+1]) then return fail; fi;
slot := cache[k+1];
if not IsBound(slot[letter[1]]) then return fail; fi;
for pr in slot[letter[1]] do
    if pr[1] = letter[2] then return ShallowCopy(pr[2]); fi;
od;
return fail;
end;
#####################################################################
CacheStore:=function(k, letter, value)
if not IsBound(cache[k+1]) then cache[k+1] := []; fi;
if not IsBound(cache[k+1][letter[1]]) then cache[k+1][letter[1]] := []; fi;
Add(cache[k+1][letter[1]], [letter[2], value]);
end;
#####################################################################
HLetter:=function(k, letter)      # h of one letter; sign irrelevant mod 2
local key, hit, dx, hdx, cycle, fill;

if k < 0 then
    Error("SGC_AttachHomotopyMod2: R!.homotopy is only defined for k >= 0");
fi;
key := [AbsInt(letter[1]), letter[2]];

hit := CacheLookup(k, key);
if hit <> fail then return hit; fi;

if k = 0 then
    # augmentation convention: eps([i,g])=1 and h_{-1}(1)=[1,1_G]
    cycle := SymDiff([key], [[1, idpos]]);
else
    dx := BoundaryLetter(k, key);
    hdx := HomotopyWord(k-1, dx);
    cycle := SymDiff([key], hdx);
fi;

fill := PartialHomotopy(k, cycle);
CacheStore(k, key, fill);
return ShallowCopy(fill);
end;
#####################################################################

R!.homotopy := function(k, letter)
    return HLetter(k, letter);
end;

R!.homotopyData := rec(
    idpos := idpos,
    cache := cache,
    boundaryWord := BoundaryWord,
    normalizeWord := NormalizeWord,
    symDiff := SymDiff
);

return R;
end;
#####################################################################

#####################################################################
SGC_HomotopyCheckMod2:=function(R, k, letter)
# verify (dh + hd)(letter) = letter - "basepoint" over GF(2)
local BoundaryWord, NormalizeWord, SymDiff, HomotopyWord, lhs, rhs, idpos;

if not IsBound(R!.homotopyData) then
    Error("SGC_HomotopyCheckMod2: run SGC_AttachHomotopyMod2 on R first");
fi;
BoundaryWord := R!.homotopyData.boundaryWord;
NormalizeWord := R!.homotopyData.normalizeWord;
SymDiff := R!.homotopyData.symDiff;
idpos := R!.homotopyData.idpos;

HomotopyWord:=function(j, word)
local out, t;
out := [];
for t in word do
    out := SymDiff(out, R!.homotopy(j, t));
od;
return out;
end;

if k = 0 then
    lhs := BoundaryWord(1, R!.homotopy(0, letter));
    rhs := SymDiff(NormalizeWord([letter]), [[1, idpos]]);
else
    lhs := SymDiff(BoundaryWord(k+1, R!.homotopy(k, letter)),
                   HomotopyWord(k-1, BoundaryWord(k, [letter])));
    rhs := NormalizeWord([letter]);
fi;

return lhs = rhs;
end;
#####################################################################

#####################################################################
SGC_CrystallographicComplexWithGeometry:=function(arg)
local G,K,pt,AG,SAG,F,V,FL,Boundaries,ind,inv,
Cells,VCells,GCells,PGCells,setV,
Dimn,BNDS, Bndy, REPS, NREPS,OREPS, Elts,Canonical, CanonicalTable,
STAB, STABREC, ACTION, Action, StandardWord,
PH, PGREPS,
n,Y,k,x,y,i;

# This is HAP's CrystallographicComplex constructor with the local cell
# geometry retained.  Keep the public non-free resolution data identical to
# HAP; the extra record is used only to build a cellular contraction.
G:=arg[1];
G:=StandardAffineCrystGroup(G);
if Length(arg)>1 then pt:=arg[2]; else
pt:=List([1..Length(One(G))-1], i->1/Primes[10+i]);
fi;

AG:=AffineCrystGroupOnRight(GeneratorsOfGroup(G));;
SAG:=StandardAffineCrystGroup(AG);;

# Data for the Voronoi-shelling contraction (SGC_AttachCellularContraction):
# the tessellation below is the Dirichlet/Voronoi diagram of the orbit
# SAG(pt), so tiles are indexed by group elements (pt is generic).  Site
# enumeration needs the translation lattice (must be the standard Z^d for a
# standard group) and one affine representative per point-group element.
if not TranslationBasis(SAG) = IdentityMat(Length(pt)) then
    Error("SGC_CrystallographicComplexWithGeometry: translation lattice of ",
          "the standardized group is not the standard lattice");
fi;
PH:=PointHomomorphism(SAG);
PGREPS:=List(Elements(Image(PH)), p->PreImagesRepresentative(PH,p));

F:=FundamentalDomainStandardSpaceGroup(pt,SAG);;
V:=Polymake(F,"VERTICES");;
setV:=Set(V);
FL:=PolymakeFaceLattice(F,true);
FL:=FL{[1..Length(FL)-1]};
ind:=List(FL,h->Minimum(Flat(h)));

for n in [1..Length(ind)] do
FL[n]:=FL[n]-ind[n]+1;
od;
ind:=List(FL,h->Length(h));

Boundaries:=[];
Boundaries[1]:=List([1..ind[1]],x->[1,0]);
for n in [2..Length(ind)] do
Boundaries[n]:=List([1..ind[n]],x->[]);
for x in [1..Length(FL[n-1])] do
for i in FL[n-1][x] do
Add(Boundaries[n][i],x);
od;
od;
Boundaries[n]:=List(Boundaries[n],b->Concatenation([Length(b)],SSortedList(b)));
od;
Boundaries[n+1]:=[Concatenation([ind[n]],[1..ind[n]])];
Boundaries[n+2]:=[];

Y:=RegularCWComplex(Boundaries);

Cells:=[];
Cells[1]:=List([1..Y!.nrCells(0)],i->[i]);
for k in [1..Dimension(Y)] do
Cells[k+1]:=List(Y!.boundaries[k+1], x->
SSortedList(  Concatenation( List(x{[2..Length(x)]}, i->Cells[k][i])  ))
);
od;

VCells:=[];
for k in [1..Length(Cells)] do
VCells[k]:= List(Cells[k], S -> SortedList(List(S,i->V[i])));
od;

inv:=function(X)
local st;
st:=Set(List(X,v->V[v]));
st:=OrbitPartInVertexSetsStandardSpaceGroup(SAG,st,setV);
return st;
end;

GCells:=[];
for k in [1..Length(Cells)] do
GCells[k]:=Classify(  Cells[k] , inv);
od;

PGCells:=[];
for k in [1..Length(GCells)] do
PGCells[k]:=[];
for x in GCells[k] do
y:=List(x, i->Position(Cells[k],i));
Add(PGCells[k],y);
od;
od;

OREPS:=[];
for k in [1..Length(PGCells)] do
OREPS[k]:=[];
for x in PGCells[k] do
for i in x do
OREPS[k][i]:=x[1];
od;
od;
od;

NREPS:=[];
for k in [1..Length(Cells)] do
REPS:=List(GCells[k],x->x[1]);
NREPS[k]:=List(REPS,x->Position(Cells[k],x));
od;

Elts:=[One(SAG)];

Canonical:=function(n,k)
local kk, ll, g, pos;
kk:=OREPS[n+1][k];
ll:=Position(NREPS[n+1],kk);
g:=RepresentativeActionOnRightOnSets(SAG,VCells[n+1][kk], VCells[n+1][k]);
if g=fail then Print("Error!!!\n"); return fail;fi;
pos:=Position(Elts,g);
if pos=fail then Add(Elts,g); pos:=Length(Elts); fi;
return [ll,pos];
end;

Dimn:=function(n);
if n<0 or n>Dimension(Y) then return 0; fi;
return Length(GCells[n+1]);
end;

BNDS:=[];
for n in [2..Length(NREPS)] do
BNDS[n-1]:=[];
for k in NREPS[n] do
x:=1*Y!.boundaries[n][k]{[2..Y!.boundaries[n][k][1]+1]};
y:=List(x, i->Canonical(n-2,i));
for i in [1..Length(Y!.orientation[n][k])] do
y[i]:=[Y!.orientation[n][k][i]*y[i][1],y[i][2]];
od;
Add(BNDS[n-1], y );
od;
od;

CanonicalTable:=[];
for n in [0..Dimension(Y)] do
CanonicalTable[n+1]:=List([1..Y!.nrCells(n)], k->Canonical(n,k));
od;

Bndy:=function(n,k);
if n<1 or n>Dimension(Y) then return []; fi;
return BNDS[n][k];
end;

STAB:=function(n,kk)
local st,k;
k:=NREPS[n+1][kk];
st:=VCells[n+1][k];
st:=StabilizerOnSetsStandardSpaceGroup(SAG,st);
st:=List(st,a->TransposedMat(a));
st:=Group(st);
return st;
end;

STABREC:=[];
for n in [1..Length(NREPS)] do
STABREC[n]:=[];
for k in [1..Length(NREPS[n])] do
STABREC[n][k]:=STAB(n-1,k);
od;
od;

STAB:=function(n,k);
return STABREC[n+1][k];
end;

StandardWord:=function(k,bnd)
local w,x,c,pos;
w:=[];
for x in bnd do
c:=CanonicalRightCosetElement(STAB(k,AbsInt(x[1])), Elts[x[2]]^-1 )^-1;
pos:=Position(Elts,c);
if pos=fail then Add(Elts,c); pos:=Length(Elts); fi;
Add(w,[x[1]*ACTION(k,AbsInt(x[1]),x[2]),pos]);
od;
return AlgebraicReduction(w);
end;

Action:=function(uu,V)
local L, M, T, d, u;
u:=TransposedMat(uu);
d:=Length(V);
L:=1*u{[1..d]};
M:=List(L,r->r{[1..d]});
T:=1*u[d+1]{[1..d]};
return V*M + T;
end;

ACTION:=function(n,k,h)
local bas, Gbas, mat,id,r,u,H,A,cg,d;

if n=0 then return 1; fi;
if Elts[h]=One(SAG) then return 1; fi;
H:=STAB(n,AbsInt(k));

id:=CanonicalRightCosetElement(H,Identity(SAG));
r:=CanonicalRightCosetElement(H,Elts[h]^-1);
r:=id^-1*r;
u:=r*Elts[h];

if not u in H then Print("ERROR\n"); fi;

k:=NREPS[n+1][AbsInt(k)];
A:=VCells[n+1][k];
cg:=Sum(A)/(Length(A));
A:=List(A,v->v-cg);
bas:=SemiEchelonMat(A).vectors;
Gbas:=List(bas,V->Action(u,V+cg)-cg);
mat:=List(Gbas,b->SolutionMat(bas,b));
d:=Determinant(mat);
if not AbsInt(d)=1 then Print("ERROR\n"); fi;
return Determinant(mat);

end;

Elts:=List(Elts,a->TransposedMat(a));
SAG:=GeneratorsOfGroup(SAG);
SAG:=List(SAG,a->TransposedMat(a));
SAG:=Group(SAG);

K:=    Objectify(HapNonFreeResolution,
           rec(
            dimension:=Dimn,
            boundary:=Bndy,
            homotopy:=fail,
            elts:=Elts,
            group:=SAG,
            stabilizer:=STAB,
            action:=ACTION,
            standardWord:=StandardWord,
            properties:=
             [["type","non-free resolution"],
              ["length",1000],
              ["characteristic", 0] ]));

K!.sgcGeometry := rec(
    regularCW := Y,
    vertices := V,
    cells := Cells,
    vertexCells := VCells,
    orbitReps := NREPS,
    orbitOfCell := OREPS,
    fundamentalDomain := F,
    canonicalTable := CanonicalTable,
    basePoint := pt,
    pointGroupReps := List(PGREPS, a->TransposedMat(a))
);

RecalculateIncidenceNumbers(K);
return K;
end;
#####################################################################

#####################################################################
# SGC_AttachCellularContraction(P)
#
# Attach an integral (Z-linear) contracting homotopy P!.homotopy to the
# non-free crystallographic complex built by
# SGC_CrystallographicComplexWithGeometry.  The convention is HAP's:
#   h(q, [i,g]) is a word in degree q+1 with  d h + h d = id  (q > 0)  and
#   d h([i,g]) = [i,g] - [1, 1_G]  in degree 0;  h(-1, *) = [[1, 1_G]].
#
# Geometry: the cells of P are the cells of the Voronoi/Dirichlet
# tessellation of R^d whose sites are the orbit points g(pt), g in G (the
# fundamental domain F is the Dirichlet cell of the generic point pt, so
# tiles correspond bijectively to group elements).  The contraction comes
# from a discrete vector field built from a distance shelling of the tiles:
#  * every tile g gets the key [ |g(pt)-pt|^2, g(pt) ], compared
#    lexicographically; the base tile is the identity, key [0, pt];
#  * a cell sigma lies in exactly the tiles whose sites realize the minimal
#    distance to sigma's barycenter (Voronoi property); sigma is owned by
#    the key-smallest such tile m(sigma);
#  * inside a tile m the owned ("new") cells are exactly those not lying in
#    any earlier tile; they are paired off by a greedy free-face collapse of
#    the fundamental polytope onto the subcomplex of old cells.  For the
#    base tile (no old cells) the collapse retains one protected vertex,
#    the first vertex orbit representative, so the augmentation basepoint
#    is HAP's [1, 1_G].
# Per-tile collapse sequences give acyclic matchings, and any V-path that
# leaves a tile strictly decreases the tile key, whose initial segments are
# finite (the orbit is a Delone set); hence the global matching is acyclic
# and proper and the homotopy recursion below terminates.
#
# The recursion is the standard vector-field contraction, in the same shape
# as HAP's cubical CrystGcomplex homotopy
# (hap/lib/ArithmeticGroups/crystGcomplex.gi):
#   h(sigma) = t*( [tau] - h(d tau - t*sigma) )   for a source sigma with
# matched coface tau containing sigma with incidence t = +-1; h = 0 on
# targets and on the critical base vertex.  All letters are normalized
# through P!.standardWord, which carries the orientation-character signs
# (P!.action) of cells with orientation-reversing stabilizer elements.
#####################################################################
SGC_AttachCellularContraction:=function(P)
local geom, Y, dimY, pt, pgReps, Elts, dvec,
      eltDict, lastSynced, EnsureElt, idpos,
      AffineAction, FloorRat,
      vertDicts, baryY, Cofaces,
      KeyOfMat, KeyLess, TilesOfPoint,
      tileSiteList, tileDataList, TileData,
      maskKeys, maskMatchings, CollapseMatching, baseVertexY,
      StandardLetter, BoundaryWord, DVF, HWord, HLetter, cache,
      n, k, j, x;

if not IsBound(P!.sgcGeometry) then
    Error("SGC_AttachCellularContraction: missing sgcGeometry; build P ",
          "with SGC_CrystallographicComplexWithGeometry");
fi;

geom := P!.sgcGeometry;
Y := geom.regularCW;
dimY := Dimension(Y);
pt := geom.basePoint;
pgReps := geom.pointGroupReps;
Elts := P!.elts;
dvec := [1..Length(pt)];
cache := List([0..dimY], q->[]);
tileSiteList := [];
tileDataList := [];
maskKeys := [];
maskMatchings := [];
baseVertexY := geom.orbitReps[1][1];

#####################################################################
FloorRat:=function(x)
if IsInt(x) then return x; fi;
if x > 0 then return Int(x); fi;
return Int(x) - 1;
end;
#####################################################################
AffineAction:=function(g,v)   # g in the P!.elts convention (transposed)
local u, d, M, T;
u := TransposedMat(g);
d := Length(v);
M := List(1*u{[1..d]}, r->r{[1..d]});
T := 1*u[d+1]{[1..d]};
return v*M + T;
end;
#####################################################################
# Dictionary mirror of P!.elts; other code (standardWord, HAP internals)
# appends to P!.elts directly, so each lookup first catches up on the delta.
eltDict := NewDictionary(P!.elts[1], true);
lastSynced := 0;
EnsureElt:=function(g)
local pos, i;
if Length(P!.elts) > lastSynced then
    for i in [lastSynced+1..Length(P!.elts)] do
        if LookupDictionary(eltDict, P!.elts[i]) = fail then
            AddDictionary(eltDict, P!.elts[i], i);
        fi;
    od;
    lastSynced := Length(P!.elts);
fi;
pos := LookupDictionary(eltDict, g);
if pos <> fail then return pos; fi;
Add(P!.elts, g);
pos := Length(P!.elts);
AddDictionary(eltDict, g, pos);
lastSynced := Length(P!.elts);
return pos;
end;
#####################################################################

idpos := EnsureElt(One(P!.group));

# Per-dimension lookup tables for the fundamental-domain complex Y:
# vertex-coordinate set -> Y-cell index, barycenters, and cofaces.
vertDicts := [];
baryY := [];
for n in [0..dimY] do
    vertDicts[n+1] := NewDictionary(geom.vertexCells[n+1][1], true);
    baryY[n+1] := [];
    for k in [1..Y!.nrCells(n)] do
        AddDictionary(vertDicts[n+1], geom.vertexCells[n+1][k], k);
        baryY[n+1][k] := Sum(geom.vertexCells[n+1][k]) /
                         Length(geom.vertexCells[n+1][k]);
    od;
od;

Cofaces := [];
for n in [0..dimY-1] do
    Cofaces[n+1] := List([1..Y!.nrCells(n)], k->[]);
    for k in [1..Y!.nrCells(n+1)] do
        x := Y!.boundaries[n+2][k];
        for j in x{[2..x[1]+1]} do
            AddSet(Cofaces[n+1][j], k);
        od;
    od;
od;
Cofaces[dimY+1] := List([1..Y!.nrCells(dimY)], k->[]);

#####################################################################
KeyOfMat:=function(g)
local s;
s := AffineAction(g, pt);
return [(s-pt)*(s-pt), s];
end;
#####################################################################
KeyLess:=function(a,b)
if a[1] <> b[1] then return a[1] < b[1]; fi;
return a[2] < b[2];
end;
#####################################################################
# All tiles containing the point b: group elements g (in the P!.elts matrix
# convention) whose site g(pt) realizes the exact minimal distance to b.
# g = (translation by lam) o (point-group representative A); enumeration is
# exact because the translation lattice is the standard Z^d.
TilesOfPoint:=function(b)
local best, cand, out, A, c, f, lo, hi, i, lam, dd, box, T, d;
d := Length(b);
best := fail;
for A in pgReps do
    c := AffineAction(A, pt);
    f := b - c;
    for lam in Cartesian(List(f, x->[FloorRat(x), FloorRat(x)+1])) do
        dd := (f-lam)*(f-lam);
        if best = fail or dd < best then best := dd; fi;
    od;
od;
cand := [];
for A in pgReps do
    c := AffineAction(A, pt);
    f := b - c;
    box := [];
    for i in dvec do
        lo := FloorRat(f[i]) - 1;
        while (f[i]-lo)^2 <= best do lo := lo - 1; od;
        hi := FloorRat(f[i]) + 2;
        while (f[i]-hi)^2 <= best do hi := hi + 1; od;
        Add(box, [lo+1..hi-1]);
    od;
    for lam in Cartesian(box) do
        dd := (f-lam)*(f-lam);
        if dd <= best then
            if dd < best then best := dd; fi;
            Add(cand, [dd, A, lam]);
        fi;
    od;
od;
out := [];
for x in cand do
    if x[1] = best then
        T := IdentityMat(d+1);
        for i in dvec do T[i][d+1] := x[3][i]; od;
        Add(out, T * x[2]);
    fi;
od;
return out;
end;
#####################################################################
# Greedy free-face collapse of Y onto the subcomplex of old cells.
# mask[n+1][k] = true iff the Y-cell k of dimension n is old.  Returns
# match[n+1][k] for new cells: a positive integer k' (paired coface, source),
# or 0 (paired face or critical, h = 0 there).
CollapseMatching:=function(mask)
local match, alive, isbase, protected, changed, freecf, nn, kk;
isbase := ForAll(mask, layer->ForAll(layer, x->not x));
if isbase then
    protected := baseVertexY;
else
    protected := 0;
fi;
match := List([0..dimY], n->[]);
alive := List([0..dimY], n->List([1..Y!.nrCells(n)], k->not mask[n+1][k]));
changed := true;
while changed do
    changed := false;
    for nn in [0..dimY-1] do
        for kk in [1..Y!.nrCells(nn)] do
            if alive[nn+1][kk] and not (nn = 0 and kk = protected) then
                freecf := Filtered(Cofaces[nn+1][kk], j->alive[nn+2][j]);
                if Length(freecf) = 1 then
                    match[nn+1][kk] := freecf[1];
                    match[nn+2][freecf[1]] := 0;
                    alive[nn+1][kk] := false;
                    alive[nn+2][freecf[1]] := false;
                    changed := true;
                fi;
            fi;
        od;
    od;
od;
if isbase then
    alive[1][protected] := false;
    match[1][protected] := 0;
fi;
for nn in [0..dimY] do
    for kk in [1..Y!.nrCells(nn)] do
        if alive[nn+1][kk] then
            Error("SGC_AttachCellularContraction: tile collapse jammed at ",
                  "dimension ", nn, ", cell ", kk);
        fi;
    od;
od;
return match;
end;
#####################################################################
# Old/new pattern and collapse matching of one tile, cached by its site.
TileData:=function(m)
local site, pos, key, mask, nn, kk, b, tiles, t, isold, data, mpos, matching;
site := AffineAction(m, pt);
pos := Position(tileSiteList, site);
if pos <> fail then return tileDataList[pos]; fi;
key := KeyOfMat(m);
mask := [];
for nn in [0..dimY] do
    mask[nn+1] := [];
    for kk in [1..Y!.nrCells(nn)] do
        b := AffineAction(m, baryY[nn+1][kk]);
        isold := false;
        for t in TilesOfPoint(b) do
            if KeyLess(KeyOfMat(t), key) then isold := true; break; fi;
        od;
        mask[nn+1][kk] := isold;
    od;
od;
if key[1] > 0 and ForAll(mask, layer->ForAll(layer, x->not x)) then
    Error("SGC_AttachCellularContraction: non-base tile with no old cells");
fi;
mpos := Position(maskKeys, mask);
if mpos = fail then
    matching := CollapseMatching(mask);
    Add(maskKeys, mask);
    Add(maskMatchings, matching);
else
    matching := maskMatchings[mpos];
fi;
data := rec(key := key, mask := mask, matching := matching);
Add(tileSiteList, site);
Add(tileDataList, data);
return data;
end;
#####################################################################
StandardLetter:=function(q, letter)
local w;
w := P!.standardWord(q, [1*letter]);
if Length(w) <> 1 then
    Error("SGC_AttachCellularContraction: standardizing one cell gave ",
          Length(w), " letters");
fi;
return w[1];
end;
#####################################################################
BoundaryWord:=function(q, word)   # integral boundary of a P-word
local out, t, x;
if q <= 0 then return []; fi;
out := [];
for t in word do
    for x in P!.boundary(q, AbsInt(t[1])) do
        if t[1] > 0 then
            Add(out, [x[1], EnsureElt(Elts[t[2]]*Elts[x[2]])]);
        else
            Add(out, [-x[1], EnsureElt(Elts[t[2]]*Elts[x[2]])]);
        fi;
    od;
od;
return AlgebraicReduction(out);
end;
#####################################################################
# The discrete vector field: for the canonical positive letter [r,e] in
# degree q return fail (target or critical cell, h = 0) or the matched
# coface as a canonicalizable letter [ll, epos] in degree q+1.
DVF:=function(q, r, e)
local rep, verts, b, tiles, m, mk, t, tk, td, minv, lverts, kkY, pairk, can, v;
if q >= dimY then return fail; fi;
rep := geom.orbitReps[q+1][r];
verts := List(geom.vertexCells[q+1][rep], v->AffineAction(Elts[e], v));
b := Sum(verts)/Length(verts);
tiles := TilesOfPoint(b);
m := tiles[1]; mk := KeyOfMat(m);
for t in tiles{[2..Length(tiles)]} do
    tk := KeyOfMat(t);
    if KeyLess(tk, mk) then m := t; mk := tk; fi;
od;
td := TileData(m);
minv := Inverse(m);
lverts := SortedList(List(verts, v->AffineAction(minv, v)));
kkY := LookupDictionary(vertDicts[q+1], lverts);
if kkY = fail then
    Error("SGC_AttachCellularContraction: cell not found in its owner tile");
fi;
if td.mask[q+1][kkY] then
    Error("SGC_AttachCellularContraction: owner tile marks the cell old");
fi;
pairk := td.matching[q+1][kkY];
if not IsBound(td.matching[q+1][kkY]) or pairk = fail then
    Error("SGC_AttachCellularContraction: unmatched cell in tile collapse");
fi;
if pairk = 0 then return fail; fi;
can := geom.canonicalTable[q+2][pairk];   # Elts[can[2]] maps rep(can[1]) to pairk
return [can[1], EnsureElt(m * Elts[can[2]])];
end;
#####################################################################
HWord:=function(q, word)
local out, t;
out := [];
for t in word do
    Append(out, HLetter(q, t));
od;
return AlgebraicReduction(out);
end;
#####################################################################
HLetter:=function(q, letter)
local std, r, e, sgn, h, xx, raw, cx, pos, i, t, c, dword;

if q < 0 then return [[1, idpos]]; fi;
if q > dimY then return []; fi;

std := StandardLetter(q, letter);
r := AbsInt(std[1]); e := std[2]; sgn := SignInt(std[1]);

if IsBound(cache[q+1][r]) and IsBound(cache[q+1][r][e]) then
    h := cache[q+1][r][e];
else
    xx := DVF(q, r, e);
    if xx = fail then
        h := [];
    else
        # raw translated boundary of the matched coface; distinct raw
        # letters standardize to distinct cells, so std occurs exactly once
        raw := [];
        for cx in P!.boundary(q+1, xx[1]) do
            Add(raw, [cx[1], EnsureElt(Elts[xx[2]]*Elts[cx[2]])]);
        od;
        pos := fail;
        for i in [1..Length(raw)] do
            cx := StandardLetter(q, raw[i]);
            if AbsInt(cx[1]) = r and cx[2] = e then
                pos := i; t := SignInt(cx[1]);
                break;
            fi;
        od;
        if pos = fail then
            Error("SGC_AttachCellularContraction: matched coface misses ",
                  "the cell; q=", q, " r=", r, " e=", e);
        fi;
        c := raw{Difference([1..Length(raw)], [pos])};
        dword := [ [xx[1], xx[2]] ];
        Append(dword, NegateWord(HWord(q, c)));
        dword := AlgebraicReduction(dword);
        if t < 0 then dword := NegateWord(dword); fi;
        h := dword;
    fi;
    if not IsBound(cache[q+1][r]) then cache[q+1][r] := []; fi;
    cache[q+1][r][e] := h;
fi;

if sgn < 0 then return NegateWord(StructuralCopy(h)); fi;
return StructuralCopy(h);
end;
#####################################################################

P!.homotopy := function(q, letter)
    return HLetter(q, letter);
end;

P!.sgcHomotopyData := rec(
    cache := cache,
    idpos := idpos,
    hword := HWord,
    boundaryWord := BoundaryWord
);

return P;
end;
#####################################################################

#####################################################################
SGC_NonFreeHomotopyCheck:=function(P, q, letter)
# verify the integral identity d h + h d = id (q > 0), resp.
# d h = id - basepoint (q = 0), on the non-free complex; chains are compared
# after P!.standardWord normalization (stabilizer-orbit normal form).
local D, H, idpos, diff;
if not IsBound(P!.sgcHomotopyData) then
    Error("SGC_NonFreeHomotopyCheck: run SGC_AttachCellularContraction first");
fi;
D := P!.sgcHomotopyData.boundaryWord;
H := P!.sgcHomotopyData.hword;
idpos := P!.sgcHomotopyData.idpos;
if q = 0 then
    diff := Concatenation(D(1, H(0, [1*letter])),
                          NegateWord([1*letter]), [[1, idpos]]);
else
    diff := Concatenation(D(q+1, H(q, [1*letter])),
                          H(q-1, D(q, [1*letter])),
                          NegateWord([1*letter]));
fi;
return P!.standardWord(q, AlgebraicReduction(diff)) = [];
end;
#####################################################################

#####################################################################
# SGC_FreeGResolutionWithHomotopy(P, N)
#
# Fork of HAP's FreeGResolution (hap/lib/Perturbations/freeRes.gi, (C) 2008
# Graham Ellis) that additionally installs a contracting homotopy on the
# total free resolution.  The resolution construction (stabilizer
# resolutions, Del_k perturbation differentials, boundaries, generator
# bookkeeping) is kept identical to HAP's so the returned boundaries match
# FreeGResolution(P, N) exactly; the only changes are
#  * the coefficient characteristic is fixed to 0 (integral),
#  * the homotopy block at the end, following HAP's FreeZGResolution
#    (hap/lib/ArithmeticGroups/freeZGRes.gi, homotopy part by Bui Anh Tuan),
#    with memoization added.  It consumes P!.homotopy, so run
#    SGC_AttachCellularContraction(P) first.
#####################################################################
SGC_FreeGResolutionWithHomotopy:=function(arg)
local
	P,N,prime,
	Dimension, DimensionRecord, DimRecs, FiltDimRecs,
        BinGp,
	Boundary,
	BoundaryP,
	Pair2Quad, Pair2QuadRec,
      	Quad2Pair,Quad2PairRec,
	HtpyGen, HtpyWord,
	StabGrps,
	StabResls,
        ResolutionFG,
	Action,
        AlgRed,
	EltsG, G, Mult, MultRecord,
	DelGen, DelWord, DelGenRec,
	PseudoBoundary,FinalBoundary,
        FilteredLength, FilteredDimension, FilteredDimensionRecord,
        VertBnd, InducedHtpyList, DelListSum, VertHtpy,
        NegateListWord, AlgRed5, HomotopyGen, HomotopyGenRec, Homotopy, BaseW0,
	L,i,k,n,q,r,s,t,bool;

SetInfoLevel(InfoWarning,0);

P:=arg[1];
N:=arg[2];
prime:=0;   # integral: GF(2) must not enter the homotopy construction

if not IsBound(P!.homotopy) or P!.homotopy = fail then
    Error("SGC_FreeGResolutionWithHomotopy: P has no contracting homotopy; ",
          "run SGC_AttachCellularContraction(P) first");
fi;

N:=Minimum(EvaluateProperty(P,"length"),N);
G:=P!.group;
bool:=not IsComponentObjectRep(One(G));
bool:=IsHapSL2Subgroup(G) or IsHapSL2OSubgroup(G) or IsBound(G!.bianchiInteger);
EltsG:=P!.elts;
BoundaryP:=P!.boundary;

BinGp:=ContractibleGcomplex("SL(2,O-2)");
BinGp:=BinGp!.stabilizer(0,4);;
BinGp:=Image(RegularActionHomomorphism(BinGp));

#############################
ResolutionFG:=function(G,n)
local x, tmp, iso,iso1,iso2,iso3,res,Q, fn;

##Added Jan 2012
if IsBound(P!.resolutions) and HasName(G) then
x:=Position(P!.resolutions[2], Name(G));
if not x=fail then return P!.resolutions[1][x]; fi;
fi;
##

###
if IsBianchiAbelianGroup(G) and Order(G)=infinity then
res:=ResolutionAbelianBianchiSubgroup(G,n);
return res;
fi;
if IsAbelian(G) and not bool then   #NOT SURE WHY I NEED TO DO THIS
res:=ResolutionAbelianGroup_alt(G,n);

return res;
fi;
###

iso:=RegularActionHomomorphism(G);
Q:=Image(iso);

if Order(Image(iso))=24 then
if IdGroup(Image(iso))=[24,3] then
iso1:=IsomorphismGroups(Q,BinGp);
res:=ResolutionFiniteGroup(BinGp,n);
res!.group:=G;
res!.elts:=List(res!.elts,x->
PreImagesRepresentative(iso,PreImagesRepresentative(iso1,x)));
return res;
fi;
fi;

res:=ResolutionGenericGroup(Q,n);
res!.group:=G;
res!.elts:=List(res!.elts,x->PreImagesRepresentative(iso,x));
return res;
###

end;
#############################

##############################################
AlgRed:= function(ww)
local x,i,v,k,u,w, vv, vvv, vvvv, pos, pos2;

w:=ww;
        for x in w do
        if x[2]<0 then x[1]:=-x[1];x[2]:=-x[2];fi;
        od;

         v:=Collected(w);
         vv:=List(v,x->x[1]);
         vvv:=List(vv,x->[AbsInt(x[1]),x[2],x[3]]);
         vvv:=SSortedList(vvv);
         vvvv:=[];
         for x in vvv do
             pos:=Position(vv,x);
             pos2:=Position(vv,[-x[1],x[2],x[3]]);
             k:=0;
             if not pos=fail then k:=k+v[pos][2]; fi;
             if not pos2=fail then k:=k-v[pos2][2]; fi;
             if not k=0 then Append(vvvv,MultiplyWord(k,[ [x[1],x[2], x[3]] ]  )); fi;

         od;
         return vvvv;

end;
##############################################

##############################################
if IsBound(P!.action) then
Action:=P!.action;
else
Action:=function(k,j,g) return 1; end;
fi;
##############################################

MultRecord:=[];
################################################################
Mult:=function(g,h)
local pos;
if not IsBound(MultRecord[g]) then MultRecord[g]:=[]; fi;
if not IsBound(MultRecord[g][h]) then
    pos:= Position(EltsG,EltsG[g]*EltsG[h]);
    if pos=fail then Add(EltsG,EltsG[g]*EltsG[h]);
    MultRecord[g][h]:= Length(EltsG);
    else MultRecord[g][h]:= pos;
    fi;
fi;
return MultRecord[g][h];
end;
################################################################

StabGrps:= List([0..Length(P)],n->
           List([1..P!.dimension(n)], k->P!.stabilizer(n,k)));

StabResls:=[];
i:=N;
#################################
for L in StabGrps do
Add(StabResls,List(L,
g->ExtendScalars(ResolutionFG(g,i),G,EltsG))
);
i:=Maximum(0,AbsInt(i-1));
od;
#################################

DimRecs:=List([0..N],i->[]);

###################################################################
Dimension:=function(k)
local dim,i,R;
dim:=0;
for i in [0..k] do
DimRecs[k+1][i+1]:=[];
for R in StabResls[i+1] do
dim:=dim+R!.dimension(k-i);
Add(DimRecs[k+1][i+1],dim);
od;
od;
return dim;
end;

DimensionRecord:=List([0..N],Dimension);

Dimension:=function(k);
return DimensionRecord[k+1];
end;
###################################################################

###################################################################
Quad2PairRec:=[];
for q in [0..N] do
Quad2PairRec[q+1]:=[];
for r in [1..Length(StabGrps[q+1])] do
Quad2PairRec[q+1][r]:=[];
for s in [0..N-q] do
Quad2PairRec[q+1][r][s+1]:=[];
od;od;od;
###################################################################

###################################################################
Pair2Quad:=function(k,n)
local qq,q,r,s,t;
#The n-th generator in degree k of our final resolution is actually the
#t-th generator in degree s of the resolution of the r-th stabilizer group
#of the q-th chain module of the non-free resolution. We need the
#function f(k,n)=[q,r,s,t] .

for qq in [0..N] do
if n <= DimRecs[k+1][qq+1][Length(DimRecs[k+1][qq+1])] then q:=qq; break; fi;
od;

r:=PositionProperty(DimRecs[k+1][q+1],x->(n<=x));

s:=k-q;

if r-1>0 then
t:=n-DimRecs[k+1][q+1][r-1];
else
if q>=1 then
t:=n-DimRecs[k+1][q][Length(  DimRecs[k+1][q] )];;
else t:=n;
fi;
fi;

Quad2PairRec[q+1][r][s+1][t]:=[k,n];
return [q,r,s,t];

end;

Pair2QuadRec:=List([1..N+1],i->[]);
for k in [0..N] do
for n in [1..Dimension(k)] do
Pair2QuadRec[k+1][n]:=Pair2Quad(k,n);
od;
od;

##############
Pair2Quad:=function(k,n)
local a;
if n>0 then
return StructuralCopy(Pair2QuadRec[k+1][n]);
else
a:=StructuralCopy(Pair2QuadRec[k+1][-n]);
a[4]:=-a[4];
return a;
fi;
end;
##############

##############
Quad2Pair:=function(q,r,s,t)
local a,pr,pt;
if r>0 then pr:=r;pt:=t;
else
pr:=-r;pt:=-t;
fi;

if pt>0 then
return StructuralCopy(Quad2PairRec[q+1][pr][s+1][pt]);
else
a:=StructuralCopy(Quad2PairRec[q+1][pr][s+1][-pt]);
a[2]:=-a[2];
return a;
fi;
end;
##############
###################################################################

###################################################################
HtpyGen:=function(q,s,r,t,g)
local y,pr,pt;
#This applies the "vertical homotopy" to the free group generator [r,t,g]
#in "dimension" [q,s]. The output is an "r-word" in "dimension" [q,s+1].

if r>0 then pr:=r;pt:=t;
else
pr:=-r;pt:=-t;
fi;

y:=StructuralCopy(StabResls[q+1][pr]!.homotopy(s,[pt,g]));
Apply(y,x->[pr,x[1],x[2]]);
return y;
end;
###################################################################

###################################################################
HtpyWord:=function(q,s,w)
local h,z,x,y;
#This applies the "vertical homotopy" to the r-word w in "dimension"
#[q,s]. The output is an r-word in "dimension" [q,s+1].

h:=[];
for y in w do
x:=[Action(q,y[1],y[3])*y[1],y[2],y[3]];
z:=HtpyGen(q,s,x[1],x[2],x[3]);
z:=List(z,a->[Action(q,a[1],a[3])*a[1],a[2],a[3]]);
Append(h,z);
od;

return AlgRed(h);
end;
###################################################################

DelGenRec:=[];
for k in [1..N+1] do
DelGenRec[k]:=[];
for q in [1..N+1] do
DelGenRec[k][q]:=[];
for s in [1..N+1] do
DelGenRec[k][q][s]:=[];
for r in [1..P!.dimension(q-1)] do
DelGenRec[k][q][s][r]:=[];
od;
od;
od;
od;
###################################################################
DelGen:=function(k,q,s,r,t)
local y,pr,pt,i;
#For k=0,1,2 ... this is the equivariant homomorphism
#Del_k:A_{q,s} ---> A_{q-k,s+k-1} applied to a free r-generator [r,t]
#in dimension [q,s].


if r>0 then pr:=r;pt:=t;
else
pr:=-r;pt:=-t;
fi;

##############
if IsBound(DelGenRec[k+1][q+1][s+1][pr][AbsInt(pt)]) then

if pt>0 then return DelGenRec[k+1][q+1][s+1][pr][pt];
else
return List(DelGenRec[k+1][q+1][s+1][pr][-pt], a->[a[1],-a[2],a[3]]);
fi;

fi;
##############

if k=0 then
if s=0 then return [];
else
y:=List(StabResls[q+1][pr]!.boundary(s,pt),x->[Action(q,r,x[2])*x[1],x[2]]);
if pt>0 then
DelGenRec[k+1][q+1][s+1][pr][pt]:= AlgRed(List(y,x->[pr,x[1],x[2]]));
return DelGenRec[k+1][q+1][s+1][pr][pt];
else DelGenRec[k+1][q+1][s+1][pr][-pt]:=AlgRed(List(y,x->[pr,-x[1],x[2]]));
return AlgRed(List(y,x->[pr,x[1],x[2]]));
fi;
fi;
fi;

if k=1 then
if s=0 then
if q=0 then return [];
fi;
y:=BoundaryP(q,pr);
if pt>0 then
DelGenRec[k+1][q+1][s+1][pr][pt]:= AlgRed(List(y,x->[x[1],1,x[2]]));
return DelGenRec[k+1][q+1][s+1][pr][pt];
else
DelGenRec[k+1][q+1][s+1][pr][-pt]:= AlgRed(List(y,x->[x[1],1,x[2]]));
return List(y,x->[x[1],-1,x[2]]);
fi;

else

if pt>0 then
DelGenRec[k+1][q+1][s+1][pr][pt]:=
AlgRed(HtpyWord(q-1,s-1,DelWord(1,q,s-1,DelGen(0,q,s,pr,-pt)))) ;
return DelGenRec[k+1][q+1][s+1][pr][pt];
else
DelGenRec[k+1][q+1][s+1][pr][-pt]:=
AlgRed(HtpyWord(q-1,s-1,DelWord(1,q,s-1,DelGen(0,q,s,pr,pt)))) ;

return
List(DelGenRec[k+1][q+1][s+1][pr][-pt], a->[a[1],-a[2],a[3]]);
fi;

fi;
fi;

y:=[];
for i in [1..k] do
Append(y,
HtpyWord(q-k,s+k-2,DelWord(i,q-k+i,s+k-i-1,DelGen(k-i,q,s,pr,-pt)))
);
od;
y:=AlgRed(y);

if pt>0 then
DelGenRec[k+1][q+1][s+1][pr][pt]:=y;
else
DelGenRec[k+1][q+1][s+1][pr][-pt]:=List(y,a->[a[1],-a[2],a[3]]);
fi;

return y;
end;
###################################################################

###################################################################
DelWord:=function(k,q,s,w)
local y,x;
#For k=0,1,2 ... this is the equivariant homomorphism
#Del_k:A_{q,s} ---> A_{q-k,s+k-1} applied to an r-word [[r,t,g],...]
#in dimension [q,s].

y:=[];
for x in w do
Append(y,List(DelGen(k,q,s,x[1],x[2]),
a->[a[1],a[2],Mult(x[3],a[3])]));
od;

return y;  #Added Jan 2013. Speeds up the calculation in some(!!) examples.
return AlgRed(y);

end;
###################################################################

###################################################################
Boundary:=function(k,n)
local q,s,r,t,x,y,z,i;
y:=Pair2Quad(k,n); q:=y[1];s:=y[3];r:=y[2];t:=y[4];

y:=[];

for i in [0..k] do
if q>=i then

z:=DelGen(i,q,s,r,t);

Append(y,
List(z,x->[Quad2Pair(q-i,x[1],s+i-1,x[2])[2],x[3]])  );
else break;
fi;
od;

return AlgebraicReduction(y);
end;
###################################################################

PseudoBoundary:=[];
for n in [1..N+1] do
PseudoBoundary[n]:=[];
od;
#######################################
FinalBoundary:=function(n,k)
local  pk;
pk:=AbsInt(k);
if  not IsBound(PseudoBoundary[n+1][pk]) then
PseudoBoundary[n+1][pk]:= Boundary(n,pk);
fi;
if k>0 then  return PseudoBoundary[n+1][k];
else  return NegateWord(PseudoBoundary[n+1][pk]); fi;
end;
#######################################



################spectral sequence requirements##################


FiltDimRecs:=[];
for k in [0..N] do
FiltDimRecs[k+1]:=[];
for i in [1..Dimension(k)] do
FiltDimRecs[k+1][i]:=Pair2Quad(k,i)[1];
od;
od;

FilteredLength:=Maximum(Flat(FiltDimRecs));

##################################################
FilteredDimension:=function(r,k);

return Length(Filtered(FiltDimRecs[k+1],x->x<=r));

end;
##################################################

###################################################################
# Contracting homotopy on the total resolution, following HAP's
# FreeZGResolution (Bui Anh Tuan).  Internal words are lists of 5-letters
# [q,s,r,t,g]: the (r,t)-generator of the stabilizer resolution in
# bidimension (q,s), translated by EltsG[g].  A sign on r or on t both mean
# negation of the letter (HAP's pr/pt convention), so the canonical form
# used for reduction and caching keeps r,t positive and carries
# SignInt(r)*SignInt(t).

##############################################
NegateListWord:=function(w)
return List(w, x->[x[1],x[2],-x[3],x[4],x[5]]);
end;
##############################################

##############################################
AlgRed5:=function(w)
local keys, coeffs, x, key, sgn, pos, out, i;
keys := [];
coeffs := [];
for x in w do
    key := [x[1], x[2], AbsInt(x[3]), AbsInt(x[4]), x[5]];
    sgn := SignInt(x[3])*SignInt(x[4]);
    pos := Position(keys, key);
    if pos = fail then
        Add(keys, key);
        Add(coeffs, sgn);
    else
        coeffs[pos] := coeffs[pos] + sgn;
    fi;
od;
out := [];
for pos in [1..Length(keys)] do
    if coeffs[pos] <> 0 then
        key := keys[pos];
        x := [key[1], key[2], key[3], key[4], key[5]];
        if coeffs[pos] < 0 then x[3] := -x[3]; fi;
        for i in [1..AbsInt(coeffs[pos])] do
            Add(out, ShallowCopy(x));
        od;
    fi;
od;
return out;
end;
##############################################

##############################################
VertBnd:=function(w)
local h,x,d;
#The vertical differential Del_0, applied to a list of 5-letters.
h:=[];
for x in w do
d:=StructuralCopy(DelGen(0,x[1],x[2],x[3],x[4]));
Apply(d,v->[x[1],x[2]-1,v[1],v[2],Mult(x[5],v[3])]);
Append(h,d);
od;
return h;
end;
##############################################

##############################################
InducedHtpyList:=function(w)
local h,y,z,pr,sgn;
#The "induced homotopy" h1: A_{q,0} -> A_{q+1,0} obtained from the
#contracting homotopy of the non-free complex P (only used at s=0):
#project to P, apply P!.homotopy, then include back into the s=0 row.
#The section of A_{q+1,0} -> P_{q+1} implied by the vertical reduction is
#i o p = 1 - Del_0 o h_vert (no h_vert o Del_0 term at s=0), so the raw
#inclusion [cell,1,elt] is normalized by that projection rather than by any
#coset/orientation convention of P!.standardWord -- the two canonical coset
#choices in HAP (standardWord vs ExtendScalars) do not agree.
h:=[];
for y in w do
pr:=AbsInt(y[3]);
sgn:=SignInt(y[3])*SignInt(y[4]);
z:=P!.homotopy(y[1],[pr,y[5]]);
z:=List(z,a->[y[1]+1,y[2],sgn*a[1],1,a[2]]);
Append(h,z);
od;
h:=AlgRed5(h);
return AlgRed5(Concatenation(h,NegateListWord(VertBnd(VertHtpy(h)))));
end;
##############################################

##############################################
DelListSum:=function(w)
local y,d,x,h,k;
#The perturbation part d(+) = Del_1 + ... + Del_q of the total
#differential, applied to a list of 5-letters.
h:=[];
for x in w do
y:=[];
for k in [1..x[1]] do
d:=StructuralCopy(DelGen(k,x[1],x[2],x[3],x[4]));
Apply(d,v->[x[1]-k,x[2]+k-1,v[1],v[2],Mult(x[5],v[3])]);
Append(y,d);
od;
Append(h,y);
od;
return h;
end;
##############################################

##############################################
VertHtpy:=function(w)
local h,x,y,v;
#The vertical homotopy h0, applied to a list of 5-letters.
h:=[];
for x in w do
v:=[x[1],x[2],Action(x[1],x[3],x[5])*x[3],x[4],x[5]];
y:=StructuralCopy(HtpyGen(v[1],v[2],v[3],v[4],v[5]));
Apply(y,a->[x[1],x[2]+1,Action(x[1],a[1],a[3])*a[1],a[2],a[3]]);
Append(h,y);
od;
return h;
end;
##############################################

HomotopyGenRec:=List([1..N+1],q->[]);
##############################################
HomotopyGen:=function(q,s,r,t,g)
local pr,pt,sgn,slot,y,h0,h1,h0dh1,e3,e,h2,v,x,ans;
#Total homotopy on the 5-letter [q,s,r,t,g]; output is a list of 5-letters
#in total degree q+s+1.  Same recursion as HAP's FreeZGResolution:
#  H = h0 - H(d+ h0) + [s=0]( h1 - h0(d+ h1) + H(d+ h0 (d+ h1)) )
#with h0 the vertical homotopy and h1 the induced homotopy.  Terminates
#because every d+ strictly lowers q.  Memoized on (q,s,|r|,|t|,g).

pr:=AbsInt(r); pt:=AbsInt(t);
sgn:=SignInt(r)*SignInt(t);

if not IsBound(HomotopyGenRec[q+1][s+1]) then HomotopyGenRec[q+1][s+1]:=[]; fi;
slot:=HomotopyGenRec[q+1][s+1];
if not IsBound(slot[pr]) then slot[pr]:=[]; fi;
if not IsBound(slot[pr][pt]) then slot[pr][pt]:=[]; fi;

if not IsBound(slot[pr][pt][g]) then

    y:=[q,s,pr,pt,g];

    if s=0 then
        h1:=InducedHtpyList([y]);
        h0dh1:=VertHtpy(DelListSum(h1));
        v:=DelListSum(h0dh1);
        e3:=[];
        for x in v do
            Append(e3,HomotopyGen(x[1],x[2],x[3],x[4],x[5]));
        od;
        e:=Concatenation(h1,NegateListWord(h0dh1),e3);
    else
        e:=[];
    fi;

    h0:=VertHtpy([y]);
    v:=DelListSum(h0);
    h2:=[];
    for x in v do
        Append(h2,HomotopyGen(x[1],x[2],x[3],x[4],x[5]));
    od;

    ans:=AlgRed5(Concatenation(h0,NegateListWord(h2),e));
    slot[pr][pt][g]:=ans;
fi;

if sgn<0 then
    return NegateListWord(slot[pr][pt][g]);
fi;
return StructuralCopy(slot[pr][pt][g]);
end;
##############################################

##############################################
BaseW0:=fail;
##############################################
Homotopy:=function(k,w)
local f,g,v,q,s,r,t,h;
#HAP homotopy interface: w = [i,g] a letter of the total resolution in
#degree k; returns a word in degree k+1.
#At k=0 the raw contraction has d h(x) = x - eps(x)*[1,b] for a fixed but
#convention-dependent basepoint element b; subtracting eps(x)*w0 with
#w0 := h([1,1]) moves the basepoint to HAP's [1, 1_G].  This does not
#disturb k=1 identities because eps vanishes on boundary words.
if w=[] then return []; fi;
f:=w[1];
g:=w[2];
v:=Pair2Quad(k,f);
q:=v[1]; r:=v[2]; s:=v[3]; t:=v[4];
h:=HomotopyGen(q,s,r,t,g);
h:=List(h,x->[Quad2Pair(x[1],x[3],x[2],x[4])[2],x[5]]);
h:=AlgebraicReduction(h);
if k=0 then
    if BaseW0=fail then
        BaseW0:=HomotopyGen(0,0,1,1,1);
        BaseW0:=List(BaseW0,x->[Quad2Pair(x[1],x[3],x[2],x[4])[2],x[5]]);
        BaseW0:=AlgebraicReduction(BaseW0);
    fi;
    if SignInt(f)>0 then
        h:=AlgebraicReduction(Concatenation(h,NegateWord(1*BaseW0)));
    else
        h:=AlgebraicReduction(Concatenation(h,1*BaseW0));
    fi;
fi;
return h;
end;
##############################################

SetInfoLevel(InfoWarning,1);

return Objectify(HapResolution,
                rec(
                dimension:=Dimension,
                filteredDimension:=FilteredDimension,
                boundary:=FinalBoundary,
                homotopy:=Homotopy,
                elts:=P!.elts,
                group:=P!.group,
                pseudoBoundary:=PseudoBoundary,
                pair2Quad:=Pair2Quad,
                quad2Pair:=Quad2Pair,
                properties:=
                   [["length",N],
                    ["filtration_length",FilteredLength],
                    ["initial_inclusion",true],
                    ["reduced",true],
                    ["type","resolution"],
                    ["characteristic",prime]  ]));

end;
#####################################################################

#####################################################################
SGC_TotalHomotopyCheck:=function(R, k, letter)
# verify HAP's integral homotopy identity on the free resolution R:
#   d h + h d = id             (k > 0)
#   d h = id - [1, 1_G]        (k = 0)
local HW, diff;
HW:=function(n, w)
local out, t;
out := [];
for t in w do
    Append(out, R!.homotopy(n, t));
od;
return AlgebraicReduction(out);
end;
if k = 0 then
    diff := Concatenation(ResolutionBoundaryOfWord(R, 1, R!.homotopy(0, 1*letter)),
                          NegateWord([1*letter]), [[1, 1]]);
else
    diff := Concatenation(ResolutionBoundaryOfWord(R, k+1, R!.homotopy(k, 1*letter)),
                          HW(k-1, ResolutionBoundaryOfWord(R, k, [1*letter])),
                          NegateWord([1*letter]));
fi;
return AlgebraicReduction(diff) = [];
end;
#####################################################################

#####################################################################
# Debug fallback only: this is the old GF(2) filling homotopy and is not the
# construction-time geometric homotopy used by SGC_ResolutionSpaceGroup.
SGC_ResolutionSpaceGroupSearchHomotopy:=function(G, n)
local R;
R := ResolutionSpaceGroup(G, n);
SGC_AttachHomotopyMod2(R);
return R;
end;
#####################################################################

#####################################################################
SGC_ResolutionSpaceGroup:=function(G, n)
# Free ZG-resolution of a crystallographic matrix group with a contracting
# homotopy, built the same way as HAP's ResolutionSpaceGroup
# (CrystallographicComplex + FreeGResolution) so the boundaries agree with
# HAP's, but with (1) the fundamental-domain geometry retained on the
# non-free complex, (2) an integral cellular contraction attached to it,
# and (3) the total perturbation homotopy installed on the free resolution.
local P, R;
P := SGC_CrystallographicComplexWithGeometry(G);
SGC_AttachCellularContraction(P);
R := SGC_FreeGResolutionWithHomotopy(P, n);
return R;
end;
#####################################################################
