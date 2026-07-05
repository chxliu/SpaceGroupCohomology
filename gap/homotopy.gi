# homotopy.gi — mod-2 contracting homotopy for resolutions without R!.homotopy.
# Read after linalg.gi and BEFORE functions.gi (see read.g).
#
# ResolutionSpaceGroup (HAP, Polymake-based) returns homotopy:=fail, but
# CR_ChainMapFromCocycle (hap/lib/Rings/cocycleChainMap.gi) needs R!.homotopy
# for every cup product. SGC_AttachHomotopyMod2 attaches one computed on
# demand: h(letter) is found by recursively reducing to a boundary-filling
# problem d(w)=z and solving it over GF(2) on a growing finite patch of the
# orbit space.
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
SGC_ResolutionSpaceGroup:=function(G, n)
# Free ZG-resolution of a crystallographic matrix group via HAP's
# Polymake-based ResolutionSpaceGroup, with a mod-2 contracting homotopy
# attached. G must be an affine crystallographic matrix group acting on the
# right (e.g. AffineCrystGroupOnRight / SpaceGroupIT / SpaceGroupBBNWZ);
# the returned R!.elts are column-convention (on-left) matrices because
# CrystallographicComplex transposes its complex before returning.
local R;
R := ResolutionSpaceGroup(G, n);
SGC_AttachHomotopyMod2(R);
return R;
end;
#####################################################################
