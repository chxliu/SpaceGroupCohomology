# linalg.gi — dispatch layer for external sparse GF(2) linear algebra.
# Read after data.gi and BEFORE functions.gi (see read.g).
#
# Every SGC_* wrapper takes a plain 0/1 integer matrix (or an SGC_SparseMat
# record) and dispatches: below SGC_LINALG_THRESHOLD (nrows*ncols), or when
# the sgclinalg binary is absent, it uses the GAP-native routine; above it,
# the matrix is written in SMS sparse-triplet format, the external binary is
# invoked via Process(), and the result is read back.
#
# Contract (relaxed parity): external results are mathematically correct and
# deterministic, but need not reproduce GAP's basis choices. Ranks and
# dimensions always agree exactly with native.

SGC_LINALG_THRESHOLD := 10^6;
SGC_LINALG_CALLS := 0;

#####################################################################
SGC_GF2ToZ:=function(v)
# GF(2) vector -> 0/1 integer vector (local copy; functions.gi has GF2ToZ
# but loads after this file).
return List(v, IntFFE);
end;
#####################################################################

#####################################################################
SGC_LinalgBinary:=function()
local path;
path := Filename(DirectoriesPackageLibrary("SpaceGroupCohomology", "bin"), "sgclinalg");
if path <> fail and IsExecutableFile(path) then return path; fi;
if IsBound(GAPInfo.SystemEnvironment.SGC_LINALG_PATH) then
    path := GAPInfo.SystemEnvironment.SGC_LINALG_PATH;
    if path <> "" and IsExecutableFile(path) then return path; fi;
fi;
return fail;
end;
#####################################################################

#####################################################################
SGC_SparseMat:=function(nrows, ncols)
return rec(isSGCSparse := true, nrows := nrows, ncols := ncols, entries := []);
end;
#####################################################################

#####################################################################
SGC_IsSparseMat:=function(M)
return IsRecord(M) and IsBound(M.isSGCSparse);
end;
#####################################################################

#####################################################################
SGC_SparseToDense:=function(S)
# Scatter triplets mod 2 into a dense 0/1 integer matrix.
local M, e;
M := NullMat(S.nrows, S.ncols);
for e in S.entries do
    M[e[1]][e[2]] := (M[e[1]][e[2]] + e[3]) mod 2;
od;
return M;
end;
#####################################################################

#####################################################################
SGC_MatDims:=function(M)
if SGC_IsSparseMat(M) then return [M.nrows, M.ncols]; fi;
if Length(M) = 0 then return [0, 0]; fi;
return [Length(M), Length(M[1])];
end;
#####################################################################

#####################################################################
SGC_ShouldOffload:=function(M)
local d;
d := SGC_MatDims(M);
return d[1] * d[2] > SGC_LINALG_THRESHOLD and SGC_LinalgBinary() <> fail;
end;
#####################################################################

#####################################################################
SGC_ToNative:=function(M)
# Whatever the input form, return a plain matrix usable on the native path.
if SGC_IsSparseMat(M) then return SGC_SparseToDense(M); fi;
return M;
end;
#####################################################################

#####################################################################
SGC_WriteSMS:=function(M, path)
# Write M (plain 0/1 integer matrix, GF(2) matrix, or sparse rec) as mod-2 SMS.
local out, d, i, j, e, v;
out := OutputTextFile(path, false);
SetPrintFormattingStatus(out, false);
d := SGC_MatDims(M);
AppendTo(out, d[1], " ", d[2], " ", 2, "\n");
if SGC_IsSparseMat(M) then
    for e in M.entries do
        AppendTo(out, e[1], " ", e[2], " ", e[3], "\n");
    od;
else
    for i in [1..d[1]] do
        for j in [1..d[2]] do
            v := M[i][j];
            if not IsInt(v) then v := IntFFE(v); fi;
            if v mod 2 <> 0 then AppendTo(out, i, " ", j, " ", 1, "\n"); fi;
        od;
    od;
fi;
AppendTo(out, "0 0 0\n");
CloseStream(out);
end;
#####################################################################

#####################################################################
SGC_ReadSMS:=function(path)
# -> plain 0/1 integer matrix (entries accumulated mod 2).
local str, lines, header, M, parts, k, i, j;
str := StringFile(path);
if str = fail then Error("SGC_ReadSMS: cannot read ", path); fi;
lines := SplitString(str, "\n");
header := List(SplitString(lines[1], " "), Int);
if header[1] = 0 then return []; fi;
M := NullMat(header[1], header[2]);
for k in [2..Length(lines)] do
    if lines[k] = "" then continue; fi;
    parts := List(SplitString(lines[k], " "), Int);
    i := parts[1]; j := parts[2];
    if i = 0 and j = 0 then break; fi;
    M[i][j] := (M[i][j] + parts[3]) mod 2;
od;
return M;
end;
#####################################################################

#####################################################################
SGC_CallBinary:=function(op, infiles, outfile)
local bin, ret;
bin := SGC_LinalgBinary();
if bin = fail then Error("SGC_CallBinary: sgclinalg binary not found"); fi;
ret := Process(DirectoryCurrent(), bin, InputTextNone(), OutputTextNone(),
               Concatenation([op], infiles, [outfile]));
if ret <> 0 then Error("sgclinalg ", op, " failed with exit code ", ret); fi;
SGC_LINALG_CALLS := SGC_LINALG_CALLS + 1;
end;
#####################################################################

#####################################################################
SGC_TmpFile:=function(nam)
return Filename(DirectoryTemporary(), nam);
end;
#####################################################################

#####################################################################
SGC_RankMod2:=function(M)
local fin, fout, words;
if not SGC_ShouldOffload(M) then
    M := SGC_ToNative(M);
    if Length(M) = 0 then return 0; fi;
    return RankMat(M*Z(2));
fi;
fin := SGC_TmpFile("A.sms"); fout := SGC_TmpFile("out.txt");
SGC_WriteSMS(M, fin);
SGC_CallBinary("rank2", [fin], fout);
words := SplitString(StringFile(fout), " \n");
return Int(words[2]);
end;
#####################################################################

#####################################################################
SGC_RowBasisMod2:=function(M)
# Row-space basis over GF(2), returned as GF(2) matrix (semi-echelonized;
# external form is fully reduced RREF — span-equal, both echelon).
local fin, fout, B;
if not SGC_ShouldOffload(M) then
    M := SGC_ToNative(M);
    return BaseMat(M*Z(2));
fi;
fin := SGC_TmpFile("A.sms"); fout := SGC_TmpFile("out.sms");
SGC_WriteSMS(M, fin);
SGC_CallBinary("rowbasis2", [fin], fout);
B := SGC_ReadSMS(fout);
return B*Z(2);
end;
#####################################################################

#####################################################################
SGC_NullspaceMod2:=function(M)
# Basis of {x : x.M = 0 mod 2} as 0/1 integer vectors
# (same output convention as BasisNullspaceModN(M, 2)).
local fin, fout;
if not SGC_ShouldOffload(M) then
    M := SGC_ToNative(M);
    return List(NullspaceMat(M*Z(2)), SGC_GF2ToZ);
fi;
fin := SGC_TmpFile("A.sms"); fout := SGC_TmpFile("out.sms");
SGC_WriteSMS(M, fin);
SGC_CallBinary("nullsp2", [fin], fout);
return SGC_ReadSMS(fout);
end;
#####################################################################

#####################################################################
SGC_SteinitzMod2:=function(K, I)
# Complement basis of rowspace(I) inside rowspace(K), as GF(2) rows —
# the role of BaseSteinitzVectors(K*Z(2), I)!.factorspace in functions.gi.
# I must be echelonized (BaseMat/SGC_RowBasisMod2 output qualifies).
# CAUTION: GAP's BaseSteinitzVectors LOOPS FOREVER if K has dependent rows,
# so the native path echelonizes K first unless it is already a basis
# (byte-identical to the old behavior for genuine-basis inputs, which is
# what every call site in functions.gi passes).
local fk, fsub, fout, C, KG;
if Length(I) = 0 then
    return SGC_RowBasisMod2(K);
fi;
if not SGC_ShouldOffload(K) then
    KG := SGC_ToNative(K)*Z(2);
    if RankMat(KG) < Length(KG) then
        KG := BaseMat(KG);
    fi;
    return BaseSteinitzVectors(KG, I*Z(2))!.factorspace;
fi;
fk := SGC_TmpFile("K.sms"); fsub := SGC_TmpFile("I.sms"); fout := SGC_TmpFile("out.sms");
SGC_WriteSMS(K, fk);
SGC_WriteSMS(I, fsub);
SGC_CallBinary("steinitz2", [fk, fsub], fout);
C := SGC_ReadSMS(fout);
return C*Z(2);
end;
#####################################################################

#####################################################################
SGC_SolverMod2:=function(M)
# Factor M once; return a closure b -> SolutionMat(M*Z(2), b*Z(2)).
# The closure accepts 0/1 integer vectors or GF(2) vectors and returns a
# GF(2) solution vector (x with x.M = b) or fail.
local d, vectors, coeffs, heads, se, fin, fout, lines, k, parts, s,
      bitsToGF2, nvec;

d := SGC_MatDims(M);
if d[1] = 0 then
    return function(b)
        if ForAll(b, x -> x = 0 or x = 0*Z(2)) then return []; fi;
        return fail;
    end;
fi;

bitsToGF2 := function(s)
    return List(s, c -> (IntChar(c) - IntChar('0')) * Z(2)^0);
end;

if not SGC_ShouldOffload(M) then
    se := SemiEchelonMatTransformation(SGC_ToNative(M)*Z(2));
    vectors := List(se.vectors, ShallowCopy);
    coeffs  := List(se.coeffs, ShallowCopy);
    heads   := se.heads;
else
    fin := SGC_TmpFile("A.sms"); fout := SGC_TmpFile("out.txt");
    SGC_WriteSMS(M, fin);
    SGC_CallBinary("semiech2", [fin], fout);
    lines := SplitString(StringFile(fout), "\n");
    parts := SplitString(lines[1], " ");
    nvec := Int(parts[2]);
    vectors := []; coeffs := [];
    for k in [1..nvec] do
        Add(vectors, bitsToGF2(SplitString(lines[1+k], " ")[2]));
    od;
    for k in [1..nvec] do
        Add(coeffs, bitsToGF2(SplitString(lines[1+nvec+k], " ")[2]));
    od;
    parts := SplitString(lines[2+2*nvec], " ");
    heads := List(parts{[2..Length(parts)]}, Int);
fi;

return function(b)
    local v, x, j;
    v := ShallowCopy(b*Z(2));
    x := ListWithIdenticalEntries(d[1], 0*Z(2));
    for j in [1..Length(heads)] do
        if v[j] <> 0*Z(2) then
            if heads[j] = 0 then return fail; fi;
            AddRowVector(v, vectors[heads[j]]);
            AddRowVector(x, coeffs[heads[j]]);
        fi;
    od;
    return x;
end;
end;
#####################################################################

#####################################################################
SGC_SparseBoundaryMat:=function(R, n, transposed)
# Mod-2 boundary matrix of the degree-n resolution module as a sparse rec,
# built directly from R!.boundary (already sparse) — never densified.
# Plain:      Dim(n) x Dim(n-1)   (rows = generators of module n)
# Transposed: Dim(n-1) x Dim(n)
local S, Dimension, Boundary, i, x, j;
Dimension := R!.dimension;
Boundary := R!.boundary;
if transposed then
    S := SGC_SparseMat(Dimension(n-1), Dimension(n));
else
    S := SGC_SparseMat(Dimension(n), Dimension(n-1));
fi;
for i in [1..Dimension(n)] do
    for x in Boundary(n, i) do
        j := AbsoluteValue(x[1]);
        if transposed then
            Add(S.entries, [j, i, 1]);
        else
            Add(S.entries, [i, j, 1]);
        fi;
    od;
od;
return S;
end;
#####################################################################
