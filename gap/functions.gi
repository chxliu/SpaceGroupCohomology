# Copyright (c) Chunxiao Liu and Weicheng Ye 2024
# Reference: Chunxiao Liu, Weicheng Ye, arXiv:2410.03607v3
#
#
# In this file (gap/functions.gi) you will find:
#
# GroupCohomologyMod2(IT) or GroupCohomologyMod2(IT,R) prints the mod-2 cohomology ring presentation.
# WPCohomologyTable(IT) or WPCohomologyTable(IT,R) prints the IWP-to-degree-3-cohomology-class table.
# WPCohomologyClass(IT,WPs) or WPCohomologyClass(IT,R,WPs) prints the degree-3 class for the selected WPs.
# SpaceGroupCohomologyRingGapInterface(IT) prints the ring presentation and the full IWP cohomology table.
#
#
# Arguments of the functions:
#
# vec: any vector;
# powvec: the vector recording powers of generators in a monomial;
# powvecs: a list of vectors recording powers of generators in a monomial;
# gensdim: the vector recording number of generators at each degree;
# gensnames: a list of generator names as strings;
# sep: a string separater between monomial expressions;
# space: a string at the end of a polynomial expression;
# R: a resolution;
# deg: degree (or the highest degree) of interest;
# u,v: cocycles as class vectors;
# p,q: degree of the cocycles u,v;
# IT: the numbering of space group;
# gens: a list of cohomology ring generators (as class vectors);
#
#
# Other support functions
# GF2ToZ(vec);
# Letter2Monomial(powvec,gensdim,gensnames);
# PrintMonomialString(powvecs,gensdim,sep [,gensnames,space]);
# CR_Mod2CocyclesAndCoboundaries(R,deg);
# Mod2CupProduct(R,u,v,p,q);
# Mod2RingGenerators(R,deg);
# Mod2RingGensAndRels(IT [,3,R,gens]);
# PointGroupTranslationExtension();
# IrreducibleWyckoffPoints();
# SGC_ResolutionForIT(IT);


# NOTE: This file is part of the SpaceGroupCohomology GAP package.
# Load via: LoadPackage("SpaceGroupCohomology");
# The HAP dependency is declared in PackageInfo.g and loaded automatically.


#####################################################################
#####################################################################

GF2ToZ:=function(v)
# Convert a GF(2) vector to a 0/1 integer vector.
# IntFFE(0*Z(2))=0 and IntFFE(Z(2))=1, so this matches the old
# elementwise "=0*Z(2) -> 0 else 1" loop for any GF(2) field-element input.
return List(v, IntFFE);
end;

#####################################################################
######################################################################


#####################################################################
#####################################################################

Letter2Monomial:=function(vec,GensDim,name)
local j,k,u,v,lett;

lett := [];
u := 1;
if name = [] then
    for j in [1..Length(GensDim)] do
        if GensDim[j] = 1 then
            if vec[u]>1 then
                Append(lett,[JoinStringsWithSeparator([["A","B","C","D","E","F"][j],"^",vec[u]],"")]);
            elif vec[u]=1 then
                Append(lett,[JoinStringsWithSeparator([["A","B","C","D","E","F"][j]],"")]);
            fi;
            u := u+1;
        else
            for k in [1..GensDim[j]] do
                if vec[u]>1 then
                    Append(lett,[JoinStringsWithSeparator([["A","B","C","D","E","F"][j],String(k),"^",vec[u]],"")]);
                elif vec[u]=1 then
                    Append(lett,[JoinStringsWithSeparator([["A","B","C","D","E","F"][j],String(k)],"")]);
                fi;
                u := u+1;
            od;
        fi;
    od;
else
    for j in [1..Length(GensDim)] do
    
        for k in [1..GensDim[j]] do
            if vec[u]>1 then
                Append(lett,[JoinStringsWithSeparator([name[u],"^",vec[u]],"")]);
            elif vec[u]=1 then
                Append(lett,[JoinStringsWithSeparator([name[u]],"")]);
            fi;
            u := u+1;
        od;
    od;
fi;
return JoinStringsWithSeparator(lett, ".");
end;

#####################################################################
######################################################################


#####################################################################
#####################################################################

PrintMonomialString:=function(arg)
local poly,vecs,GensDim,sep,gennames,space;

vecs := arg[1];
GensDim := arg[2];
sep := arg[3];
gennames := [];
if Length(arg) >=4 then
    gennames := arg[4];
fi;
space:="  ";
if Length(arg) >=5 then
    space := arg[5];
fi;

poly:=List(vecs,x->Letter2Monomial(x,GensDim,gennames));
Print(JoinStringsWithSeparator(poly,sep),space);
return 0;
end;

######################################################################
#####################################################################


#####################################################################
#####################################################################

CR_Mod2CocyclesAndCoboundaries:=function(arg)
local
	R, n, toggle, Dimension, Boundary,
	M1Ts, M2Ts, CCsolver,
	kerdim, imgdim, cohdim, Mod2Cohomologydim,
    BasisKerd1, BasisImaged2, Rels, CobandCoc,
	#Smith, SmithRecord, TorsionCoefficients,
	ColMat, InvColMat,
	RemoveRowsMat, InsertRowsList,
	CycleToClass, ClassToCycle,
	CohomologyBasisMatrix,
	i, j, k, x, sum;


R:=arg[1];
n:=arg[2];
Dimension:=R!.dimension;
Boundary:=R!.boundary;
toggle := true;
if Length(arg) > 2 then       #optional 3rd arg now actually selects behavior;
    toggle := arg[3];         #all existing callers pass `true`, so behavior is unchanged.
fi;

if n <0 then return false; fi;
if n=0 then return [0]; fi;

#####################################################################


	################SPARSE BOUNDARY MATRICES######################
#M1 is Dim(n) x Dim(n-1);  M2 is Dim(n+1) x Dim(n).  Both are built directly
#from the sparse R!.boundary as transposed sparse records (gap/linalg.gi) and
#never densified when huge; below SGC_LINALG_THRESHOLD the wrappers densify
#and run the exact routines the old dense code used (BasisNullspaceModN,
#BaseMat), so small-group output is unchanged.

if Dimension(n) = 0 then
    BasisKerd1:=[];
    BasisImaged2:=[];
else
    if Dimension(n+1)>0 then
        M2Ts:=SGC_SparseBoundaryMat(R, n+1, true);   #= TransposedMat(M2), Dim(n) x Dim(n+1)
    else
        M2Ts:=SGC_SparseMat(Dimension(n), 1);        #transposed zero row, as the old code built
    fi;
    BasisKerd1:=SGC_NullspaceMod2(M2Ts);

    M1Ts:=SGC_SparseBoundaryMat(R, n, true);         #= TransposedMat(M1), Dim(n-1) x Dim(n)
    BasisImaged2:=SGC_RowBasisMod2(M1Ts);
fi;

#Print(BasisKerd1);
#Print(BasisImaged2);

imgdim:=Length(BasisImaged2);
kerdim:=Length(BasisKerd1);
cohdim:=kerdim-imgdim;


CobandCoc:=ShallowCopy(BasisImaged2);
if cohdim > 0 then
    if imgdim = 0 then
        CobandCoc := BasisKerd1 * Z(2);
    else
        Append(CobandCoc, SGC_SteinitzMod2(BasisKerd1, BasisImaged2));
    fi;
fi;

for i in [1..kerdim] do
    CobandCoc[i]:=GF2ToZ(CobandCoc[i]);
od;


if toggle=false then 
return rec(
		#cocyclesBasis:=BasisKerd1,
		#boundariesCoefficients:=Rels,
		#torsionCoefficients:=fail,
        cocyclesBasis:=CobandCoc,
        Mod2Cohomologydim:=cohdim,
		cocycleToClass:=fail,
		classToCocycle:=fail );
fi;


if cohdim > 0 then
    CohomologyBasisMatrix := CobandCoc{[kerdim-cohdim+1 .. kerdim]};
fi;

#Factor CobandCoc once and reuse for every CycleToClass call (the old code
#re-eliminated the same matrix via SolutionMat on each call). CobandCoc has
#independent rows, so the solution is unique and identical to SolutionMat's.
if cohdim > 0 then
    CCsolver := SGC_SolverMod2(CobandCoc);
fi;

#####################################################################
CycleToClass:=function(v)
local u;

if cohdim = 0 then
return [];
fi;
u:=CCsolver(v){[kerdim-cohdim+1 .. kerdim]};
return GF2ToZ(u);

end;
#####################################################################

#####################################################################
ClassToCycle:=function(u)
local w;

if cohdim = 0 then
    return List([1..Dimension(n)], x->0);
fi;

w := List(u * CohomologyBasisMatrix, x -> x mod 2);
return w;
end;
#####################################################################

return 	rec(
		#cocyclesBasis:=BasisKerd1,
	 	#boundariesCoefficients:=Rels,
        #torsionCoefficients:=TorsionCoefficients,
        cocyclesBasis:=CobandCoc,
        Mod2Cohomologydim:=cohdim,
	 	cocycleToClass:=CycleToClass,
	 	classToCocycle:=ClassToCycle );

end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################

Mod2CupProducts:=function(arg)
local
    R, u, vs, p, q, P, Q, N,
    uCocycle,
    vCocycle,
    uvCocycle,
    uChainMap,
    DimensionR,
    products, v, ww, i, w, x, sw;

    ####################BEGIN TO READ THE INPUT##################
R:=arg[1];
DimensionR:=R!.dimension;
u:=arg[2];
vs:=arg[3];
p:=arg[4];
q:=arg[5];

if Length(arg)>5 then P:=arg[6];
else
P:=CR_Mod2CocyclesAndCoboundaries(R,p,true);
fi;

if Length(arg)>6 then Q:=arg[7];
else
Q:=CR_Mod2CocyclesAndCoboundaries(R,q,true);
fi;

if Length(arg)>7 then N:=arg[8];
else
N:=CR_Mod2CocyclesAndCoboundaries(R,p+q,true);
fi;
    #####################FINISHED REAQDING THE INPUT#############

uCocycle:=P.classToCocycle(u);
uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,p,q);
ww:=[];
for i in [1..DimensionR(p+q)] do
    ww[i]:=uChainMap([[i,1]]);
od;
products:=[];
for v in vs do
    vCocycle:=Q.classToCocycle(v);
    uvCocycle:=[];
    for i in [1..DimensionR(p+q)] do
        w:=ww[i];
        sw:=0;
        for x in w do
            #sw:=sw+ SignInt(x[1])*vCocycle[AbsoluteValue(x[1])];
            sw:=sw+ vCocycle[AbsoluteValue(x[1])];
        od;
        uvCocycle[i]:=sw mod 2;
    od;
    Add(products,N.cocycleToClass(uvCocycle));
od;
return products;
end;
#####################################################################
#####################################################################

#####################################################################
#####################################################################

Mod2CupProduct:=function(arg)
local batchArgs;
batchArgs:=ShallowCopy(arg);
batchArgs[3]:=[arg[3]];
return CallFuncList(Mod2CupProducts,batchArgs)[1];
end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
# CohomologyBasis: standard basis (one vector per torsion entry). Shared by
# Mod2RingGenerators and Mod2RingGensAndRels (previously duplicated in both).
CohomologyBasis:=function(Torsion)
local i, v, Basis;
Basis:=[];
for i in [1..Length(Torsion)] do
    v:=List([1..Length(Torsion)], j->0);
    v[i]:=Torsion[i];
    Append(Basis, [v]);
od;
return Basis;
end;
#####################################################################
# AppendGenuineGens: append the genuine degree-k generators GenK to the
# degree-k basis (CupBaseK/CupBaseKLett), skipping any that are already
# expressible in the current basis. Factors the identical "append generators"
# block that was repeated for k = 2,3,4,5,6. offset = #generators of degree < k.
AppendGenuineGens:=function(SolverCache, CupBaseK, CupBaseKLett, GenK, GensLett, offset, degree)
local iu, sol;
for iu in [(offset+1)..(offset+Length(GenK))] do
    if CupBaseK = [] then
        Append(CupBaseK,[GenK[iu-offset]]);
        Append(CupBaseKLett,[GensLett[iu]]);
    else
        sol :=SGC_CachedSolveIn(SolverCache,CupBaseK,GenK[iu-offset]);
        if sol = fail then
            Append(CupBaseK,[GenK[iu-offset]]);
            Append(CupBaseKLett,[GensLett[iu]]);
        else
            Print("Error: containing fake degree-", degree, " generator(s)!!\n");
        fi;
    fi;
od;
end;
#####################################################################

Mod2RingGenerators:=function(arg)
local
        R, n, GG, IT,
        Gens, GensBasis, GensBasis1ton, Cups, Cupped, cupped, CuppedBasis, spacedim,
        uCocycle, vCocycle, uvCocycle, ww, uChainMap,
        sol, CB, CohBases,
        BasisP, BasisQ,
        i, j, p, q, u, v, ln, iu, iv, w, x, sw;

#This function computes, for a given n, the generators at degree 1,2,...,n.
#e.g. Mod2RingGenerators(IT=76,n=4);

n:=arg[2];
spacedim:=3;

if Length(arg)=3 then
    spacedim:=arg[3];
fi;

if IsInt(arg[1]) then
    IT := arg[1];
    GG := SpaceGroupIT(spacedim,IT);
    R := SGC_ResolutionSpaceGroup(GG,n+1);

else
    if IsGroup(arg[1]) then
        GG := arg[1];
        R := SGC_ResolutionSpaceGroup(GG,n+1);
    else
        R:=arg[1];
    fi;
fi;


CB:=[];
for p in [1..n] do
    CB[p]:=CR_Mod2CocyclesAndCoboundaries(R,p,true);  #CR_Mod2CocyclesAndCoboundaries gives all the cocycles followed by all the coboundaries as vectors
od;

CohBases:=List([1..n],p->CohomologyBasis(List([1..CB[p].Mod2Cohomologydim],i->1)));

GensBasis1ton :=[CohBases[1]];  # Record degree-1 generators  #Make an identity matrix, consisting of basis vectors for cocycles


for j in [2..n] do        # Then deal with degree-j generators for j=2,3,...,n

    Cups:=CohBases[j];    #Make an identity matrix, consisting of basis vectors for cocycles

    Cupped :=[];

    for p in [QuoInt(j+1,2)..(j-1)] do  #QuoInt(j,2) is gap command for the usual Int(j/2)
        q:=j-p;
        BasisP:=CohBases[p];
        BasisQ:=CohBases[q];

        iu :=1;
        for u in BasisP do

            uCocycle:=CB[p].classToCocycle(u);
            uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,p,q);
            ww:=[];
            for i in [1..(R!.dimension(j))] do
                Append(ww, [uChainMap([[i,1]])]);
            od;

            iv :=1;
            for v in BasisQ do

                if ((p > q) or (p=q and iv>=iu)) then

                    vCocycle:=CB[q].classToCocycle(v);

                    uvCocycle:=[];
                    for i in [1..(R!.dimension(j))] do
                        w:=ww[i];
                        sw:=0;
                        for x in w do
                            sw:=sw + vCocycle[AbsoluteValue(x[1])];
                        od;
                        uvCocycle[i]:=sw mod 2;
                    od;

                    cupped := CB[j].cocycleToClass(uvCocycle);
    
                    Append(Cupped,[cupped*Z(2)]);

                    #cupped :=Mod2CupProduct(R,u,v,p,q,CB[p],CB[q],CB[j]);

                fi;

                iv := iv+1;
            od;

            iu := iu+1;
        od;

    od;

    
    if Cups = [] then        #This is when cohomology dim = 0
        GensBasis1ton[j]:= [];
    else
        if Cupped = [] then
            Gens := Cups*Z(2);
        else
            CuppedBasis := List(SGC_RowBasisMod2(Cupped),ShallowCopy);
            Gens := SGC_SteinitzMod2(Cups,CuppedBasis);
        fi;

        GensBasis :=[];

        for i in [1..Length(Gens)] do
            GensBasis[i]:=GF2ToZ(Gens[i]);
        od;

        GensBasis1ton[j]:=GensBasis;
    fi;

od;

return GensBasis1ton;
end;

#####################################################################
#####################################################################

#Central special-group degree policy. Groups with degree-4 ring generators need
#relations through degree 8; groups with degree-6 ring generators need relations
#through degree 12; every other group needs relations through degree 6.
#The resolution is always built to (max relation degree)+1.
SGC_Degree4GeneratorGroups := [108, 109, 120, 130, 136, 140, 142, 197, 204, 230];
SGC_Degree6GeneratorGroups := [219, 226, 228];
SGC_MaxRelationDegreeForIT := function(IT)
if (IT in SGC_Degree4GeneratorGroups) = true then
    return 8;
elif (IT in SGC_Degree6GeneratorGroups) = true then
    return 12;
else
    return 6;
fi;
end;

#Enumerate every generator-exponent vector of total degree exactly deg, as sums
#of rows of gensLett (one unit vector per ring generator, degrees in genDeg).
#Multisets are walked in non-decreasing index order, so this is the recursive
#equivalent of the hand-nested relation-reduction loops used for r<=8, without
#their nesting-depth limit (needed for the r=9..12 search, where the full
#cartesian-product style would also be combinatorially infeasible).
SGC_MonomialsOfDegree := function(gensLett, genDeg, deg)
local result, recurse;
result := [];
recurse := function(startIdx, acc)
local idx, next;
    if acc*genDeg = deg then
        Append(result,[acc]);
        return;
    fi;
    for idx in [startIdx..Length(gensLett)] do
        next := acc + gensLett[idx];
        if next*genDeg <= deg then
            recurse(idx, next);
        fi;
    od;
end;
recurse(1, List(genDeg,x->0));
return result;
end;

#####################################################################
#####################################################################

SGC_NewRelationReducer:=function(gensLett, genDegrees, relationsByDegree, degree)
local monomials, vectorLength, monomialKeys, keyDictionary, monomialKey,
      encode, span, relationProductCount, k, relDegree, m, relation, row,
      T, i;
monomials := SGC_MonomialsOfDegree(gensLett, genDegrees, degree);
vectorLength := Length(monomials)+1;
monomialKeys := [];
monomialKey := function(term)
local key, radix, j;
    if Length(term) <> Length(genDegrees) then
        Error("SGC_NewRelationReducer: term has wrong width");
    fi;
    key := 0;
    radix := 1;
    for j in [1..Length(term)] do
        key := key + term[j]*radix;
        radix := radix*(degree+1);
    od;
    return key;
end;
if Length(monomials) = 0 then
    keyDictionary := fail;
else
    keyDictionary := NewDictionary(0,true);
    for i in [1..Length(monomials)] do
        monomialKeys[i] := monomialKey(monomials[i]);
        AddDictionary(keyDictionary,monomialKeys[i],i);
    od;
fi;
encode := function(terms)
local row, term, pos;
    row := List([1..vectorLength],x->0);
    for term in terms do
        if term*genDegrees <> degree then
            Error("SGC_NewRelationReducer: term has wrong weighted degree");
        fi;
        if keyDictionary = fail then
            Error("SGC_NewRelationReducer: no target monomial for term");
        fi;
        pos := LookupDictionary(keyDictionary,monomialKey(term));
        if pos = fail then
            Error("SGC_NewRelationReducer: no target monomial for term");
        fi;
        row[pos] := 1;
    od;
    return row;
end;
span := SGC_NewSpanTesterMod2(SGC_SparseMat(0,vectorLength));
relationProductCount := 0;
for k in [1..degree-2] do
    relDegree := degree-k;
    if relDegree <= Length(relationsByDegree)
       and IsBound(relationsByDegree[relDegree]) then
        for m in SGC_MonomialsOfDegree(gensLett,genDegrees,k) do
            for relation in relationsByDegree[relDegree] do
                row := encode(List(relation,term->m+term));
                span.addIfIndependent(row);
                relationProductCount := relationProductCount+1;
            od;
        od;
    fi;
od;
T := rec(degree := degree, monomials := monomials, vectorLength := vectorLength,
         relationProductCount := relationProductCount, rank := span.rank(),
         encode := encode);
T.contains := function(row)
    return span.contains(row);
end;
T.addIfIndependent := function(row)
local answer;
    answer := span.addIfIndependent(row);
    T.rank := span.rank();
    return answer;
end;
return T;
end;

#####################################################################
#####################################################################

Mod2RingGensAndRels:=function(arg)
local
        R,n,GG,IT,Gen1,Gen2,Gen3,Gen4,Gen5,Gen6,GensGAP6,spacedim,GenDimAll,GenDegAll,
        Gens, GensLett, zeroH, Cupped, CupRelsLett, CupTemp, CupTempLett,
        CupBase2all,CupBase2,CupBase3,CupBase4,CupBase5,CupBase6,CupBase7,CupBase8,
        CupBase9,CupBase10,CupBase11,CupBase12,
        CupBase1Lett, CupBase2Lett,CupBase3Lett,CupBase4Lett,CupBase5Lett,CupBase6Lett,CupBase7Lett,CupBase8Lett,
        CupBase9Lett,CupBase10Lett,CupBase11Lett,CupBase12Lett,
        CupRel2Lett, CupRel3Lett, CupRel4Lett, CupRel5Lett, CupRel6Lett, CupRel7Lett, CupRel8Lett,
        CupRel9Lett, CupRel10Lett, CupRel11Lett, CupRel12Lett,
        CupRelByDeg, CupStep, RelReducer, SolverCache,
        cupped,M1,sum,row,BasisImaged2,cuppedRaw,CupBase6RawCobandCoc,CupBase6Raw,
        Lett1, Lett2, mono, IO,
        uCocycle, vCocycle, uvCocycle, uChainMap, ww,
        sol, sol1, solrel, cc, CB,
        BasisP, BasisQ, SmithRecord, IToPosition,
        returnData, maxRelationDegree, relations, generatorNames, ringData,
        #NonNegativeVec,
        i,j,p,q,r,s,t,u,v,w,x,y,z, ln, ip,iq,ir,is,it,iu,iv,iw,ix,iy,iz,sw;

#Arguments: arg[1] = IT (# of space group); arg[2] = spacedim (2 or 3, default 3);
#           arg[3] = R (resolution); arg[4] = Gens (generators);
#           arg[5] = true requests structured, silent output (internal use);
#           arg[6] = optional explicit relation-degree cap (internal use);
#           arg[7] = optional generator names (internal use for an untagged R).
#Without arg[6], the highest relation degree is derived from the resolution length
#as n := Length(Size(R))-1.
#e.g.: Mod2RingGensAndRels(89);                 #builds its own resolution+generators
#e.g.: Mod2RingGensAndRels(89,3);               #3 = spacedim
#e.g.: Mod2RingGensAndRels(89,3,R89);           #(R received, generators recomputed)
#e.g.: Mod2RingGensAndRels(89,3,R89,Gens);      #standard call from SpaceGroupCohomologyRingGapInterface


if Length(arg)=0 or not ((IsInt(arg[1]) and arg[1] >= 1 and arg[1] <= 230)
                         or (arg[1] = fail and Length(arg) >= 7)) then
    Error("Mod2RingGensAndRels: IT must be 1..230\n");
fi;

returnData := Length(arg) >= 5 and arg[5] = true;
maxRelationDegree := fail;
if Length(arg) >= 6 then
    if not IsInt(arg[6]) or not arg[6] in [6,8,12] then
        Error("Mod2RingGensAndRels: relation-degree cap must be 6, 8, or 12\n");
    fi;
    maxRelationDegree := arg[6];
fi;

if Length(arg)=1 then
    spacedim:=3;
fi;

if Length(arg)>=2 then
    spacedim:=arg[2];
fi;

IT := arg[1];  #Group number in International Table for Crystallography (ITC)
#Print("Begin Group No. ", IT, ":\n");

if Length(arg)<=2 then
    GG := SpaceGroupIT(spacedim,IT);
    #Choose the resolution length (n is re-derived from the resolution below at "n:=Length(Size(R))-1").
    #This mirrors the policy in SpaceGroupCohomologyRingGapInterface via SGC_MaxRelationDegreeForIT.
    n := SGC_MaxRelationDegreeForIT(IT);
    R := SGC_ResolutionSpaceGroup(GG,n+1);                #Construct resolution;
    Gens:=Mod2RingGenerators(R,4,spacedim);               #Calculate generators
elif Length(arg) = 3 then                                 #resolution given, generators not
    R := arg[3];                                          #receive resolution from input;
    Gens:=Mod2RingGenerators(R,4,spacedim);               #recompute generators from the given resolution;
else                                                      #Length(arg) >= 4
    R := arg[3];                                          #receive resolution from input;
    Gens:=arg[4];                                         #receive generators from input;
fi;

Gen1:=Gens[1];
Gen2:=Gens[2];
Gen3:=Gens[3];
Gen4:=Gens[4];
Gen5:=[]; #based on the posteriori fact that no space group has degree-5 generators.
Gen6:=[];

#Degree-6 ring generators are selected exactly the way Gen4 is: Mod2RingGenerators
#completes a basis of H^6 from classes outside the decomposable (cup-product) span.
#CR_Mod2CocyclesAndCoboundaries(R,6) inside needs the degree-7 module, hence the
#resolution-length guard (Category C resolutions are built to length 13 anyway).
if Length(Gens) >= 6 and IsBound(Gens[6]) then
    Gen6 := Gens[6];
elif (IT in SGC_Degree6GeneratorGroups) = true and Length(Size(R)) >= 7 then
    GensGAP6 := Mod2RingGenerators(R,6,spacedim);
    Gen6 := GensGAP6[6];
fi;

GenDimAll:=[Length(Gen1),Length(Gen2),Length(Gen3),Length(Gen4),Length(Gen5),Length(Gen6)];
GenDegAll:=Concatenation(List([1..Length(Gen1)],x->1),List([1..Length(Gen2)],x->2),List([1..Length(Gen3)],x->3),List([1..Length(Gen4)],x->4),List([1..Length(Gen5)],x->5),List([1..Length(Gen6)],x->6));


#####################################################################
IToPosition:=function(v)
local v0,k;

v0:=[];
for k in [1..Length(v)] do
    if v[k] = 1 then
    Append(v0, [ k ]);
    fi;
od;
return v0;
end;
######################################################################
#local ps,k;
#
#end;
######################################################################


GensLett:=CohomologyBasis(List([1..(Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5)+Length(Gen6))],i->1));
#GensLett records the powers, later to be used to convert generators/relations to letters.

SolverCache:=SGC_NewSolverCache();
CupRelByDeg:=[];


n:=Length(Size(R))-1;     #This is the highest degree at which the relations are calculated for groups No. <= 220.
                          #for 221, Length(Size(R)) = 6, n = 5, and we will use the "raw code" to calculate relations at degree 6.
if maxRelationDegree <> fail then
    if maxRelationDegree > n then
        Error("Mod2RingGensAndRels: resolution only supports relations through degree ",
              n, ", but degree ", maxRelationDegree, " was requested\n");
    fi;
    n := maxRelationDegree;
fi;
if IT = fail and n < 6 then
    Error("Mod2RingGensAndRels: an untagged resolution must support relations through degree 6\n");
fi;


CB:=[];
for p in [1..n] do
CB[p]:=CR_Mod2CocyclesAndCoboundaries(R,p,true);
od;

#Cache the all-zero class vector at each degree once. The coboundary test below
#("is u-cup-v zero in cohomology?") otherwise rebuilds the zero vector.
#on every cup product. Degrees 1..n are exactly those used (and computable) on every path.
zeroH:=[];
for p in [1..n] do
    zeroH[p]:=List([1..CB[p].Mod2Cohomologydim],x->0);
od;


####################### r = 1 ##########################


CupBase1Lett:= List([1..Length(Gen1)],x->GensLett[x]);


####################### r = 2 ##########################

CupRelsLett := [];

CupBase2 := [];
CupBase2Lett := [];
CupRel2Lett := [];
Cupped := List(Gen1,u->Mod2CupProducts(R,u,Gen1,1,1,CB[1],CB[1],CB[2]));


iu :=1;                                           #First step: calculate u-cup-v for u != v
for u in Gen1 do
    iv :=1;
    for v in Gen1 do
        Lett1 := GensLett[iu]+GensLett[iv];
        if iv>iu then
            cupped := Cupped[iu][iv];
    
            if cupped = zeroH[2] then   #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel2Lett,[[Lett1]]);
            else
                if CupBase2 = [] then             #if no basis yet then push in the genuine cocycle u-cup-v
                    Append(CupBase2,[cupped]);
                    Append(CupBase2Lett,[Lett1]);
                else
                    sol :=SGC_CachedSolveIn(SolverCache,CupBase2,cupped);
                    if sol = fail then                        #if u-cup-v is a genuine new cocycle
                        Append(CupBase2,[cupped]);
                        Append(CupBase2Lett,[Lett1]);
                    else                       #if u-cup-v is expressable by other cocycles in basis
                        solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase2Lett[x]));
                        if (Lett1 in CupBase2Lett) = false then
                            Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                            Append(CupRel2Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
                        fi;
                    fi;
                fi;
            fi;
        fi;
        iv := iv+1;
    od;
    iu := iu+1;
od;


iu :=1;                                    #Second step: calculate u-cup-u
for u in Gen1 do
    Lett1 := GensLett[iu]+GensLett[iu];
    cupped := Cupped[iu][iu];
    
        if cupped = zeroH[2] then   #if u-cup-u is a coboundary
            Append(CupRelsLett,[Lett1]);
            Append(CupRel2Lett,[[Lett1]]);
        else
            if CupBase2 = [] then             #if no basis yet then push in the genuine cocycle u-cup-u
                Append(CupBase2,[cupped]);
                Append(CupBase2Lett,[Lett1]);
            else
                sol :=SGC_CachedSolveIn(SolverCache,CupBase2,cupped);
                if sol = fail then                                  #if u-cup-u is a genuine new cocycle
                    Append(CupBase2,[cupped]);
                    Append(CupBase2Lett,[Lett1]);
                else                              #if u-cup-u is expressable by other cocycles in basis
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase2Lett[x]));
                    if (Lett1 in CupBase2Lett) = false then
                        Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                        Append(CupRel2Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
                    fi;
                fi;
            fi;
        fi;
    iu := iu+1;
od;

AppendGenuineGens(SolverCache, CupBase2, CupBase2Lett, Gen2, GensLett, Length(Gen1), 2);
CupRelByDeg[2] := CupRel2Lett;


if Length(CupBase2Lett) = CB[2].Mod2Cohomologydim then
    Print("");#Print("dim(H^2)=", CB[2].Mod2Cohomologydim,", ");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^2) = ", Length(CupBase2Lett) - CB[2].Mod2Cohomologydim,"\n");
fi;

##Begin printing the relations at deg 2:
if Length(CupRel2Lett) > 0 then
fi;
##End printing the relations at deg 2.


####################### r = 3 ##########################


CupBase3 :=[];
CupBase3Lett :=[];
CupRel3Lett := [];
CupTemp := [];
CupTempLett := [];
RelReducer := SGC_NewRelationReducer(GensLett,GenDegAll,CupRelByDeg,3);
if IsBoundGlobal("SGC_RELATION_REDUCER_OBSERVER") then
    ValueGlobal("SGC_RELATION_REDUCER_OBSERVER")(rec(degree := RelReducer.degree, relationProductCount := RelReducer.relationProductCount, rank := RelReducer.rank));
fi;


#Begins: degree-2 cup with degree-1-gen
#
#
iu :=1;
for u in CupBase2 do
    Cupped := Mod2CupProducts(R,u,Gen1,2,1,CB[2],CB[1],CB[3]);
    iv :=1;
    for v in Gen1 do
        cupped :=Cupped[iv];      #then calculate u-cup-v
        Lett1 := CupBase2Lett[iu] + GensLett[iv];
                
        ####
        ####


        #### Start: determine whether u-cup-v goes to the basis
        ####
        solrel := [];
        if (cupped = zeroH[3]) then         #if u-cup-v is a coboundary
            solrel := [Lett1];
        else
            if CupBase3 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
                Append(CupBase3,[cupped]);
                Append(CupBase3Lett,[Lett1]);
            else
                sol :=SGC_CachedSolveIn(SolverCache,CupBase3,cupped);
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase3,[cupped]);
                    Append(CupBase3Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase3Lett[x]));
                fi;
            fi;
        fi;
        ####
        #### End:  determine whether u-cup-v goes to the basis
                    
                    
        #### Start: determine whether u-cup-v gives a new relation
        ####
        if (Lett1 in CupBase3Lett) = false then
                    
            if RelReducer.addIfIndependent(RelReducer.encode(solrel)) then
                Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                Append(CupRel3Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
            fi;
        fi;
        ####
        #### End: determine whether u-cup-v gives a new relation
        
        
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Finished: degree-2 cup with degree-1-gen


AppendGenuineGens(SolverCache, CupBase3, CupBase3Lett, Gen3, GensLett, Length(Gen1)+Length(Gen2), 3);
CupRelByDeg[3] := CupRel3Lett;


#


if Length(CupBase3Lett) = CB[3].Mod2Cohomologydim then
    Print("");#Print("dim(H^3)=", CB[3].Mod2Cohomologydim,", ");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^3) = ", Length(CupBase3Lett) - CB[3].Mod2Cohomologydim,"\n");
fi;


##Begin printing the relations at deg 3:
if Length(CupRel3Lett) > 0 then
fi;
##End printing the relations at deg 3.


####################### r = 4 ##########################

CupBase4 :=[];
CupBase4Lett :=[];
CupRel4Lett := [];
CupTemp := [];
CupTempLett := [];
RelReducer := SGC_NewRelationReducer(GensLett,GenDegAll,CupRelByDeg,4);
if IsBoundGlobal("SGC_RELATION_REDUCER_OBSERVER") then
    ValueGlobal("SGC_RELATION_REDUCER_OBSERVER")(rec(degree := RelReducer.degree, relationProductCount := RelReducer.relationProductCount, rank := RelReducer.rank));
fi;


#Step-1 begins here: degree-3 cup with degree-1-gen
#
#
iu :=1;
for u in CupBase3 do
    Cupped := Mod2CupProducts(R,u,Gen1,3,1,CB[3],CB[1],CB[4]);
    iv :=1;
    for v in Gen1 do
        cupped :=Cupped[iv];      #then calculate u-cup-v
        Lett1 := CupBase3Lett[iu] + GensLett[iv];
        
        ####
        ####
    

        #### Start: determine whether u-cup-v goes to the basis
        ####
        solrel := [];
        if cupped = zeroH[4] then         #if u-cup-v is a coboundary
            solrel := [Lett1];
            
        elif CupBase4 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
            Append(CupBase4,[cupped]);
            Append(CupBase4Lett,[Lett1]);
        else
            sol :=SGC_CachedSolveIn(SolverCache,CupBase4,cupped);
            if sol = fail then                  #if u-cup-v is a genuine new cocycle
                Append(CupBase4,[cupped]);
                Append(CupBase4Lett,[Lett1]);
            else                                #if u-cup-v is expressable by other cocycles in basis
                solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase4Lett[x]));
            fi;
        fi;
        ####
        #### End:  determine whether u-cup-v goes to the basis
            
            
        #### Start: determine whether u-cup-v gives a new relation
        ####
        if (Lett1 in CupBase4Lett) = false then
                    
            if RelReducer.addIfIndependent(RelReducer.encode(solrel)) then
                Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                Append(CupRel4Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
            fi;
        fi;
        ####
        #### End: determine whether u-cup-v gives a new relation
        
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-1 finished: degree-3 cup with degree-1-gen


#


#Step-2 begins here: degree-2-gen cup with degree-2-gen
#
#
iu :=1;
for u in Gen2 do
    iv :=1;
    for v in Gen2 do
        if iv>=iu then
        
            cupped :=Mod2CupProduct(R,u,v,2,2,CB[2],CB[2],CB[4]);      #then calculate u-cup-v
        
            Lett1 := GensLett[Length(Gen1)+iu] + GensLett[Length(Gen1)+iv];
        
            if cupped = zeroH[4] then         #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel4Lett,[[Lett1]]);
            elif CupBase4 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
                Append(CupBase4,[cupped]);
                Append(CupBase4Lett,[Lett1]);
            else
                sol :=SGC_CachedSolveIn(SolverCache,CupBase4,cupped);
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase4,[cupped]);
                    Append(CupBase4Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase4Lett[x]));
                    if (Lett1 in CupBase4Lett) = false then
                        Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                        Append(CupRel4Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
                    fi;
                fi;
            fi;
        fi;
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-2 finished: degree-2-gen cup with degree-2-gen

AppendGenuineGens(SolverCache, CupBase4, CupBase4Lett, Gen4, GensLett, Length(Gen1)+Length(Gen2)+Length(Gen3), 4);
CupRelByDeg[4] := CupRel4Lett;


if Length(CupBase4Lett) = CB[4].Mod2Cohomologydim then
    Print("");#Print("dim(H^4)=", CB[4].Mod2Cohomologydim,", ");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^4) = ", Length(CupBase4Lett) - CB[4].Mod2Cohomologydim,"\n");
fi;


##Begin printing the relations at deg 4:
if Length(CupRel4Lett) > 0 then
fi;
##End printing the relations at deg 4.


####################### r = 5 ##########################

CupBase5 :=[];
CupBase5Lett :=[];
CupRel5Lett := [];
CupTemp := [];
CupTempLett := [];
RelReducer := SGC_NewRelationReducer(GensLett,GenDegAll,CupRelByDeg,5);
if IsBoundGlobal("SGC_RELATION_REDUCER_OBSERVER") then
    ValueGlobal("SGC_RELATION_REDUCER_OBSERVER")(rec(degree := RelReducer.degree, relationProductCount := RelReducer.relationProductCount, rank := RelReducer.rank));
fi;


#Step-1 begins here: degree-4 cup with degree-1-gen
#
#
iu :=1;
for u in CupBase4 do

    #### implementing Mod2CupProduct(R,u,v,4,1,CB[4],CB[1],CB[5]) -- part 1: ####
    ####
    uCocycle:=CB[4].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,4,1);
    ww:=[];
    for i in [1..(R!.dimension(5))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,4,1,CB[4],CB[1],CB[5])  -- part 1 ended ####

    iv :=1;
    for v in Gen1 do
    
        
        #### implementing Mod2CupProduct(R,u,v,4,1,CB[4],CB[1],CB[5]) -- part 2: ####
        ####
        vCocycle:=CB[1].classToCocycle(v);
        uvCocycle:=[];
        for i in [1..(R!.dimension(5))] do
            w:=ww[i];
            sw:=0;
            for x in w do
                sw:=sw + vCocycle[AbsoluteValue(x[1])];
            od;
            uvCocycle[i]:=sw mod 2;
        od;
        cupped := CB[5].cocycleToClass(uvCocycle);
        ####
        #### implementing Mod2CupProduct(R,u,v,4,1,CB[4],CB[1],CB[5]) -- part 2 ended ####

        Lett1 := CupBase4Lett[iu] + GensLett[iv];


        ####
        
        ####

        
        #### Start: determine whether u-cup-v goes to the basis
        ####
        solrel := [];
        if cupped = zeroH[5] then         #if u-cup-v is a coboundary
            solrel := [Lett1];
            
        elif CupBase5 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
            Append(CupBase5,[cupped]);
            Append(CupBase5Lett,[Lett1]);
        else
            sol :=SGC_CachedSolveIn(SolverCache,CupBase5,cupped);
            if sol = fail then                  #if u-cup-v is a genuine new cocycle
                Append(CupBase5,[cupped]);
                Append(CupBase5Lett,[Lett1]);
            else                                #if u-cup-v is expressable by other cocycles in basis
                solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase5Lett[x]));
            fi;
        fi;
        ####
        #### End:  determine whether u-cup-v goes to the basis
            
            
        #### Start: determine whether u-cup-v gives a new relation
        ####
        if (Lett1 in CupBase5Lett) = false then
            if RelReducer.addIfIndependent(RelReducer.encode(solrel)) then
                Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                Append(CupRel5Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
            fi;
        fi;
        ####
        #### End: determine whether u-cup-v gives a new relation
        
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-1 finished: degree-4 cup with degree-1-gen


#Step-2 begins here: degree-3-gen cup with degree-2-gen
#
#
iu :=1;
for u in Gen3 do

    #### implementing Mod2CupProduct(R,u,v,3,2,CB[3],CB[2],CB[5]) -- part 1: ####
    ####
    uCocycle:=CB[3].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,3,2);
    ww:=[];
    for i in [1..(R!.dimension(5))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,3,2,CB[3],CB[2],CB[5])  -- part 1 ended ####

    iv :=1;
    for v in Gen2 do
        
            
        #### implementing Mod2CupProduct(R,u,v,3,2,CB[3],CB[2],CB[5]) -- part 2: ####
        ####
        vCocycle:=CB[2].classToCocycle(v);
        uvCocycle:=[];
        for i in [1..(R!.dimension(5))] do
            w:=ww[i];
            sw:=0;
            for x in w do
                sw:=sw + vCocycle[AbsoluteValue(x[1])];
            od;
            uvCocycle[i]:=sw mod 2;
        od;
        cupped := CB[5].cocycleToClass(uvCocycle);
        ####
        #### implementing Mod2CupProduct(R,u,v,3,2,CB[3],CB[2],CB[5]) -- part 2 ended ####
        
        Lett1 := GensLett[Length(Gen1)+Length(Gen2)+iu] + GensLett[Length(Gen1)+iv];
                
        
        #### Start: determine whether u-cup-v goes to the basis
        ####
        solrel := [];
        if cupped = zeroH[5] then         #if u-cup-v is a coboundary
            solrel := [Lett1];
            
        elif CupBase5 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
            Append(CupBase5,[cupped]);
            Append(CupBase5Lett,[Lett1]);
        else
            sol :=SGC_CachedSolveIn(SolverCache,CupBase5,cupped);
            if sol = fail then                  #if u-cup-v is a genuine new cocycle
                Append(CupBase5,[cupped]);
                Append(CupBase5Lett,[Lett1]);
            else                                #if u-cup-v is expressable by other cocycles in basis
                solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase5Lett[x]));
            fi;
        fi;
        ####
        #### End:  determine whether u-cup-v goes to the basis
            
            
        #### Start: determine whether u-cup-v gives a new relation
        ####
        if (Lett1 in CupBase5Lett) = false then
            if RelReducer.addIfIndependent(RelReducer.encode(solrel)) then
                Append(CupRelsLett,[Lett1]);  #Record the letter that appears on the LHS of the relation equation
                Append(CupRel5Lett,[solrel]); #Record the letters that appear in the relation equation. Note that all but the first letters are in the cohomology basis (that we choose).
            fi;
        fi;
        ####
        #### End: determine whether u-cup-v gives a new relation
        
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-2 finished: degree-3-gen cup with degree-2-gen


AppendGenuineGens(SolverCache, CupBase5, CupBase5Lett, Gen5, GensLett, Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4), 5);
CupRelByDeg[5] := CupRel5Lett;


#


if Length(CupBase5Lett) = CB[5].Mod2Cohomologydim then
    Print("");#Print("dim(H^5)=", CB[5].Mod2Cohomologydim,", ");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^5) = ", Length(CupBase5Lett) - CB[5].Mod2Cohomologydim,"\n");
fi;


##Begin printing the relations at deg 5:
if Length(CupRel5Lett) >0 then
fi;
##End printing the relations at deg 5.


####################### r = 6 ##########################

CupBase6 :=[];
CupBase6Lett :=[];
CupRel6Lett := [];
CupTemp := [];
CupTempLett := [];
M1:=[];
CupBase6Raw:=[];

#Below starts the "Raw Method" for the degree 6 relations for groups 221-230, part1
#
if n = 5 then

    #M1 is Dim(R6) x Dim(R5);

    for i in [1..R!.dimension(6)] do
        row:=List([1..R!.dimension(5)],j->0);
        for x in R!.boundary(6,i) do          #scatter the sparse boundary directly into its target column
            j:=AbsoluteValue(x[1]);
            row[j]:= row[j] + SignInt(x[1]);
        od;
        for j in [1..R!.dimension(5)] do
            row[j]:= RemInt(row[j],2);
        od;
    M1[i]:=row;
    od;
    
    if M1 = [ ] then
        BasisImaged2:=[];
    else
        BasisImaged2:=BaseMat(TransposedMat(M1)*Z(2));
    fi;

    CupBase6RawCobandCoc:=[];
    for i in [1..Length(BasisImaged2)] do
        Append(CupBase6RawCobandCoc,[BasisImaged2[i]]);
    od;

fi;

#
#Above ends the "Raw Method" for the degree 6 relations for groups 221-230, part1


RelReducer := SGC_NewRelationReducer(GensLett,GenDegAll,CupRelByDeg,6);
if IsBoundGlobal("SGC_RELATION_REDUCER_OBSERVER") then
    ValueGlobal("SGC_RELATION_REDUCER_OBSERVER")(rec(degree := RelReducer.degree, relationProductCount := RelReducer.relationProductCount, rank := RelReducer.rank));
fi;


# Below is the good method when resolution at degree 7 has been constructed:
#
#
#
if n >= 6 then

#Step-1 begins here: degree-5 cup with degree-1-gen
#
#
iu :=1;
for u in CupBase5 do

    #### implementing Mod2CupProduct(R,u,v,5,1,CB[5],CB[1],CB[6]) -- part 1: ####
    ####
    uCocycle:=CB[5].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,5,1);
    ww:=[];
    for i in [1..(R!.dimension(6))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,5,1,CB[5],CB[1],CB[6])  -- part 1 ended ####

    iv :=1;
    for v in Gen1 do
    
    
        #### implementing Mod2CupProduct(R,u,v,5,1,CB[5],CB[1],CB[6]) -- part 2: ####
        ####
        vCocycle:=CB[1].classToCocycle(v);
        uvCocycle:=[];
        for i in [1..(R!.dimension(6))] do
            w:=ww[i];
            sw:=0;
            for x in w do
                sw:=sw + vCocycle[AbsoluteValue(x[1])];
            od;
            uvCocycle[i]:=sw mod 2;
        od;
        cupped := CB[6].cocycleToClass(uvCocycle);
        ####
        #### implementing Mod2CupProduct(R,u,v,5,1,CB[5],CB[1],CB[6]) -- part 2 ended ####

        Lett1 := CupBase5Lett[iu] + GensLett[iv];

        
        #### Start: determine whether u-cup-v goes to the basis
        ####
        solrel := [];
        if cupped = zeroH[6] then         #if u-cup-v is a coboundary
            solrel := [Lett1];
            
        elif CupBase6 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
            Append(CupBase6,[cupped]);
            Append(CupBase6Lett,[Lett1]);
        else
            sol :=SGC_CachedSolveIn(SolverCache,CupBase6,cupped);
            if sol = fail then                  #if u-cup-v is a genuine new cocycle
                Append(CupBase6,[cupped]);
                Append(CupBase6Lett,[Lett1]);
            else                                #if u-cup-v is expressable by other cocycles in basis
                solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase6Lett[x]));
            fi;
        fi;
        ####
        #### End:  determine whether u-cup-v goes to the basis
            
            
        #### Start: determine whether u-cup-v gives a new relation
        ####
        if (Lett1 in CupBase6Lett) = false then
                    
            if RelReducer.addIfIndependent(RelReducer.encode(solrel)) then
                Append(CupRelsLett,[Lett1]);
                Append(CupRel6Lett,[solrel]);
            fi;
        fi;
        ####
        #### End: determine whether u-cup-v gives a new relation

    
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-1 finished: degree-5 cup with degree-1-gen


#Step-2 begins here: (degree-2 cup with degree-2) cup with degree-2-gen, or degree-4 cup with degree-2-gen
#
#
iu :=1;
for u in CupBase4 do
    if List([1..Length(Gen1)],x->CupBase4Lett[iu][x]) = List([1..Length(Gen1)],x->0) then   #if u is of the form Bi-cup-Bj or D
        #### implementing Mod2CupProduct(R,u,v,4,2,CB[4],CB[2],CB[6]) -- part 1: ####
        ####
        uCocycle:=CB[4].classToCocycle(u);
        uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,4,2);
        ww:=[];
        for i in [1..(R!.dimension(6))] do
            Append(ww, [uChainMap([[i,1]])]);
        od;
        ####
        #### implementing Mod2CupProduct(R,u,v,4,2,CB[4],CB[2],CB[6])  -- part 1 ended ####

        iv :=1;
        for v in Gen2 do
    
        
            #### implementing Mod2CupProduct(R,u,v,4,2,CB[4],CB[2],CB[6]) -- part 2: ####
            ####
            vCocycle:=CB[2].classToCocycle(v);
            uvCocycle:=[];
            for i in [1..(R!.dimension(6))] do
                w:=ww[i];
                sw:=0;
                for x in w do
                    sw:=sw + vCocycle[AbsoluteValue(x[1])];
                od;
                uvCocycle[i]:=sw mod 2;
            od;
            cupped := CB[6].cocycleToClass(uvCocycle);
            ####
            #### implementing Mod2CupProduct(R,u,v,4,2,CB[4],CB[2],CB[6]) -- part 2 ended ####
        

            Lett1 := CupBase4Lett[iu] + GensLett[Length(Gen1)+iv];

        
            #### Start: determine whether u-cup-v goes to the basis
            ####
            solrel := [];
            if cupped = zeroH[6] then         #if u-cup-v is a coboundary
                solrel := [Lett1];
              
            elif CupBase6 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
                Append(CupBase6,[cupped]);
                Append(CupBase6Lett,[Lett1]);
            else
                sol :=SGC_CachedSolveIn(SolverCache,CupBase6,cupped);
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase6,[cupped]);
                    Append(CupBase6Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase6Lett[x]));
                fi;
            fi;
            ####
            #### End:  determine whether u-cup-v goes to the basis
            
            
            #### Start: determine whether u-cup-v gives a new relation
            ####
            if (Lett1 in CupBase6Lett) = false then
               
                if RelReducer.addIfIndependent(RelReducer.encode(solrel)) then
                    Append(CupRelsLett,[Lett1]);
                    Append(CupRel6Lett,[solrel]);
                fi;
            fi;
            ####
            #### End: determine whether u-cup-v gives a new relation
        
            iv := iv+1;
        od;
    fi;
    iu := iu+1;
od;
#
#
#Step-2 finished: (degree-2 cup with degree-2) cup with degree-2-gen, or degree-4 cup with degree-2-gen


#Step-3 begins here: degree-3-gen cup with degree-3-gen
#
#
iu :=1;
for u in Gen3 do

    #### implementing Mod2CupProduct(R,u,v,3,3,CB[3],CB[3],CB[6]) -- part 1: ####
    ####
    uCocycle:=CB[3].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,3,3);
    ww:=[];
    for i in [1..(R!.dimension(6))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,3,3,CB[3],CB[3],CB[6])  -- part 1 ended ####

    iv :=1;
    for v in Gen3 do
        if iv>=iu then
            Lett1 := GensLett[Length(Gen1)+Length(Gen2)+iu] + GensLett[Length(Gen1)+Length(Gen2)+iv];
            #Of course now u-cup-v does not contain patterns in CupRelsLett:
        
            
            #### implementing Mod2CupProduct(R,u,v,3,3,CB[3],CB[3],CB[6]) -- part 2: ####
            ####
            vCocycle:=CB[3].classToCocycle(v);
            uvCocycle:=[];
            for i in [1..(R!.dimension(6))] do
                w:=ww[i];
                sw:=0;
                for x in w do
                    sw:=sw + vCocycle[AbsoluteValue(x[1])];
                od;
                uvCocycle[i]:=sw mod 2;
            od;
            cupped := CB[6].cocycleToClass(uvCocycle);
            ####
            #### implementing Mod2CupProduct(R,u,v,3,3,CB[3],CB[3],CB[6]) -- part 2 ended ####
        
            if cupped = zeroH[6] then         #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel6Lett,[[Lett1]]);
            elif CupBase6 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
                Append(CupBase6,[cupped]);
                Append(CupBase6Lett,[Lett1]);
            else
                sol :=SGC_CachedSolveIn(SolverCache,CupBase6,cupped);
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase6,[cupped]);
                    Append(CupBase6Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase6Lett[x]));
                    if (Lett1 in CupBase6Lett) = false then
                        Append(CupRelsLett,[Lett1]);
                        Append(CupRel6Lett,[solrel]);
                    fi;
                fi;
            fi;
        fi;
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-3 finished: degree-3-gen cup with degree-3-gen


AppendGenuineGens(SolverCache, CupBase6, CupBase6Lett, Gen6, GensLett, Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5), 6);


if Length(CupBase6) = CB[6].Mod2Cohomologydim then
    Print("");#Print("dim(H^6)=", CB[6].Mod2Cohomologydim,".\n");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^6) = ", Length(CupBase6) - CB[6].Mod2Cohomologydim,"\n");
fi;

###########################
###########################
###########################
else    ##Below starts the "Raw Method" for the degree 6 relations for groups 221-230, part2
###########################
###########################
###########################

#Step-1 begins here: degree-5 cup with degree-1-gen
#
#
iu :=1;
for u in CupBase5 do

    #### implementing Mod2CupProduct(R,u,v,5,1,CB[5],CB[1],CB[6]) -- part 1: ####
    ####
    uCocycle:=CB[5].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,5,1);
    ww:=[];
    for i in [1..(R!.dimension(6))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,5,1,CB[5],CB[1],CB[6])  -- part 1 ended ####

    iv :=1;
    for v in Gen1 do
    
    
        #### implementing Mod2CupProduct(R,u,v,5,1,CB[5],CB[1],CB[6]) -- part 2: ####
        ####
        vCocycle:=CB[1].classToCocycle(v);
        uvCocycle:=[];
        for i in [1..(R!.dimension(6))] do
            w:=ww[i];
            sw:=0;
            for x in w do
                sw:=sw + vCocycle[AbsoluteValue(x[1])];
            od;
            uvCocycle[i]:=sw mod 2;
        od;
        cuppedRaw := uvCocycle;
        ####
        #### implementing Mod2CupProduct(R,u,v,5,1,CB[5],CB[1],CB[6]) -- part 2 ended ####

        Lett1 := CupBase5Lett[iu] + GensLett[iv];

        
        #### Start: determine whether u-cup-v goes to the basis
        ####
        solrel := [];
        if (SGC_CachedSolveIn(SolverCache,BasisImaged2,cuppedRaw) = fail) = false then        #if u-cup-v is a coboundary
            solrel := [Lett1];
            
        elif CupBase6Raw = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
            Append(CupBase6Raw,[cuppedRaw]);
            Append(CupBase6Lett,[Lett1]);
            Append(CupBase6RawCobandCoc,[cuppedRaw]);
        else
            sol :=SGC_CachedSolveIn(SolverCache,CupBase6RawCobandCoc,cuppedRaw);
            if sol = fail then                  #if u-cup-v is a genuine new cocycle
                Append(CupBase6Raw,[cuppedRaw]);
                Append(CupBase6Lett,[Lett1]);
                Append(CupBase6RawCobandCoc,[cuppedRaw]);
            else                                #if u-cup-v is expressable by other cocycles in basis
                sol1:=List([(Length(BasisImaged2)+1)..Length(sol)],x->sol[x]);
                solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol1)),x->CupBase6Lett[x]));
            fi;
        fi;
        ####
        #### End:  determine whether u-cup-v goes to the basis
            
            
        #### Start: determine whether u-cup-v gives a new relation
        ####
        if (Lett1 in CupBase6Lett) = false then
                    
            if RelReducer.addIfIndependent(RelReducer.encode(solrel)) then
                Append(CupRelsLett,[Lett1]);
                Append(CupRel6Lett,[solrel]);
            fi;
        fi;
        ####
        #### End: determine whether u-cup-v gives a new relation

    
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-1 finished: degree-5 cup with degree-1-gen


#Step-2 begins here: (degree-2 cup with degree-2) cup with degree-2-gen, or degree-4 cup with degree-2-gen
#
#
iu :=1;
for u in CupBase4 do
    if List([1..Length(Gen1)],x->CupBase4Lett[iu][x]) = List([1..Length(Gen1)],x->0) then   #if u is of the form Bi-cup-Bj or D
        #### implementing Mod2CupProduct(R,u,v,4,2,CB[4],CB[2],CB[6]) -- part 1: ####
        ####
        uCocycle:=CB[4].classToCocycle(u);
        uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,4,2);
        ww:=[];
        for i in [1..(R!.dimension(6))] do
            Append(ww, [uChainMap([[i,1]])]);
        od;
        ####
        #### implementing Mod2CupProduct(R,u,v,4,2,CB[4],CB[2],CB[6])  -- part 1 ended ####

        iv :=1;
        for v in Gen2 do
    
        
            #### implementing Mod2CupProduct(R,u,v,4,2,CB[4],CB[2],CB[6]) -- part 2: ####
            ####
            vCocycle:=CB[2].classToCocycle(v);
            uvCocycle:=[];
            for i in [1..(R!.dimension(6))] do
                w:=ww[i];
                sw:=0;
                for x in w do
                    sw:=sw + vCocycle[AbsoluteValue(x[1])];
                od;
                uvCocycle[i]:=sw mod 2;
            od;
            cuppedRaw := uvCocycle;
            ####
            #### implementing Mod2CupProduct(R,u,v,4,2,CB[4],CB[2],CB[6]) -- part 2 ended ####
        

            Lett1 := CupBase4Lett[iu] + GensLett[Length(Gen1)+iv];

        
            #### Start: determine whether u-cup-v goes to the basis
            ####
            solrel := [];
            if (SGC_CachedSolveIn(SolverCache,BasisImaged2,cuppedRaw) = fail) = false then        #if u-cup-v is a coboundary
                solrel := [Lett1];
            
            elif CupBase6Raw = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
                Append(CupBase6Raw,[cuppedRaw]);
                Append(CupBase6Lett,[Lett1]);
                Append(CupBase6RawCobandCoc,[cuppedRaw]);
            else
                sol :=SGC_CachedSolveIn(SolverCache,CupBase6RawCobandCoc,cuppedRaw);
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase6Raw,[cuppedRaw]);
                    Append(CupBase6Lett,[Lett1]);
                    Append(CupBase6RawCobandCoc,[cuppedRaw]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    sol1:=List([(Length(BasisImaged2)+1)..(Length(sol))],x->sol[x]);
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol1)),x->CupBase6Lett[x]));
                fi;
            fi;
            ####
            #### End:  determine whether u-cup-v goes to the basis
            
            
            #### Start: determine whether u-cup-v gives a new relation
            ####
            if (Lett1 in CupBase6Lett) = false then
               
                if RelReducer.addIfIndependent(RelReducer.encode(solrel)) then
                    Append(CupRelsLett,[Lett1]);
                    Append(CupRel6Lett,[solrel]);
                fi;
            fi;
            ####
            #### End: determine whether u-cup-v gives a new relation
        
            iv := iv+1;
        od;
    fi;
    iu := iu+1;
od;
#
#
#Step-2 finished: (degree-2 cup with degree-2) cup with degree-2-gen, or degree-4 cup with degree-2-gen


#Step-3 begins here: degree-3-gen cup with degree-3-gen
#
#
iu :=1;
for u in Gen3 do

    #### implementing Mod2CupProduct(R,u,v,3,3,CB[3],CB[3],CB[6]) -- part 1: ####
    ####
    uCocycle:=CB[3].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,3,3);
    ww:=[];
    for i in [1..(R!.dimension(6))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,3,3,CB[3],CB[3],CB[6])  -- part 1 ended ####

    iv :=1;
    for v in Gen3 do
        if iv>=iu then
            Lett1 := GensLett[Length(Gen1)+Length(Gen2)+iu] + GensLett[Length(Gen1)+Length(Gen2)+iv];
            #Of course now u-cup-v does not contain patterns in CupRelsLett:
        
            
            #### implementing Mod2CupProduct(R,u,v,3,3,CB[3],CB[3],CB[6]) -- part 2: ####
            ####
            vCocycle:=CB[3].classToCocycle(v);
            uvCocycle:=[];
            for i in [1..(R!.dimension(6))] do
                w:=ww[i];
                sw:=0;
                for x in w do
                    sw:=sw + vCocycle[AbsoluteValue(x[1])];
                od;
                uvCocycle[i]:=sw mod 2;
            od;
            cuppedRaw := uvCocycle;
            ####
            #### implementing Mod2CupProduct(R,u,v,3,3,CB[3],CB[3],CB[6]) -- part 2 ended ####
            
            #### Start:  determine whether u-cup-v goes to the basis
            ####
            if (SGC_CachedSolveIn(SolverCache,BasisImaged2,cuppedRaw) = fail) = false then        #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel6Lett,[[Lett1]]);
            
            elif CupBase6Raw = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
                Append(CupBase6Raw,[cuppedRaw]);
                Append(CupBase6Lett,[Lett1]);
                Append(CupBase6RawCobandCoc,[cuppedRaw]);
            else
                sol :=SGC_CachedSolveIn(SolverCache,CupBase6RawCobandCoc,cuppedRaw);
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase6Raw,[cuppedRaw]);
                    Append(CupBase6Lett,[Lett1]);
                    Append(CupBase6RawCobandCoc,[cuppedRaw]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    sol1:=List([(Length(BasisImaged2)+1)..(Length(sol))],x->sol[x]);
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol1)),x->CupBase6Lett[x]));
                    if (Lett1 in CupBase6Lett) = false then
                        Append(CupRelsLett,[Lett1]);
                        Append(CupRel6Lett,[solrel]);
                    fi;
                fi;
            fi;
            ####
            #### End:  determine whether u-cup-v goes to the basis
            
        fi;
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-3 finished: degree-3-gen cup with degree-3-gen

fi;
#
#
#
#         ##Above ends the "Raw Method" for the degree 6 relations for groups 221-230, part2

CupRelByDeg[6] := CupRel6Lett;

##Begin printing the relations at deg 6:
if Length(CupRel6Lett) >0 then
fi;
##End printing the relations at deg 6.

####################################################################################
#The part of the code below only runs when relations at r=7 and r=8 are needed, which requires resolution at degree 9 constructed:
#
#
#
#The "if" below is the one that starts working on the r=7 and r=8 relations, which is executed only when Gen4 is non-empty!

if n >= 8 then

####################### r = 7 ##########################

CupBase7 :=[];
CupBase7Lett :=[];
CupRel7Lett := [];
CupTemp := [];
CupTempLett := [];
RelReducer := SGC_NewRelationReducer(GensLett,GenDegAll,CupRelByDeg,7);
if IsBoundGlobal("SGC_RELATION_REDUCER_OBSERVER") then
    ValueGlobal("SGC_RELATION_REDUCER_OBSERVER")(rec(degree := RelReducer.degree, relationProductCount := RelReducer.relationProductCount, rank := RelReducer.rank));
fi;


#Step-1 begins here: degree-6 cup with degree-1-gen
#
#
iu :=1;
for u in CupBase6 do

    #### implementing Mod2CupProduct(R,u,v,6,1,CB[6],CB[1],CB[7]) -- part 1: ####
    ####
    uCocycle:=CB[6].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,6,1);
    ww:=[];
    for i in [1..(R!.dimension(7))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,6,1,CB[6],CB[1],CB[7])  -- part 1 ended ####

    iv :=1;
    for v in Gen1 do
    
        
        #### implementing Mod2CupProduct(R,u,v,6,1,CB[6],CB[1],CB[7]) -- part 2: ####
        ####
        vCocycle:=CB[1].classToCocycle(v);
        uvCocycle:=[];
        for i in [1..(R!.dimension(7))] do
            w:=ww[i];
            sw:=0;
            for x in w do
                sw:=sw + vCocycle[AbsoluteValue(x[1])];
            od;
            uvCocycle[i]:=sw mod 2;
        od;
        cupped := CB[7].cocycleToClass(uvCocycle);
        ####
        #### implementing Mod2CupProduct(R,u,v,6,1,CB[6],CB[1],CB[7]) -- part 2 ended ####

        Lett1 := CupBase6Lett[iu] + GensLett[iv];

        
        #### Start: determine whether u-cup-v goes to the basis
        ####
        solrel := [];
        if cupped = zeroH[7] then         #if u-cup-v is a coboundary
            solrel := [Lett1];
            
        else
            if CupBase7 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
                Append(CupBase7,[cupped]);
                Append(CupBase7Lett,[Lett1]);
            else
                sol :=SGC_CachedSolveIn(SolverCache,CupBase7,cupped);
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase7,[cupped]);
                    Append(CupBase7Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase7Lett[x]));
                fi;
            fi;
        fi;
        ####
        #### End:  determine whether u-cup-v goes to the basis
            
            
        #### Start: determine whether u-cup-v gives a new relation
        ####
        if (Lett1 in CupBase7Lett) = false then
                    
            if RelReducer.addIfIndependent(RelReducer.encode(solrel)) then
                Append(CupRelsLett,[Lett1]);
                Append(CupRel7Lett,[solrel]);
            fi;
        fi;
        ####
        #### End: determine whether u-cup-v gives a new relation

    
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-1 finished: degree-6 cup with degree-1-gen


#Step-2 begins here: (degree-3 cup with degree-2) cup with degree-2-gen
#
#
iu :=1;
for u in CupBase5 do
    if List([1..Length(Gen1)],x->CupBase5Lett[iu][x]) = List([1..Length(Gen1)],x->0) then   #if u is of the form Bi-cup-Cj
        #### implementing Mod2CupProduct(R,u,v,5,2,CB[5],CB[2],CB[7]) -- part 1: ####
        ####
        uCocycle:=CB[5].classToCocycle(u);
        uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,5,2);
        ww:=[];
        for i in [1..(R!.dimension(7))] do
            Append(ww, [uChainMap([[i,1]])]);
        od;
        ####
        #### implementing Mod2CupProduct(R,u,v,5,2,CB[5],CB[2],CB[7])  -- part 1 ended ####

        iv :=1;
        for v in Gen2 do
    
        
            #### implementing Mod2CupProduct(R,u,v,5,2,CB[5],CB[2],CB[7]) -- part 2: ####
            ####
            vCocycle:=CB[2].classToCocycle(v);
            uvCocycle:=[];
            for i in [1..(R!.dimension(7))] do
                w:=ww[i];
                sw:=0;
                for x in w do
                    sw:=sw + vCocycle[AbsoluteValue(x[1])];
                od;
                uvCocycle[i]:=sw mod 2;
            od;
            cupped := CB[7].cocycleToClass(uvCocycle);
            ####
            #### implementing Mod2CupProduct(R,u,v,5,2,CB[5],CB[2],CB[7]) -- part 2 ended ####
        

            Lett1 := CupBase5Lett[iu] + GensLett[Length(Gen1)+iv];

        
            #### Start: determine whether u-cup-v goes to the basis
            ####
            solrel := [];
            if cupped = zeroH[7] then         #if u-cup-v is a coboundary
                solrel := [Lett1];
              
            else
                if CupBase7 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
                    Append(CupBase7,[cupped]);
                    Append(CupBase7Lett,[Lett1]);
                else
                    sol :=SGC_CachedSolveIn(SolverCache,CupBase7,cupped);
                    if sol = fail then                  #if u-cup-v is a genuine new cocycle
                        Append(CupBase7,[cupped]);
                        Append(CupBase7Lett,[Lett1]);
                    else                                #if u-cup-v is expressable by other cocycles in basis
                        solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase7Lett[x]));
                    fi;
                fi;
            fi;
            ####
            #### End:  determine whether u-cup-v goes to the basis
            
            
            #### Start: determine whether u-cup-v gives a new relation
            ####
            if (Lett1 in CupBase7Lett) = false then
               
                if RelReducer.addIfIndependent(RelReducer.encode(solrel)) then
                    Append(CupRelsLett,[Lett1]);
                    Append(CupRel7Lett,[solrel]);
                fi;
            fi;
            ####
            #### End: determine whether u-cup-v gives a new relation
        
            iv := iv+1;
        od;
    fi;
    iu := iu+1;
od;
#
#
#Step-2 finished: (degree-3 cup with degree-2) cup with degree-2-gen


#Step-3 begins here: degree-4-gen cup with degree-3-gen
#
#
iu :=1;
for u in Gen4 do

    #### implementing Mod2CupProduct(R,u,v,4,3,CB[4],CB[3],CB[7]) -- part 1: ####
    ####
    uCocycle:=CB[4].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,4,3);
    ww:=[];
    for i in [1..(R!.dimension(7))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,4,3,CB[4],CB[3],CB[7])  -- part 1 ended ####

    iv :=1;
    for v in Gen3 do
        if iv>=iu then
            Lett1 := GensLett[Length(Gen1)+Length(Gen2)+Length(Gen3)+iu] + GensLett[Length(Gen1)+Length(Gen2)+iv];
            #Of course now u-cup-v does not contain patterns in CupRelsLett:
        
            
            #### implementing Mod2CupProduct(R,u,v,4,3,CB[4],CB[3],CB[7]) -- part 2: ####
            ####
            vCocycle:=CB[3].classToCocycle(v);
            uvCocycle:=[];
            for i in [1..(R!.dimension(7))] do
                w:=ww[i];
                sw:=0;
                for x in w do
                    sw:=sw + vCocycle[AbsoluteValue(x[1])];
                od;
                uvCocycle[i]:=sw mod 2;
            od;
            cupped := CB[7].cocycleToClass(uvCocycle);
            ####
            #### implementing Mod2CupProduct(R,u,v,4,3,CB[4],CB[3],CB[7]) -- part 2 ended ####
        
            if cupped = zeroH[7] then         #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel7Lett,[[Lett1]]);
            elif CupBase7 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
                Append(CupBase7,[cupped]);
                Append(CupBase7Lett,[Lett1]);
            else
                sol :=SGC_CachedSolveIn(SolverCache,CupBase7,cupped);
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase7,[cupped]);
                    Append(CupBase7Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase7Lett[x]));
                    if (Lett1 in CupBase7Lett) = false then
                        Append(CupRelsLett,[Lett1]);
                        Append(CupRel7Lett,[solrel]);
                    fi;
                fi;
            fi;
        fi;
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-3 finished: degree-4-gen cup with degree-3-gen

#Note that there's no Gen7!!


if Length(CupBase7) = CB[7].Mod2Cohomologydim then
    Print("");#Print("dim(H^7)=", CB[7].Mod2Cohomologydim,".\n");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^7) = ", Length(CupBase7) - CB[7].Mod2Cohomologydim,"\n");
fi;
CupRelByDeg[7] := CupRel7Lett;


####################### r = 8 ##########################

CupBase8 :=[];
CupBase8Lett :=[];
CupRel8Lett := [];
CupTemp := [];
CupTempLett := [];
RelReducer := SGC_NewRelationReducer(GensLett,GenDegAll,CupRelByDeg,8);
if IsBoundGlobal("SGC_RELATION_REDUCER_OBSERVER") then
    ValueGlobal("SGC_RELATION_REDUCER_OBSERVER")(rec(degree := RelReducer.degree, relationProductCount := RelReducer.relationProductCount, rank := RelReducer.rank));
fi;


#Step-1 begins here: degree-7 cup with degree-1-gen
#
#
iu :=1;
for u in CupBase7 do

    #### implementing Mod2CupProduct(R,u,v,7,1,CB[7],CB[1],CB[8]) -- part 1: ####
    ####
    uCocycle:=CB[7].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,7,1);
    ww:=[];
    for i in [1..(R!.dimension(8))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,7,1,CB[7],CB[1],CB[8])  -- part 1 ended ####

    iv :=1;
    for v in Gen1 do
    
        
        #### implementing Mod2CupProduct(R,u,v,7,1,CB[7],CB[1],CB[8]) -- part 2: ####
        ####
        vCocycle:=CB[1].classToCocycle(v);
        uvCocycle:=[];
        for i in [1..(R!.dimension(8))] do
            w:=ww[i];
            sw:=0;
            for x in w do
                sw:=sw + vCocycle[AbsoluteValue(x[1])];
            od;
            uvCocycle[i]:=sw mod 2;
        od;
        cupped := CB[8].cocycleToClass(uvCocycle);
        ####
        #### implementing Mod2CupProduct(R,u,v,7,1,CB[7],CB[1],CB[8]) -- part 2 ended ####

        Lett1 := CupBase7Lett[iu] + GensLett[iv];

        
        #### Start: determine whether u-cup-v goes to the basis
        ####
        solrel := [];
        if cupped = zeroH[8] then         #if u-cup-v is a coboundary
            solrel := [Lett1];
            
        else
            if CupBase8 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
                Append(CupBase8,[cupped]);
                Append(CupBase8Lett,[Lett1]);
            else
                sol :=SGC_CachedSolveIn(SolverCache,CupBase8,cupped);
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase8,[cupped]);
                    Append(CupBase8Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase8Lett[x]));
                fi;
            fi;
        fi;
        ####
        #### End:  determine whether u-cup-v goes to the basis
            
            
        #### Start: determine whether u-cup-v gives a new relation
        ####
        if (Lett1 in CupBase8Lett) = false then
                    
            if RelReducer.addIfIndependent(RelReducer.encode(solrel)) then
                Append(CupRelsLett,[Lett1]);
                Append(CupRel8Lett,[solrel]);
            fi;
        fi;
        ####
        #### End: determine whether u-cup-v gives a new relation

    
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-1 finished: degree-7 cup with degree-1-gen


#Step-2 begins here: (degree-3 cup with degree-3, or degree-2 cup with degree-2 cup with degree-2) cup with degree-2-gen
#
#
iu :=1;
for u in CupBase6 do
    if List([1..Length(Gen1)],x->CupBase6Lett[iu][x]) = List([1..Length(Gen1)],x->0) then   #if u is of the form Ci-cup-Cj or Bi-cup-Bj-cup-Bk
        #### implementing Mod2CupProduct(R,u,v,6,2,CB[6],CB[2],CB[8]) -- part 1: ####
        ####
        uCocycle:=CB[6].classToCocycle(u);
        uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,6,2);
        ww:=[];
        for i in [1..(R!.dimension(8))] do
            Append(ww, [uChainMap([[i,1]])]);
        od;
        ####
        #### implementing Mod2CupProduct(R,u,v,6,2,CB[6],CB[2],CB[8])  -- part 1 ended ####

        iv :=1;
        for v in Gen2 do
    
        
            #### implementing Mod2CupProduct(R,u,v,6,2,CB[6],CB[2],CB[8]) -- part 2: ####
            ####
            vCocycle:=CB[2].classToCocycle(v);
            uvCocycle:=[];
            for i in [1..(R!.dimension(8))] do
                w:=ww[i];
                sw:=0;
                for x in w do
                    sw:=sw + vCocycle[AbsoluteValue(x[1])];
                od;
                uvCocycle[i]:=sw mod 2;
            od;
            cupped := CB[8].cocycleToClass(uvCocycle);
            ####
            #### implementing Mod2CupProduct(R,u,v,6,2,CB[6],CB[2],CB[8]) -- part 2 ended ####
        

            Lett1 := CupBase6Lett[iu] + GensLett[Length(Gen1)+iv];

        
            #### Start: determine whether u-cup-v goes to the basis
            ####
            solrel := [];
            if cupped = zeroH[8] then         #if u-cup-v is a coboundary
                solrel := [Lett1];
              
            else
                if CupBase8 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
                    Append(CupBase8,[cupped]);
                    Append(CupBase8Lett,[Lett1]);
                else
                    sol :=SGC_CachedSolveIn(SolverCache,CupBase8,cupped);
                    if sol = fail then                  #if u-cup-v is a genuine new cocycle
                        Append(CupBase8,[cupped]);
                        Append(CupBase8Lett,[Lett1]);
                    else                                #if u-cup-v is expressable by other cocycles in basis
                        solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase8Lett[x]));
                    fi;
                fi;
            fi;
            ####
            #### End:  determine whether u-cup-v goes to the basis
            
            
            #### Start: determine whether u-cup-v gives a new relation
            ####
            if (Lett1 in CupBase8Lett) = false then
               
                if RelReducer.addIfIndependent(RelReducer.encode(solrel)) then
                    Append(CupRelsLett,[Lett1]);
                    Append(CupRel8Lett,[solrel]);
                fi;
            fi;
            ####
            #### End: determine whether u-cup-v gives a new relation
        
            iv := iv+1;
        od;
    fi;
    iu := iu+1;
od;
#
#
#Step-2 finished: (degree-3 cup with degree-3, degree-2 cup with degree-2 cup with degree-2) cup with degree-2-gen


#Step-3 begins here: degree-4-gen cup with degree-4-gen
#
#
iu :=1;
for u in Gen4 do

    #### implementing Mod2CupProduct(R,u,v,4,4,CB[4],CB[4],CB[8]) -- part 1: ####
    ####
    uCocycle:=CB[4].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,4,4);
    ww:=[];
    for i in [1..(R!.dimension(8))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,4,4,CB[4],CB[4],CB[8])  -- part 1 ended ####

    iv :=1;
    for v in Gen4 do
        if iv>=iu then
            Lett1 := GensLett[Length(Gen1)+Length(Gen2)+Length(Gen3)+iu] + GensLett[Length(Gen1)+Length(Gen2)+Length(Gen3)+iv];
            #Of course now u-cup-v does not contain patterns in CupRelsLett:
        
            
            #### implementing Mod2CupProduct(R,u,v,4,4,CB[4],CB[4],CB[8]) -- part 2: ####
            ####
            vCocycle:=CB[4].classToCocycle(v);
            uvCocycle:=[];
            for i in [1..(R!.dimension(8))] do
                w:=ww[i];
                sw:=0;
                for x in w do
                    sw:=sw + vCocycle[AbsoluteValue(x[1])];
                od;
                uvCocycle[i]:=sw mod 2;
            od;
            cupped := CB[8].cocycleToClass(uvCocycle);
            ####
            #### implementing Mod2CupProduct(R,u,v,4,4,CB[4],CB[4],CB[8]) -- part 2 ended ####
        
            if cupped = zeroH[8] then         #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel8Lett,[[Lett1]]);
            elif CupBase8 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
                Append(CupBase8,[cupped]);
                Append(CupBase8Lett,[Lett1]);
            else
                sol :=SGC_CachedSolveIn(SolverCache,CupBase8,cupped);
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase8,[cupped]);
                    Append(CupBase8Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase8Lett[x]));
                    if (Lett1 in CupBase8Lett) = false then
                        Append(CupRelsLett,[Lett1]);
                        Append(CupRel8Lett,[solrel]);
                    fi;
                fi;
            fi;
        fi;
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step-3 finished: degree-4-gen cup with degree-4-gen

#Note that there's no Gen8!!


if Length(CupBase8) = CB[8].Mod2Cohomologydim then
    Print("");#Print("dim(H^8)=", CB[8].Mod2Cohomologydim,".\n");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^8) = ", Length(CupBase8) - CB[8].Mod2Cohomologydim,"\n");
fi;
CupRelByDeg[8] := CupRel8Lett;


fi;
#the above "fi;" is the one that ends working on the r=7 and r=8 relations, which is executed
#whenever the resolution has length >= 9 (degree-4 and degree-6 generator groups).
#
#
#####################################################################################
#
#Relations at degree 9..12 exist only for the degree-6 generator groups
#(SGC_Degree6GeneratorGroups = [219, 226, 228]), whose resolutions are built to
#length 13; every other group never reaches this block. For these groups
#Gen4 = Gen5 = [], so of the generator degrees {1,2,3,6} each degree-r monomial
#is reached by peeling off its lowest-degree generator factor: Step 1 cups
#CupBase(r-1) with Gen1, Step 2 cups CupBase(r-2) without Gen1 factors with
#Gen2, Step 3 cups CupBase(r-3) without Gen1/Gen2 factors with Gen3, and the
#only monomial with Gen6 factors alone is Gen6*Gen6 at r=12 (Step 4). A Step-3
#source basis restricted to Gen3/Gen6 factors is non-empty only when its degree
#r-3 is a multiple of 3, i.e. only at r=9 and r=12; r=10 and r=11 have no Step 3.
#The per-step logic reproduces the r=7/r=8 Steps verbatim, parameterized by
#degrees and offsets instead of hand-copied per degree split.

if n >= 12 then

#One relation-search Step (the r=7/r=8 Step-1/2/3 pattern, parameterized):
#cup every element of CupBaseA (degree dega, restricted so its first
#restrictWidth letter components vanish; restrictWidth 0 means unrestricted)
#with every generator in GenB (degree degb, letters at GensLett[offB+..]),
#extend CupBaseR/CupBaseRLett (degree degr = dega+degb) with genuinely new
#classes, and record relations that are not reducible from lower degrees.
CupStep := function(CupBaseA, CupBaseALett, dega, GenB, offB, degb, restrictWidth, CupBaseR, CupBaseRLett, CupRelRLett, degr)
local iu, iv, i, u, v, w, x, sw, uCocycle, uChainMap, ww, vCocycle, uvCocycle,
      cupped, Lett1, solrel, sol;
iu := 1;
for u in CupBaseA do
    if restrictWidth = 0 or List([1..restrictWidth],x->CupBaseALett[iu][x]) = List([1..restrictWidth],x->0) then

        #### implementing Mod2CupProduct(R,u,v,dega,degb,CB[dega],CB[degb],CB[degr]) -- part 1: ####
        uCocycle:=CB[dega].classToCocycle(u);
        uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,dega,degb);
        ww:=[];
        for i in [1..(R!.dimension(degr))] do
            Append(ww, [uChainMap([[i,1]])]);
        od;

        iv := 1;
        for v in GenB do

            #### implementing Mod2CupProduct -- part 2: ####
            vCocycle:=CB[degb].classToCocycle(v);
            uvCocycle:=[];
            for i in [1..(R!.dimension(degr))] do
                w:=ww[i];
                sw:=0;
                for x in w do
                    sw:=sw + vCocycle[AbsoluteValue(x[1])];
                od;
                uvCocycle[i]:=sw mod 2;
            od;
            cupped := CB[degr].cocycleToClass(uvCocycle);

            Lett1 := CupBaseALett[iu] + GensLett[offB+iv];

            #### Start: determine whether u-cup-v goes to the basis
            solrel := [];
            if cupped = zeroH[degr] then         #if u-cup-v is a coboundary
                solrel := [Lett1];
            else
                if CupBaseR = [] then        #if no basis yet then push in the genuine cocycle u-cup-v to basis
                    Append(CupBaseR,[cupped]);
                    Append(CupBaseRLett,[Lett1]);
                else
                    sol :=SGC_CachedSolveIn(SolverCache,CupBaseR,cupped);
                    if sol = fail then                  #if u-cup-v is a genuine new cocycle
                        Append(CupBaseR,[cupped]);
                        Append(CupBaseRLett,[Lett1]);
                    else                                #if u-cup-v is expressable by other cocycles in basis
                        solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBaseRLett[x]));
                    fi;
                fi;
            fi;
            #### End:  determine whether u-cup-v goes to the basis

            #### Start: determine whether u-cup-v gives a new relation
            if (Lett1 in CupBaseRLett) = false then
                if RelReducer.addIfIndependent(RelReducer.encode(solrel)) then
                    Append(CupRelsLett,[Lett1]);
                    Append(CupRelRLett,[solrel]);
                fi;
            fi;
            #### End: determine whether u-cup-v gives a new relation

            iv := iv+1;
        od;
    fi;
    iu := iu+1;
od;
end;

####################### r = 9 ##########################

CupBase9 :=[];
CupBase9Lett :=[];
CupRel9Lett := [];
RelReducer := SGC_NewRelationReducer(GensLett,GenDegAll,CupRelByDeg,9);
if IsBoundGlobal("SGC_RELATION_REDUCER_OBSERVER") then
    ValueGlobal("SGC_RELATION_REDUCER_OBSERVER")(rec(degree := RelReducer.degree, relationProductCount := RelReducer.relationProductCount, rank := RelReducer.rank));
fi;
#Step 1: degree-8 basis cup with degree-1 generators
CupStep(CupBase8, CupBase8Lett, 8, Gen1, 0, 1, 0, CupBase9, CupBase9Lett, CupRel9Lett, 9);
#Step 2: degree-7 basis without Gen1 factors cup with degree-2 generators
CupStep(CupBase7, CupBase7Lett, 7, Gen2, Length(Gen1), 2, Length(Gen1), CupBase9, CupBase9Lett, CupRel9Lett, 9);
#Step 3: degree-6 basis without Gen1/Gen2 factors (raw Gen6 and Gen3*Gen3 monomials) cup with degree-3 generators
CupStep(CupBase6, CupBase6Lett, 6, Gen3, Length(Gen1)+Length(Gen2), 3, Length(Gen1)+Length(Gen2), CupBase9, CupBase9Lett, CupRel9Lett, 9);

#Note that there's no Gen9!!

if Length(CupBase9) = CB[9].Mod2Cohomologydim then
    Print("");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^9) = ", Length(CupBase9) - CB[9].Mod2Cohomologydim,"\n");
fi;
CupRelByDeg[9] := CupRel9Lett;

####################### r = 10 #########################

CupBase10 :=[];
CupBase10Lett :=[];
CupRel10Lett := [];
RelReducer := SGC_NewRelationReducer(GensLett,GenDegAll,CupRelByDeg,10);
if IsBoundGlobal("SGC_RELATION_REDUCER_OBSERVER") then
    ValueGlobal("SGC_RELATION_REDUCER_OBSERVER")(rec(degree := RelReducer.degree, relationProductCount := RelReducer.relationProductCount, rank := RelReducer.rank));
fi;
#Step 1: degree-9 basis cup with degree-1 generators
CupStep(CupBase9, CupBase9Lett, 9, Gen1, 0, 1, 0, CupBase10, CupBase10Lett, CupRel10Lett, 10);
#Step 2: degree-8 basis without Gen1 factors cup with degree-2 generators
CupStep(CupBase8, CupBase8Lett, 8, Gen2, Length(Gen1), 2, Length(Gen1), CupBase10, CupBase10Lett, CupRel10Lett, 10);
#No Step 3: degree 7 is not a multiple of 3, so no Gen3/Gen6-only monomial has degree 7.

if Length(CupBase10) = CB[10].Mod2Cohomologydim then
    Print("");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^10) = ", Length(CupBase10) - CB[10].Mod2Cohomologydim,"\n");
fi;
CupRelByDeg[10] := CupRel10Lett;

####################### r = 11 #########################

CupBase11 :=[];
CupBase11Lett :=[];
CupRel11Lett := [];
RelReducer := SGC_NewRelationReducer(GensLett,GenDegAll,CupRelByDeg,11);
if IsBoundGlobal("SGC_RELATION_REDUCER_OBSERVER") then
    ValueGlobal("SGC_RELATION_REDUCER_OBSERVER")(rec(degree := RelReducer.degree, relationProductCount := RelReducer.relationProductCount, rank := RelReducer.rank));
fi;
#Step 1: degree-10 basis cup with degree-1 generators
CupStep(CupBase10, CupBase10Lett, 10, Gen1, 0, 1, 0, CupBase11, CupBase11Lett, CupRel11Lett, 11);
#Step 2: degree-9 basis without Gen1 factors cup with degree-2 generators
CupStep(CupBase9, CupBase9Lett, 9, Gen2, Length(Gen1), 2, Length(Gen1), CupBase11, CupBase11Lett, CupRel11Lett, 11);
#No Step 3: degree 8 is not a multiple of 3, so no Gen3/Gen6-only monomial has degree 8.

if Length(CupBase11) = CB[11].Mod2Cohomologydim then
    Print("");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^11) = ", Length(CupBase11) - CB[11].Mod2Cohomologydim,"\n");
fi;
CupRelByDeg[11] := CupRel11Lett;

####################### r = 12 #########################

CupBase12 :=[];
CupBase12Lett :=[];
CupRel12Lett := [];
RelReducer := SGC_NewRelationReducer(GensLett,GenDegAll,CupRelByDeg,12);
if IsBoundGlobal("SGC_RELATION_REDUCER_OBSERVER") then
    ValueGlobal("SGC_RELATION_REDUCER_OBSERVER")(rec(degree := RelReducer.degree, relationProductCount := RelReducer.relationProductCount, rank := RelReducer.rank));
fi;
#Step 1: degree-11 basis cup with degree-1 generators
CupStep(CupBase11, CupBase11Lett, 11, Gen1, 0, 1, 0, CupBase12, CupBase12Lett, CupRel12Lett, 12);
#Step 2: degree-10 basis without Gen1 factors cup with degree-2 generators
CupStep(CupBase10, CupBase10Lett, 10, Gen2, Length(Gen1), 2, Length(Gen1), CupBase12, CupBase12Lett, CupRel12Lett, 12);
#Step 3: degree-9 basis without Gen1/Gen2 factors cup with degree-3 generators
CupStep(CupBase9, CupBase9Lett, 9, Gen3, Length(Gen1)+Length(Gen2), 3, Length(Gen1)+Length(Gen2), CupBase12, CupBase12Lett, CupRel12Lett, 12);

#Step 4 begins here: degree-6-gen cup with degree-6-gen (mirrors the r=8 Gen4-cup-Gen4 Step)
#
#
iu :=1;
for u in Gen6 do

    #### implementing Mod2CupProduct(R,u,v,6,6,CB[6],CB[6],CB[12]) -- part 1: ####
    ####
    uCocycle:=CB[6].classToCocycle(u);
    uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,6,6);
    ww:=[];
    for i in [1..(R!.dimension(12))] do
        Append(ww, [uChainMap([[i,1]])]);
    od;
    ####
    #### implementing Mod2CupProduct(R,u,v,6,6,CB[6],CB[6],CB[12])  -- part 1 ended ####

    iv :=1;
    for v in Gen6 do
        if iv>=iu then
            Lett1 := GensLett[Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5)+iu] + GensLett[Length(Gen1)+Length(Gen2)+Length(Gen3)+Length(Gen4)+Length(Gen5)+iv];
            #Of course now u-cup-v does not contain patterns in CupRelsLett:


            #### implementing Mod2CupProduct(R,u,v,6,6,CB[6],CB[6],CB[12]) -- part 2: ####
            ####
            vCocycle:=CB[6].classToCocycle(v);
            uvCocycle:=[];
            for i in [1..(R!.dimension(12))] do
                w:=ww[i];
                sw:=0;
                for x in w do
                    sw:=sw + vCocycle[AbsoluteValue(x[1])];
                od;
                uvCocycle[i]:=sw mod 2;
            od;
            cupped := CB[12].cocycleToClass(uvCocycle);
            ####
            #### implementing Mod2CupProduct(R,u,v,6,6,CB[6],CB[6],CB[12]) -- part 2 ended ####

            if cupped = zeroH[12] then         #if u-cup-v is a coboundary
                Append(CupRelsLett,[Lett1]);
                Append(CupRel12Lett,[[Lett1]]);
            elif CupBase12 = [] then        #if no basis yet then push in the genuine cocycle u-cup-v
                Append(CupBase12,[cupped]);
                Append(CupBase12Lett,[Lett1]);
            else
                sol :=SGC_CachedSolveIn(SolverCache,CupBase12,cupped);
                if sol = fail then                  #if u-cup-v is a genuine new cocycle
                    Append(CupBase12,[cupped]);
                    Append(CupBase12Lett,[Lett1]);
                else                                #if u-cup-v is expressable by other cocycles in basis
                    solrel:=Concatenation([Lett1], List(IToPosition(GF2ToZ(sol)),x->CupBase12Lett[x]));
                    if (Lett1 in CupBase12Lett) = false then
                        Append(CupRelsLett,[Lett1]);
                        Append(CupRel12Lett,[solrel]);
                    fi;
                fi;
            fi;
        fi;
        iv := iv+1;
    od;
    iu := iu+1;
od;
#
#
#Step 4 finished: degree-6-gen cup with degree-6-gen

#Note that there's no Gen12!!

if Length(CupBase12) = CB[12].Mod2Cohomologydim then
    Print("");
else
    Print("!!!! No match!!!! dim(Chosen basis) - dim(H^12) = ", Length(CupBase12) - CB[12].Mod2Cohomologydim,"\n");
fi;
CupRelByDeg[12] := CupRel12Lett;

fi;
#the above "fi;" ends the r = 9..12 relation search for the degree-6 generator groups.
#
#
#
#
#####################################################################################

#
####   Begin printing cohomology ring   ####
#
#CupRelByDeg was populated incrementally as each relation degree completed.

relations := [];
for r in [2..12] do
    if IsBound(CupRelByDeg[r]) and Length(CupRelByDeg[r]) > 0 then
        Add(relations,[r,CupRelByDeg[r]]);
    fi;
od;

if Length(arg) >= 7 then
    generatorNames := arg[7];
else
    generatorNames := GENNAMES[IT];
fi;
if Length(generatorNames) <> Sum(GenDimAll) then
    Error("Mod2RingGensAndRels: generator-name count does not match generator count for IT=", IT, "\n");
fi;

ringData := rec(
    IT := IT,
    generators := rec(
        names := ShallowCopy(generatorNames),
        degrees := ShallowCopy(GenDegAll),
        classes := [Gen1,Gen2,Gen3,Gen4,Gen5,Gen6]
    ),
    relations := relations,
    relationsByDegree := CupRelByDeg,
    bases := [CupBase1Lett,CupBase2Lett,CupBase3Lett,CupBase4Lett],
    generatorDimensions := ShallowCopy(GenDimAll),
    generatorDegrees := ShallowCopy(GenDegAll),
    maxRelationDegree := n
);

if not returnData then
    mono:=List(List([1..Sum(GenDimAll)],x->GensLett[x]),x->Letter2Monomial(x,GenDimAll,generatorNames));
    Print("Z2[", JoinStringsWithSeparator(mono,","), "]");

    mono := 0;
    for r in [2..12] do
        if IsBound(CupRelByDeg[r]) and Length(CupRelByDeg[r]) > 0 then
            if mono = 0 then           #mono = 1 means at least one relation already printed
                Print("/<R", r);
            else
                Print(",R", r);
            fi;
            mono := 1;
        fi;
    od;
    if mono = 1 then
        Print(">\n");
    else
        Print("\n");
    fi;

    for r in [2..12] do
        if IsBound(CupRelByDeg[r]) and Length(CupRelByDeg[r]) > 0 then
            Print(Concatenation("R", String(r), ":  "));
            List(CupRelByDeg[r],x->PrintMonomialString(x,GenDimAll,"+",generatorNames));
            Print("\n");
        fi;
    od;

fi;

####   End printing cohomology ring   ####
#

if returnData then
    return ringData;
fi;
return ringData.bases;
end;
#####################################################################
#####################################################################

SGC_PrintMod2RingData:=function(Data)
local mono, pair, letters;

letters:=CohomologyBasis(List([1..Length(Data.generators.names)],i->1));
mono:=List(letters,x->Letter2Monomial(x,Data.generatorDimensions,Data.generators.names));
Print("Z2[", JoinStringsWithSeparator(mono,","), "]");
if Length(Data.relations) > 0 then
    Print("/<",JoinStringsWithSeparator(List(Data.relations,x->Concatenation("R",String(x[1]))),","),">\n");
else
    Print("\n");
fi;
for pair in Data.relations do
    Print(Concatenation("R", String(pair[1]), ":  "));
    List(pair[2],x->PrintMonomialString(x,Data.generatorDimensions,"+",Data.generators.names));
    Print("\n");
od;
end;

#####################################################################
#####################################################################


#####################################################################
#####################################################################

PointGroupTranslationExtension:=function(arg)
local
    Gs,arithmeticNo,ZZZ,Gpt,R,C,Homo,
    T1,T2,T3,elem,eleml,conjT1,conjT2,conjT3,
    i,j,k,l,rep,rep1,flag;
    

Gs:=[Group([(1,2)]),Group([(1,2),(3,4)]),Group([(1,2),(3,4),(5,6)]),Group([(1,3)(2,4),(1,2,3,4)]),Group([(1,3)(2,4),(1,2,3,4),(5,6)]),Group([(1,3)(2,4),(1,3)(5,6),(1,4)(2,3)(5,6)]),Group([(1,3)(2,4),(1,3)(5,6),(1,4)(2,3)(5,6),(7,8)]),Group([(1,2,3)]),Group([(1,2,3),(4,5)]),Group([(1,2,3),(2,3)(4,5)]),Group([(1,2,3),(2,3)(4,5),(6,7)]),Group([(1,2,3),(4,5),(6,7)]),Group([(1,2,3),(2,3)(4,5),(6,7),(8,9)]),Group([(1,2)(3,4),(1,3)(2,4),(3,2,1)]),Group([(1,2)(3,4),(1,3)(2,4),(3,2,1),(5,6)]),Group([(1,2)(3,4),(1,3)(2,4),(3,2,1),(1,2)]),Group([(1,2)(3,4),(1,3)(2,4),(3,2,1),(1,2),(5,6)])];

arithmeticNo:=[[[2],[3,4],[5],[6,7],[8,9]],[[10,11,13,14],[12,15],[16,17,18,19],[20,21],[22],[23,24],[25,26,27,28,29,30,31,32,33,34],[35,36,37],[38,39,40,41],[42,43],[44,45,46]],[[47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62],[63,64,65,66,67,68],[69,70],[71,72,73,74]],[[75,76,77,78],[79,80],[81],[82]],[[83,84,85,86],[87,88]],[[89,90,91,92,93,94,95,96],[97,98],[99,100,101,102,103,104,105,106],[107,108,109,110],[111,112,113,114],[115,116,117,118],[119,120],[121,122]],[[123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138],[139,140,141,142]],[[143,144,145],[146]],[[147],[148],[168,169,170,171,172,173],[174]],[[149],[150],[151],[152],[153],[154],[155],[156],[157],[158],[159],[160,161]],[[162,163],[164,165],[166,167],[177,178,179,180,181,182],[183,184,185,186],[187,188],[189,190]],[[175,176]],[[191,192,193,194]],[[195],[196],[197],[198],[199]],[[200,201],[202,203],[204],[205],[206]],[[207,208],[209,210],[211],[212,213],[214],[215],[216],[217],[218],[219],[220]],[[221,222,223,224],[225,226,227,228],[229,230]]];

ZZZ:=GL(3,Integers);;

T1:=[[1,0,0,1],[0,1,0,0],[0,0,1,0],[0,0,0,1]]; #standard translation T1
T2:=[[1,0,0,0],[0,1,0,1],[0,0,1,0],[0,0,0,1]]; #standard translation T2
T3:=[[1,0,0,0],[0,1,0,0],[0,0,1,1],[0,0,0,1]]; #standard translation T3

Print("Group Extension Information for the 73 Arithmetic Classes Carrying Representation rho:\n");
Print("H^2_rho(PG,Z^3):\n");


for i in [1..Length(arithmeticNo)] do           #Length(arithmeticNo) == 17
    Gpt:=Gs[i];
    R:=ResolutionFiniteGroup(Gpt,4);
    for j in [1..Length(arithmeticNo[i])] do
        flag := arithmeticNo[i][j][1];
        
        for k in arithmeticNo[i][j] do
            for l in [1..Length(PGGens230[k])] do
                elem := List([1..3],x->List([1..3],y->PGGens230[flag][l][x][y]));
                eleml := List([1..3],x->List([1..3],y->PGGens230[k][l][x][y]));
                if (elem = eleml) = false then
                    Print("point group generators not chosen consistently with others in the same arithmetic class!\n");
                fi;
            od;
        od;
        
        rep := [];
        for elem in PGGens230[flag] do
            conjT1 := elem * T1 * elem^(-1);
            conjT2 := elem * T2 * elem^(-1);
            conjT3 := elem * T3 * elem^(-1);
            Append(rep,[[[conjT1[1,4],conjT2[1,4],conjT3[1,4]],[conjT1[2,4],conjT2[2,4],conjT3[2,4]],[conjT1[3,4],conjT2[3,4],conjT3[3,4]]]]);
        od;
        
        
        Homo := GroupHomomorphismByImages(Gpt,ZZZ,GeneratorsOfGroup(Gpt),rep);
        C:=HomToIntegralModule(R,Homo);
        if Length(arithmeticNo[i][j]) = 1 then
            Print("No. ", arithmeticNo[i][j][1],": ", Cohomology(C,2),"  ");
        elif List([1..Length(arithmeticNo[i][j])],x->arithmeticNo[i][j][x]) = List([1..Length(arithmeticNo[i][j])],x->arithmeticNo[i][j][1]+x-1) then
            Print("No. ", arithmeticNo[i][j][1],"-", arithmeticNo[i][j][Length(arithmeticNo[i][j])],": ", Cohomology(C,2),"  ");
        else
            Print("No. ", JoinStringsWithSeparator(List([1..Length(arithmeticNo[i][j])],x->String(arithmeticNo[i][j][x])),"&"),": ", Cohomology(C,2),"  ");
        fi;
    od;
    Print("\n");
od;


return true;
end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
IrreducibleWyckoffPoints:=function(arg)
local
    IT, Dim, SG, Rec, IWPs, x,
    WyckoffPosRelations, WyckoffGraphRecord;
    
#####################################################################
WyckoffPosRelations := function( W )
    local S, T, d, len, gens, G, L, m, O, o, i, j, k, Si, Sj, index, lst;

    S := WyckoffSpaceGroup( W[1] );
    T := TranslationBasis( S );
    d := DimensionOfMatrixGroup( S ) - 1;
    len  := Length( W );
    gens := GeneratorsOfGroup( S );
    gens := Filtered( gens, g -> g{[1..d]}{[1..d]} <> IdentityMat( d ) );
    if IsAffineCrystGroupOnLeft( S ) then
        gens := List( gens, TransposedMat );
    fi;
    G := GroupByGenerators( gens, One( S ) );
    L := List( W, w -> rec( translation := WyckoffTranslation( w ),
                            basis       := WyckoffBasis( w ),
                            spaceGroup  := S ) );

    m := NullMat( len, len );
    for i in [1..len] do
        O := Orbit( G, L[i], ImageAffineSubspaceLattice );
        for j in [1..len] do
            Sj := WyckoffStabilizer( W[j] );
            Si := WyckoffStabilizer( W[i] );
            index := Size(Sj) / Size(Si);
            if Length(L[j].basis) < Length(L[i].basis) and IsInt(index) then
                lst := Filtered(O,o->IsSubspaceAffineSubspaceLattice(o,L[j]));
                m[j][i] := Length( lst );
            fi;
        od;
    od;

    for i in Reversed([1..Length(W)]) do
        for j in Reversed([1..i-1]) do
            if m[j][i]<>0 then
                for k in [1..j-1] do
                    if m[k][j]<>0 then m[k][i]:=0; fi;
                od;
            fi;
        od;
    od;

    return m;

end;

#############################################################################
##
#F  WyckoffGraphRecord( <lst> ) . . . . . . . Create record for Wyckoff graph
##
WyckoffGraphRecord := function( lst )

    local L, m, R, i, level, j;

    L := List( lst, w -> rec( wypos := w,
                              dim   := Length( WyckoffBasis(w) ),
                              size  := Size( WyckoffStabilizer(w) ),
                              class := w!.class ) );
    Sort( L, function(a,b) return a.size > b.size; end );

    m := WyckoffPosRelations( List( L, x -> x.wypos ) );

    R := rec( levels   := [],
              classes  := [],
              vertices := [],
              edges    := [] );

    for i in [1..Length(L)] do
        level := [ L[i].dim, L[i].size ];
        AddSet( R.levels, level );
        AddSet( R.classes, [ L[i].class, level ] );
        Add( R.vertices, [ L[i].wypos, level, L[i].class ] );
        for j in [1..i-1] do
            if m[j][i]<>0 then Add( R.edges, [ i, j, m[j][i] ] ); fi;
        od;
    od;

    return R;

end;
#####################################################################

IT := arg[1];
if Length(arg) = 1 then
    SG := SpaceGroupIT(3,IT);
else
    SG := SpaceGroupIT(arg[2],IT);  #arg[2] must give the dimension (either 2 or 3);
fi;

Rec := WyckoffGraphRecord(WyckoffPositions(SG));


IWPs := Difference([1..Length(Rec.vertices)],Set(List(Rec.edges,x->x[1]))); #IWPs stores the index of the IWPs;

for x in IWPs do
    Print(Rec.vertices[x][1],": ", ["point","line","plane","volume"][1 + Rec.vertices[x][2][1]],"\n");
od;

Print("Number of IWPs for Group No.", IT, ": ",Length(IWPs),"\n");

return Length(IWPs);
end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################

SGC_ResolutionForIT:=function(IT)
local T1,T2,T3,PGGen,G,Gp;

if not IsInt(IT) or IT < 1 or IT > 230 then
    Error("SGC_ResolutionForIT: IT must be 1..230\n");
fi;
if not IsBound(PGGens230) then
    Error("SGC_ResolutionForIT: PGGens230 is not loaded\n");
fi;

T1:=[[1,0,0,1],[0,1,0,0],[0,0,1,0],[0,0,0,1]];
T2:=[[1,0,0,0],[0,1,0,1],[0,0,1,0],[0,0,0,1]];
T3:=[[1,0,0,0],[0,1,0,0],[0,0,1,1],[0,0,0,1]];
PGGen:=PGGens230[IT];

G:=Group(Concatenation([T1,T2,T3],PGGen));
Gp:=AffineCrystGroupOnRight(GeneratorsOfGroup(TransposedMatrixGroup(G)));

return SGC_ResolutionSpaceGroup(Gp,SGC_MaxRelationDegreeForIT(IT)+1);
end;

#####################################################################
#####################################################################

SGC_CohomologyData:=function(arg)
local
    IT,
    C2, C2p, M, P, C3,
    PGGen, PGGen33, PGMat33, PGMat, PGMatinv, PGind, PGMat33Dict, PGindDict,
    o0,o1,o2,o3,o4,o5,
    R,CB,
    GenName_standard, Rdim,lst1,lst2,lst3,lst4,lst5,lst6,lst1to4,
    MatToPow,GapToPow,GapToPowCache,Fbarhomotopyindinv,Invofg,Prodg1g2Pow,IndToElem,
    Homotopydeg1,Homotopydeg2,Homotopydeg3,Homotopydeg4,
    func,funcs,receive,FuncVal,TopoInvdeg3,
    Gen1, Gen2, Gen3, Gen4, Gen6, GensGAP, requiredGeneratorDegree, GensDim1to4, GensDeg1to4,
    Gen3Failed, Decomp3, dimH1, dimH2, known, rk0, gcand, savedBreakOnError,
    RingData, BasesLett, Base3Lett,
    i,j,k,p,x,y;

#####################################################################
IndToElem:=function(lst,Lst)         #given lst of 0&1's, extract the elements in Lst corresponding to 1's
local i,l,out;

if (Length(lst) = Length(Lst)) = false then
    Error("IndToElem: index vector and list have different lengths (",
          Length(lst), " vs ", Length(Lst), ")\n");
fi;
l:=Length(lst);
out := [];
for i in [1..l] do
    if lst[i] = 1 then
        Append(out,[Lst[i]]);
    fi;
od;
return out;
end;
#####################################################################
MatToPow:=function(mat)            #given 4x4 matrix, output the list of powers of group generators
local i, mat33, trans;

mat33:=List([1..3],i->List([1..3],j->mat[i,j]));

i:=LookupDictionary(PGMat33Dict,mat33);
if i = fail then
    Error("MatToPow: the 3x3 rotation part is not among the enumerated point-group matrices.\n");
fi;

trans:=mat*PGMatinv[i];

return Concatenation(List([1..3],x->trans[x,4]),PGind[i]);
end;
#####################################################################
GapToPow:=function(i)            #given the index i s.t. mat:=R!.elts[i], output the list of powers of group generators
if not IsBound(GapToPowCache[i]) then
    GapToPowCache[i]:=MatToPow(R!.elts[i]);
fi;
return GapToPowCache[i];
end;
#####################################################################
Invofg:=function(v)
local vpg, transmat, p;

transmat := [[1,0,0,v[1]],[0,1,0,v[2]],[0,0,1,v[3]],[0,0,0,1]];
vpg := List([4..(Length(v))],x->v[x]);
p := LookupDictionary(PGindDict,vpg);
if p = fail then
    Error("Invofg: point-group part ", vpg, " not found in PGind.\n");
fi;

return MatToPow((transmat * PGMat[p])^(-1));
end;
#####################################################################
Prodg1g2Pow:=function(v1,v2)
local vpg1, vpg2, transmat1, transmat2, prod, p1, p2;

transmat1 := [[1,0,0,v1[1]],[0,1,0,v1[2]],[0,0,1,v1[3]],[0,0,0,1]];
transmat2 := [[1,0,0,v2[1]],[0,1,0,v2[2]],[0,0,1,v2[3]],[0,0,0,1]];
vpg1 := List([4..(Length(v1))],x->v1[x]);
vpg2 := List([4..(Length(v2))],x->v2[x]);
p1 := LookupDictionary(PGindDict,vpg1);
p2 := LookupDictionary(PGindDict,vpg2);
if p1 = fail or p2 = fail then
    Error("Prodg1g2Pow: point-group part not found in PGind (", vpg1, " or ", vpg2, ").\n");
fi;

prod := transmat1 * PGMat[p1] * transmat2 * PGMat[p2];

return MatToPow(prod);
end;
#####################################################################
Fbarhomotopyindinv:=function(i,lst)            #This is the function Fbarhomotopyindinv in Mathematica
#ResolutionSpaceGroup words are LEFT-module (letter [j,g] = Elts[g].e_j), so
#the bar-resolution contracting homotopy PREPENDS the translation:
#h(g.[t1|..|tk]) = [g|t1|..|tk]. (The old ResolutionAlmostCrystalGroup code
#appended, matching its right-module convention.)
return List(lst,x->Concatenation([i],x));
end;
#####################################################################
FuncVal:=function(lett,v)                #Given a degree-3 monomial and argument (either g1,g1,g1 or g1,g2,g2 or g1,g2,g3), evaluate the cocycle.
local deg,i,j,jval,k,val,lett1;
deg := GensDeg1to4*lett;

#Only degree-3 monomials reach this function (TopoInvdeg3 always passes
#Base3Lett letters). The old degree-1/2 branches were dead and the degree-2
#one was wrong (it indexed funcs[2] with a degree-1 position); anyone wiring
#up lower-degree invariant evaluation later should get a clear error here.
if deg <> 3 then
    Error("FuncVal: only degree-3 monomials are supported, got degree ", deg, "\n");
fi;

lett1 := ShallowCopy(lett);

for i in [1..Length(lett)] do            #finding the first generator that exists in lett
    if lett1[i] > 0 then
        break;
    fi;
od;

if i > GensDim1to4[1]+GensDim1to4[2] then                              #if a degree-3 generator
    val := funcs[3][i-GensDim1to4[1]-GensDim1to4[2]](v[1],v[2],v[3]);
else                                                                   #if not a degree-3 gen., then must be a cup prod.
    lett1[i] := lett1[i] - 1;
    for j in [1..Length(lett)] do
        if lett1[j] > 0 then
            jval := lett1[j];        #label of generator stored in j; power of this generator stored in jval
            break;
        fi;
    od;
    if j > GensDim1to4[1] then                                         #if a degree-1 gen. cup a degree-2 gen.
        val := funcs[1][i](v[1]) * funcs[2][j-GensDim1to4[1]](v[2],v[3]);
    else                                                               #if not, then must be cup of three degree-1 gens.
        lett1[j] := lett1[j] - 1;
        k := Position(lett1,1);
        val := funcs[1][i](v[1]) * funcs[1][j](v[2]) * funcs[1][k](v[3]);
    fi;
fi;
return val;
end;
#####################################################################
TopoInvdeg3:=function(arg) #usage: TopoInvdeg3(list_of_group_elements,list_of_letters,[matrices giving linear combination of letters])
local gs,letters,solrels,vallist;

gs := arg[1];               #List of group elements [g1] or [g1,g2] or [g1,g2,g3] at which the cocycles are evaluated
letters := arg[2];          #List of letters representing the monomials
vallist := fail;            #stays fail unless a topological-invariant formula matches below


if Length(arg) = 2 then
    solrels := fail;
else
    solrels := arg[3];
fi;

if Length(gs) = 1 then
    if Prodg1g2Pow(gs[1],gs[1]) = gs[1]*0 then
        vallist := List(letters,x->FuncVal(x,[gs[1],gs[1],gs[1]]));                                    #Topo inv varphi1
    else
        Print("varphi(g,g,g) is not a topological invariant!!!!\n");
    fi;
elif Length(gs) = 2 then
    if (Prodg1g2Pow(gs[2],gs[2]) = gs[2]*0) and (Prodg1g2Pow(gs[1],gs[2]) = Prodg1g2Pow(gs[2],gs[1])) then #Topo inv varphi2
        vallist := List(letters,x->FuncVal(x,[gs[1],gs[2],gs[2]])+FuncVal(x,[gs[2],gs[1],gs[2]])+FuncVal(x,[gs[2],gs[2],gs[1]]));
    else
        Print("varphi(g1,g2,g2) is not a topological invariant!!!!\n");
        Print(gs[1],gs[2],"\n");
    fi;

elif Length(gs) = 3 then

    if (Prodg1g2Pow(gs[1],gs[2]) = Prodg1g2Pow(gs[2],gs[1])) and (Prodg1g2Pow(gs[1],gs[3]) = Prodg1g2Pow(gs[3],gs[1])) and (Prodg1g2Pow(gs[2],gs[3]) = Prodg1g2Pow(gs[3],gs[2])) then     #Topo inv varphi3
        
        vallist := List(letters,x->FuncVal(x,[gs[1],gs[2],gs[3]])+FuncVal(x,[gs[1],gs[3],gs[2]])+FuncVal(x,[gs[2],gs[1],gs[3]])+FuncVal(x,[gs[2],gs[3],gs[1]])+FuncVal(x,[gs[3],gs[1],gs[2]])+FuncVal(x,[gs[3],gs[2],gs[1]]));
    
    elif ((Prodg1g2Pow(gs[2],gs[1]) = Prodg1g2Pow(Invofg(gs[1]),gs[2])) and (Prodg1g2Pow(gs[1],gs[3]) = Prodg1g2Pow(gs[3],gs[1])) and (Prodg1g2Pow(gs[2],gs[3]) = Prodg1g2Pow(gs[3],gs[2]))) then     #Topo inv tildevarphi for No. 7,26,36,39,46,57,62
        
        vallist := List(letters,x->FuncVal(x,[gs[3],Prodg1g2Pow(gs[1],gs[2]),Prodg1g2Pow(gs[1],Invofg(gs[2]))])+FuncVal(x,[gs[3],gs[1],gs[2]])+FuncVal(x,[gs[3],gs[1],Invofg(gs[2])])+FuncVal(x,[gs[3],gs[2],Invofg(gs[2])])+FuncVal(x,[Prodg1g2Pow(gs[1],gs[2]),gs[3],Prodg1g2Pow(gs[1],Invofg(gs[2]))])+FuncVal(x,[gs[1],gs[3],gs[2]])+FuncVal(x,[gs[1],gs[3],Invofg(gs[2])])+FuncVal(x,[gs[2],gs[3],Invofg(gs[2])])+FuncVal(x,[Prodg1g2Pow(gs[1],gs[2]),Prodg1g2Pow(gs[1],Invofg(gs[2])),gs[3]])+FuncVal(x,[gs[1],gs[2],gs[3]])+FuncVal(x,[gs[1],Invofg(gs[2]),gs[3]])+FuncVal(x,[gs[2],Invofg(gs[2]),gs[3]]));

    elif ((Prodg1g2Pow(gs[3],gs[2]) = Prodg1g2Pow(Invofg(gs[1]),gs[3])) and (Prodg1g2Pow(gs[3],gs[1]) = Prodg1g2Pow(gs[2],gs[3])) and (Prodg1g2Pow(gs[1],gs[2]) = Prodg1g2Pow(gs[2],gs[1]))) then     #Topo inv hatvarphi for No. 76 & 78
        
        vallist := List(letters,x->FuncVal(x,[gs[1],gs[2],Prodg1g2Pow(Invofg(gs[1]),gs[3])])+FuncVal(x,[gs[2],gs[1],Prodg1g2Pow(Invofg(gs[1]),gs[3])])+FuncVal(x,[gs[1],Prodg1g2Pow(Invofg(gs[1]),gs[3]),gs[1]])+FuncVal(x,[gs[2],gs[3],gs[2]])+FuncVal(x,[gs[3],gs[1],gs[2]])+FuncVal(x,[gs[3],gs[2],gs[1]]));

    elif ((Prodg1g2Pow(gs[3],gs[1]) = Prodg1g2Pow(Invofg(gs[1]),gs[3])) and (Prodg1g2Pow(gs[3],gs[2]) = Prodg1g2Pow(Invofg(gs[2]),gs[3])) and (Prodg1g2Pow(gs[1],gs[2]) = Prodg1g2Pow(gs[2],gs[1]))) then     #Topo inv hatvarphi for No. 4
        
        vallist := List(letters,x->FuncVal(x,[gs[1],gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[3])])+FuncVal(x,[gs[2],gs[1],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[3])])+FuncVal(x,[gs[1],Prodg1g2Pow(Invofg(gs[1]),gs[3]),gs[2]])+FuncVal(x,[gs[2],Prodg1g2Pow(Invofg(gs[2]),gs[3]),gs[1]])+FuncVal(x,[gs[3],gs[1],gs[2]])+FuncVal(x,[gs[3],gs[2],gs[1]]));

    elif ((Prodg1g2Pow(gs[3],gs[1]) = Prodg1g2Pow(Invofg(gs[2]),gs[3])) and (Prodg1g2Pow(gs[3],gs[2]) = Prodg1g2Pow(Invofg(gs[1]),gs[3])) and (Prodg1g2Pow(gs[1],gs[2]) = Prodg1g2Pow(gs[2],gs[1]))) then     #Topo inv hatvarphi for No. 9,161
        
        vallist := List(letters,x->FuncVal(x,[gs[1],gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[3])])+FuncVal(x,[gs[2],gs[1],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[3])])+FuncVal(x,[gs[1],Prodg1g2Pow(Invofg(gs[1]),gs[3]),gs[1]])+FuncVal(x,[gs[2],Prodg1g2Pow(Invofg(gs[2]),gs[3]),gs[2]])+FuncVal(x,[gs[3],gs[1],gs[2]])+FuncVal(x,[gs[3],gs[2],gs[1]]));
    else
        Print("varphi(g1,g2,g3) is not a topological invariant!!!!\n");
    fi;

elif Length(gs) = 4 then
    if (Prodg1g2Pow(gs[1],gs[2]) = Prodg1g2Pow(gs[2],gs[1])) and (Prodg1g2Pow(gs[3],gs[1]) = Prodg1g2Pow(Invofg(gs[1]),gs[3])) and (Prodg1g2Pow(gs[2],gs[3]) = Prodg1g2Pow(gs[3],gs[2])) and (Prodg1g2Pow(gs[1],gs[4]) = Prodg1g2Pow(gs[4],gs[1])) and (Prodg1g2Pow(gs[4],gs[2]) = Prodg1g2Pow(Invofg(gs[2]),gs[4])) and (Prodg1g2Pow(gs[3],gs[3]) = gs[2]) and (Prodg1g2Pow(gs[4],gs[4]) = gs[1]) then     #Topo inv varphi3 for No. 19 and 198
        
        vallist := List(letters,x->FuncVal(x,[gs[2],gs[1],Invofg(gs[1])]) + FuncVal(x,[gs[1],gs[2],Invofg(gs[1])]) + FuncVal(x,[gs[1],Invofg(gs[1]),gs[2]]) + FuncVal(x,[gs[1],gs[2],Invofg(gs[2])])+ FuncVal(x,[gs[2],gs[1],Invofg(gs[2])]) + FuncVal(x,[gs[2],Invofg(gs[2]),gs[1]]) + FuncVal(x,[gs[1],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[3]),Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[3])]) + FuncVal(x,[Prodg1g2Pow(Invofg(gs[2]),gs[3]),gs[1],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[3])]) + FuncVal(x,[Prodg1g2Pow(Invofg(gs[2]),gs[3]),Prodg1g2Pow(Invofg(gs[2]),gs[3]),gs[1]]) + FuncVal(x,[gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4]),Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])]) + FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])]) + FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2]]));
    
    elif (Prodg1g2Pow(gs[1],gs[2]) = Prodg1g2Pow(gs[2],gs[1])) and (Prodg1g2Pow(gs[3],gs[1]) = Prodg1g2Pow(Invofg(gs[1]),gs[3])) and (Prodg1g2Pow(gs[3],gs[2]) = Prodg1g2Pow(Invofg(gs[2]),gs[3])) and (Prodg1g2Pow(gs[1],gs[4]) = Prodg1g2Pow(gs[4],gs[1])) and (Prodg1g2Pow(gs[4],gs[2]) = Prodg1g2Pow(Invofg(gs[2]),gs[4])) and (Prodg1g2Pow(gs[4],gs[4]) = gs[1]) and (Prodg1g2Pow(gs[4],gs[3]) = Prodg1g2Pow(Prodg1g2Pow(gs[1],gs[2]),Prodg1g2Pow(gs[3],gs[4]))) then     #Topo inv varphi3 for No. 29

        vallist := List(letters,x->FuncVal(x,[gs[2],gs[1],Invofg(gs[1])])+FuncVal(x,[gs[1],gs[2],Invofg(gs[1])])+FuncVal(x,[gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4]),Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])])+FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])])+FuncVal(x,[gs[1],Invofg(gs[1]),gs[2]])+FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2]])+FuncVal(x,[gs[2],gs[1],Prodg1g2Pow(gs[3],gs[4])])+FuncVal(x,[gs[1],gs[2],Prodg1g2Pow(gs[3],gs[4])])+FuncVal(x,[gs[1],Prodg1g2Pow(gs[3],gs[4]),gs[2]])+FuncVal(x,[gs[3],Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2]])+FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[3],gs[2]])+FuncVal(x,[gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4]),Prodg1g2Pow(Invofg(gs[2]),gs[3])])+FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2],Prodg1g2Pow(Invofg(gs[2]),gs[3])])+FuncVal(x,[gs[2],Prodg1g2Pow(Invofg(gs[2]),gs[3]),Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])])+FuncVal(x,[gs[3],gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])]));
    
    elif (Prodg1g2Pow(gs[1],gs[2]) = Prodg1g2Pow(gs[2],gs[1])) and (Prodg1g2Pow(gs[3],gs[1]) = Prodg1g2Pow(Invofg(gs[1]),gs[3])) and (Prodg1g2Pow(gs[3],gs[2]) = Prodg1g2Pow(Invofg(gs[2]),gs[3])) and (Prodg1g2Pow(gs[1],gs[4]) = Prodg1g2Pow(gs[4],gs[1])) and (Prodg1g2Pow(gs[4],gs[2]) = Prodg1g2Pow(Invofg(gs[2]),gs[4])) and (Prodg1g2Pow(gs[4],gs[4]) = gs[1]) and (Prodg1g2Pow(gs[4],gs[3]) = Prodg1g2Pow(gs[1],Prodg1g2Pow(gs[3],gs[4]))) then     #Topo inv varphi3 for No. 33
        
        vallist :=
        List(letters,x->FuncVal(x,[gs[2],gs[1],Invofg(gs[1])])+FuncVal(x,[gs[1],gs[2],Invofg(gs[1])])+FuncVal(x,[gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4]),Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])])+FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])])+FuncVal(x,[gs[1],Invofg(gs[1]),gs[2]])+FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2]])+FuncVal(x,[gs[2],gs[1],Prodg1g2Pow(gs[3],gs[4])])+FuncVal(x,[gs[1],gs[2],Prodg1g2Pow(gs[3],gs[4])])+FuncVal(x,[gs[2],Prodg1g2Pow(gs[3],gs[4]),gs[2]])+FuncVal(x,[gs[1],Prodg1g2Pow(gs[3],gs[4]),gs[2]])+FuncVal(x,[gs[3],Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2]])+FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[3],gs[2]])+FuncVal(x,[gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4]),Prodg1g2Pow(Invofg(gs[2]),gs[3])])+FuncVal(x,[Prodg1g2Pow(Invofg(gs[1]),gs[4]),gs[2],Prodg1g2Pow(Invofg(gs[2]),gs[3])])+FuncVal(x,[gs[2],Prodg1g2Pow(Invofg(gs[2]),gs[3]),Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])])+FuncVal(x,[gs[3],gs[2],Prodg1g2Pow(Prodg1g2Pow(Invofg(gs[1]),Invofg(gs[2])),gs[4])]));
    fi;
else
    Print("Wrong in checking topological invariant: Number of group elements is not between 1 and 4!!\n");
fi;
if vallist = fail then      #no formula matched (or input flagged "not a topological invariant") -- fail loudly, not with an opaque unbound-variable error
    Error("TopoInvdeg3: no topological-invariant formula matched the given group element(s): ", gs, "\n");
fi;
if solrels = fail then
    return GF2ToZ(vallist*Z(2));
fi;
return GF2ToZ((solrels*vallist)*Z(2));
end;
#####################################################################


    ####################BEGIN TO READ THE INPUT##################

if (Length(arg) <> 1 and Length(arg) <> 2)
   or not IsInt(arg[1]) or arg[1] < 1 or arg[1] > 230 then
    Error("cohomology API: expected IT or (IT,R), with IT in 1..230\n");
fi;

IT:=arg[1];

#The data globals below come from gap/data.gi. Fail with a clear message if
#that file was not Read in first.
if not (IsBound(PGGens230) and IsBound(funcs230) and IsBound(IWP) and IsBound(GENNAMES)) then
    Error("cohomology API: required data not loaded -- Read the data file defining PGGens230, funcs230, IWP and GENNAMES before calling.\n");
fi;
if Length(arg) = 2 then
    if not (IsRecord(arg[2]) or IsComponentObjectRep(arg[2])) or not IsBound(arg[2]!.dimension)
       or not IsBound(arg[2]!.boundary) or not IsBound(arg[2]!.elts) then
        Error("cohomology API: R must be a HAP resolution\n");
    fi;
    if Length(Size(arg[2]))-1 < SGC_MaxRelationDegreeForIT(IT) then
        Error("cohomology API: R does not reach the required relation degree ",
              SGC_MaxRelationDegreeForIT(IT), " for IT=", IT, "\n");
    fi;
fi;

PGMat33 := [];
PGMat := [];     #forward 4x4 point-group representatives; PGMat[i] = PGMatinv[i]^(-1),
                 #precomputed so the hot paths (Invofg, Prodg1g2Pow, LSM scan) never invert
PGMatinv := [];
PGind := [];

PGGen:=PGGens230[IT];
funcs:=funcs230[IT];

PGGen33 := List([1..Length(PGGen)],k->List([1..3],i->List([1..3],j->PGGen[k][i,j])));

if Length(PGGen) = 0 then
    Append(PGind,[[]]);
    PGMat33:=[[[1,0,0],[0,1,0],[0,0,1]]];
    PGMat:=[[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]];
    PGMatinv:=[[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]];
elif Length(PGGen) = 1 then
    for o1 in [0..(Order(PGGen33[1])-1)] do
        Append(PGind,[[o1]]);
        Append(PGMat33,[PGGen33[1]^o1]);
        Append(PGMat,[PGGen[1]^o1]);
        Append(PGMatinv,[PGMat[Length(PGMat)]^(-1)]);
    od;
elif Length(PGGen) = 2 then
    for o1 in [0..(Order(PGGen33[1])-1)] do
        for o2 in [0..1] do  #normal form: IT 75-88 have Order(PGGen[2])=4 but PGGen[2]^2=PGGen[1]
            Append(PGind,[[o1,o2]]);
            Append(PGMat33,[PGGen33[1]^o1*PGGen33[2]^o2]);
            Append(PGMat,[PGGen[1]^o1*PGGen[2]^o2]);
            Append(PGMatinv,[PGMat[Length(PGMat)]^(-1)]);
        od;
    od;
elif Length(PGGen) = 3 then
    for o1 in [0..(Order(PGGen33[1])-1)] do
        for o2 in [0..1] do  #normal form: IT 75-88 have Order(PGGen[2])=4 but PGGen[2]^2=PGGen[1]
            for o3 in [0..(Order(PGGen33[3])-1)] do
                Append(PGind,[[o1,o2,o3]]);
                Append(PGMat33,[PGGen33[1]^o1*PGGen33[2]^o2*PGGen33[3]^o3]);
                Append(PGMat,[PGGen[1]^o1*PGGen[2]^o2*PGGen[3]^o3]);
                Append(PGMatinv,[PGMat[Length(PGMat)]^(-1)]);
            od;
        od;
    od;
elif Length(PGGen) = 4 then
    for o1 in [0..(Order(PGGen33[1])-1)] do
        for o2 in [0..1] do  #normal form: IT 75-88 have Order(PGGen[2])=4 but PGGen[2]^2=PGGen[1]
            for o3 in [0..(Order(PGGen33[3])-1)] do
                for o4 in [0..(Order(PGGen33[4])-1)] do
                    Append(PGind,[[o1,o2,o3,o4]]);
                    Append(PGMat33,[PGGen33[1]^o1*PGGen33[2]^o2*PGGen33[3]^o3*PGGen33[4]^o4]);
                    Append(PGMat,[PGGen[1]^o1*PGGen[2]^o2*PGGen[3]^o3*PGGen[4]^o4]);
                    Append(PGMatinv,[PGMat[Length(PGMat)]^(-1)]);
                od;
            od;
        od;
    od;
elif Length(PGGen) = 5 then
    for o1 in [0..(Order(PGGen33[1])-1)] do
        for o2 in [0..1] do  #normal form: IT 75-88 have Order(PGGen[2])=4 but PGGen[2]^2=PGGen[1]
            for o3 in [0..(Order(PGGen33[3])-1)] do
                for o4 in [0..(Order(PGGen33[4])-1)] do
                    for o5 in [0..(Order(PGGen33[5])-1)] do
                        Append(PGind,[[o1,o2,o3,o4,o5]]);
                        Append(PGMat33,[PGGen33[1]^o1*PGGen33[2]^o2*PGGen33[3]^o3*PGGen33[4]^o4*PGGen33[5]^o5]);
                        Append(PGMat,[PGGen[1]^o1*PGGen[2]^o2*PGGen[3]^o3*PGGen[4]^o4*PGGen[5]^o5]);
                        Append(PGMatinv,[PGMat[Length(PGMat)]^(-1)]);
                    od;
                od;
            od;
        od;
    od;
else
    Print("Number of Point Group Generators Exceeds 5 -- WRONG!!!");
fi;

#These dictionaries deliberately retain the first occurrence, matching Position.
#They are private to this Context's closures and therefore tied to its exact
#point-group enumeration and resolution.
PGindDict:=NewDictionary(PGind[1],true);
PGMat33Dict:=NewDictionary(PGMat33[1],true);
for i in [1..Length(PGind)] do
    if LookupDictionary(PGindDict,PGind[i]) = fail then
        AddDictionary(PGindDict,PGind[i],i);
    fi;
    if LookupDictionary(PGMat33Dict,PGMat33[i]) = fail then
        AddDictionary(PGMat33Dict,PGMat33[i],i);
    fi;
od;


if Length(arg) = 1 then
    R:=SGC_ResolutionForIT(IT);
else
    R:=arg[2];
fi;
GapToPowCache:=[];
#Resolution for the group now available.
#R!.elts are column-convention affine matrices in the group generated by
#[T1,T2,T3] and PGGen (CrystallographicComplex standardizes -- a no-op here,
#TranslationBasis is the identity -- and transposes back before returning),
#so GapToPow feeds them to MatToPow directly. MatToPow errors loudly if an
#element ever falls outside the enumerated point-group matrices.


Homotopydeg1:=List([1..R!.dimension(1)],x->List(R!.boundary(1,x),y->[y[2]]));
Homotopydeg2:=List([1..R!.dimension(2)],x->Concatenation(List(R!.boundary(2,x),y->Fbarhomotopyindinv(y[2],Homotopydeg1[AbsInt(y[1])]))));
Homotopydeg3:=List([1..R!.dimension(3)],x->Concatenation(List(R!.boundary(3,x),y->Fbarhomotopyindinv(y[2],Homotopydeg2[AbsInt(y[1])]))));


CB:=[];
for p in [1..4] do
CB[p]:=CR_Mod2CocyclesAndCoboundaries(R,p,true);
od;


Gen1:=[];
for func in funcs[1] do
    Append(Gen1,[CB[1].cocycleToClass(List([1..R!.dimension(1)],x->RemInt(Sum(List(Homotopydeg1[x],y->func(GapToPow(y[1])))),2)))]);
od;


Gen2:=[];
for func in funcs[2] do
    Append(Gen2,[CB[2].cocycleToClass(List([1..R!.dimension(2)],x->RemInt(Sum(List(Homotopydeg2[x],y->func(GapToPow(y[1]),GapToPow(y[2])))),2)))]);
od;

#Explicit degree-3 cocycle functions are evaluated on the resolution via the
#contracting homotopy. For a few groups (the R-suffixed generator names of
#IT 225/227/229) the stored function does not transport to a cocycle on this
#resolution -- historically those slots carried vectors hardcoded in the old
#ResolutionAlmostCrystalGroup basis. Such slots are now filled deterministically
#with GAP-computed ring generators, chosen independent of the degree-3
#decomposable classes so they are genuine generators in the current basis.
#The LSM section is unaffected: it evaluates the explicit functions directly at
#topological-invariant configurations and has its own full-rank consistency check.
Gen3:=[];
Gen3Failed:=[];
#BreakOnError must be off while probing, or a failed transport would drop a
#batch-mode session into the break loop instead of unwinding to the catch.
savedBreakOnError := BreakOnError;
for func in funcs[3] do
    BreakOnError := false;
    receive := CALL_WITH_CATCH(function()
        return CB[3].cocycleToClass(List([1..R!.dimension(3)],x->RemInt(Sum(List(Homotopydeg3[x],y->func(GapToPow(y[1]),GapToPow(y[2]),GapToPow(y[3])))),2)));
    end, []);
    BreakOnError := savedBreakOnError;
    if receive[1] = true and receive[2] <> fail then
        Append(Gen3,[receive[2]]);
    else
        Append(Gen3,[fail]);
        Append(Gen3Failed,[Length(Gen3)]);
    fi;
od;

#Compute generic generators once, and only to the highest degree actually
#needed. Explicit Gen1/Gen2 and successfully transported Gen3 remain authoritative.
requiredGeneratorDegree:=0;
if (IT in SGC_Degree6GeneratorGroups) = true then
    requiredGeneratorDegree:=6;
elif (IT in SGC_Degree4GeneratorGroups) = true or Length(Gen3Failed) > 0 then
    requiredGeneratorDegree:=4;
fi;
if requiredGeneratorDegree > 0 then
    GensGAP:=Mod2RingGenerators(R,requiredGeneratorDegree,3);
fi;

if Length(Gen3Failed) > 0 then
    #Degree-3 decomposable classes: cup products H^1 x H^2 (this includes all
    #triple products of degree-1 classes, since H^1 x H^1 lands in H^2).
    dimH1 := Length(CB[1].cocycleToClass(List([1..R!.dimension(1)],x->0)));
    dimH2 := Length(CB[2].cocycleToClass(List([1..R!.dimension(2)],x->0)));
    Decomp3 := [];
    for x in CohomologyBasis(List([1..dimH1],i->1)) do
        for y in CohomologyBasis(List([1..dimH2],i->1)) do
            Append(Decomp3,[Mod2CupProduct(R,x,y,1,2,CB[1],CB[2],CB[3])]);
        od;
    od;
    #Steinitz-style completion: the known span is the decomposables plus the
    #successfully transported explicit classes; each failed slot takes the first
    #GAP generator that enlarges that span (Mod2RingGenerators guarantees its
    #degree-3 output spans H^3 modulo the decomposables, so enough exist).
    known := Concatenation(Decomp3, Filtered(Gen3, x -> x <> fail));
    for k in Gen3Failed do
        rk0 := SGC_RankMod2(known);
        for gcand in GensGAP[3] do
            if Gen3[k] = fail and SGC_RankMod2(Concatenation(known,[gcand])) > rk0 then
                Gen3[k] := gcand;
                Append(known,[gcand]);
            fi;
        od;
        if Gen3[k] = fail then
            Error("SpaceGroupCohomologyRingGapInterface: cannot complete degree-3 generator slot ", k, " for IT=", IT);
        fi;
    od;
fi;

#The explicit degree-3 cocycle functions must give linearly independent classes;
#a dependent set means funcs230/data.gi is inconsistent with this resolution.
if Length(Gen3) > 0 and Length(SGC_RowBasisMod2(Gen3)) < Length(Gen3) then
    Error("SpaceGroupCohomologyRingGapInterface: dependent degree-3 explicit generators for IT=", IT);
fi;


Gen4:=[];
if (IT in SGC_Degree4GeneratorGroups) = true then
    Gen4 := GensGAP[4];
fi;
Gen6:=[];
if (IT in SGC_Degree6GeneratorGroups) = true then
    Gen6 := GensGAP[6];
fi;

#GAP-computed degree-4 generators are a Steinitz completion and must be independent;
#a dependent set means Mod2RingGenerators returned an inconsistent basis.
if Length(Gen4) > 0 and Length(SGC_RowBasisMod2(Gen4)) < Length(Gen4) then
    Error("SpaceGroupCohomologyRingGapInterface: dependent degree-4 generators for IT=", IT);
fi;


RingData := Mod2RingGensAndRels(IT,3,R,[Gen1,Gen2,Gen3,Gen4,[],Gen6],true,
                                SGC_MaxRelationDegreeForIT(IT));
BasesLett := RingData.bases;

GensDim1to4 := [Length(Gen1),Length(Gen2),Length(Gen3),Length(Gen4)];
GensDeg1to4:=Concatenation(List([1..Length(Gen1)],x->1),List([1..Length(Gen2)],x->2),List([1..Length(Gen3)],x->3),List([1..Length(Gen4)],x->4));

Base3Lett := BasesLett[3];

return rec(
    isSGCCohomologyContext := true,
    IT := IT,
    resolution := R,
    ring := RingData,
    gen3Failed := Gen3Failed,
    base3Letters := Base3Lett,
    generatorDimensions := GensDim1to4,
    generatorDegrees := GensDeg1to4,
    pointGroupIndices := PGind,
    pointGroupMatrices := PGMat,
    topoInvdeg3 := TopoInvdeg3,
    productPower := Prodg1g2Pow,
    indToElem := IndToElem
);
end;

#####################################################################
#####################################################################

SGC_WPCohomologyData:=function(Context)
local
    IT, Base3Lett, GensDim1to4, PGind, PGMat,
    TopoInvdeg3, Prodg1g2Pow, IndToElem,
    LSMLett, Mat, mat2, vec, LSMMat, CountLSM,
    g1, g2, g2Square, Candidates, candidate1, candidate2,
    pointIndex, x1, y1, z1, identityAffine, Span, addCandidate,
    scanCandidates, rank,
    table, coordinates, entry, label, key, target, equalPosition,
    zeroCoordinate, j, x;

if not (IsRecord(Context) or IsComponentObjectRep(Context))
   or not IsBound(Context.isSGCCohomologyContext)
   or Context.isSGCCohomologyContext <> true then
    Error("SGC_WPCohomologyData: expected a prepared cohomology Context\n");
fi;

IT:=Context.IT;
Base3Lett:=Context.base3Letters;
GensDim1to4:=Context.generatorDimensions;
PGind:=Context.pointGroupIndices;
PGMat:=Context.pointGroupMatrices;
TopoInvdeg3:=Context.topoInvdeg3;
Prodg1g2Pow:=Context.productPower;
IndToElem:=Context.indToElem;

#Print all the elements of the mod-2 cohomology at degree 3:


############################### BELOW ARE LSM RELATED CODES ###############################


Mat:=[];


#First: record all the LSM TIs listed in the IWP table from gap/data.gi.
#
#
CountLSM := [];
for x in IWP[IT] do
    if (x[2] = []) = false then
        Append(Mat,[TopoInvdeg3(x[2],Base3Lett)]);
        Append(CountLSM,[x[2]]);
    fi;
od;


Span:=SGC_NewSpanTesterMod2(Mat);
identityAffine:=[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]];

addCandidate:=function(candidateVector)
if Span.addIfIndependent(candidateVector) then
    Add(Mat,candidateVector);
fi;
return Span.rank()=Length(Base3Lett);
end;


#Second: find all the non-LSM TIs, which are of one of the following four types:
#
#
scanCandidates:=function()
if Span.rank()=Length(Base3Lett) then
    return;
fi;

#Build the bounded candidates once, in the old PGind -> x -> y -> z order.
#The affine matrix is used only for order/trace classification and the inner
#commutativity test; Context.productPower remains authoritative for formulas.
Candidates:=[];
for pointIndex in [1..Length(PGind)] do
    for x1 in [-2..2] do
        for y1 in [-2..2] do
            for z1 in [-2..2] do
                Add(Candidates,rec(
                    powerVector:=Concatenation([x1,y1,z1],PGind[pointIndex]),
                    affineMatrix:=[[1,0,0,x1],[0,1,0,y1],[0,0,1,z1],[0,0,0,1]]*PGMat[pointIndex],
                    pointIndex:=pointIndex
                ));
            od;
        od;
    od;
od;

for candidate2 in Candidates do
    if (IT <= 220) or (PGind[candidate2.pointIndex][3] = 0) then
        g2:=candidate2.powerVector;
        mat2:=candidate2.affineMatrix;
        if mat2^2=identityAffine then
            if Trace(mat2)=0 then                    #C2 rotation
                vec:=TopoInvdeg3([g2],Base3Lett);
                if addCandidate(vec) then
                    return;
                fi;
            elif Trace(mat2)=2 then                  #Mirror
                vec:=TopoInvdeg3([g2],Base3Lett);
                if addCandidate(vec) then
                    return;
                fi;
                for candidate1 in Candidates do
                    g1:=candidate1.powerVector;
                    if g1<>g1*0
                       and mat2*candidate1.affineMatrix=candidate1.affineMatrix*mat2 then
                        vec:=TopoInvdeg3([g1,g2],Base3Lett);
                        if addCandidate(vec) then
                            return;
                        fi;
                    fi;
                od;
            fi;
        elif mat2^4=identityAffine and Trace(mat2)=2 then         #C4 rotation
            g2Square:=Prodg1g2Pow(g2,g2);
            vec:=TopoInvdeg3([g2,g2Square],Base3Lett);
            if addCandidate(vec) then
                return;
            fi;
        fi;
    fi;
od;
end;

scanCandidates();
rank:=Span.rank();

if rank = Length(Base3Lett) and rank = Length(Mat) then

    LSMMat := List(TransposedMat(InverseMatMod(Mat,2)));
    LSMLett := List([1..Length(CountLSM)],x->LSMMat[x]);
else
    Error("WPCohomologyTable: full rank not achieved: ", rank,
          "!=", Length(Base3Lett), " or ", rank,
          "!=", Length(Mat), "\n");
fi;

table:=rec();
coordinates:=rec();
zeroCoordinate:=List([1..Length(Base3Lett)],x->0);
j:=1;
for entry in IWP[IT] do
    label:=entry[1];
    equalPosition:=Position(label,'=');
    if equalPosition = fail then
        key:=label;
        target:=fail;
    else
        key:=label{[1..equalPosition-1]};
        target:=label{[equalPosition+1..Length(label)]};
    fi;

    if entry[2] <> [] then
        coordinates.(key):=ShallowCopy(LSMLett[j]);
        j:=j+1;
    elif target = "empty" or target = fail then
        coordinates.(key):=ShallowCopy(zeroCoordinate);
    else
        if not IsBound(coordinates.(target)) then
            Error("WPCohomologyTable: alias target ", target,
                  " is not defined before ", key, " for IT=", IT, "\n");
        fi;
        coordinates.(key):=ShallowCopy(coordinates.(target));
    fi;
    table.(key):=IndToElem(coordinates.(key),Base3Lett);
od;

return rec(
    isSGCWPCohomologyData:=true,
    context:=Context,
    table:=table,
    coordinates:=coordinates,
    rawIWP:=IWP[IT],
    lsmLetters:=LSMLett
);
end;

#####################################################################
#####################################################################

SGC_PrintWPCohomologyData:=function(Context,WPData)
local entry,j;

Print("LSM:\n");
j:=1;
for entry in WPData.rawIWP do
    Print(entry[1]," ");
    if entry[2] <> [] then
        PrintMonomialString(Context.indToElem(WPData.lsmLetters[j],Context.base3Letters),
                            Context.generatorDimensions,"+",
                            Context.ring.generators.names,"\n");
        j:=j+1;
    else
        Print("\n");
    fi;
od;
end;

#####################################################################
#####################################################################

SGC_PrintCohomologyClass:=function(Context,Class)
if Length(Class) = 0 then
    Print("0\n");
else
    PrintMonomialString(Class,Context.ring.generatorDimensions,"+",
                        Context.ring.generators.names,"\n");
fi;
end;

#####################################################################
#####################################################################

GroupCohomologyMod2:=function(arg)
local Context, R;

if Length(arg) = 1 and (IsRecord(arg[1]) or IsComponentObjectRep(arg[1]))
   and IsBound(arg[1].isSGCCohomologyContext)
   and arg[1].isSGCCohomologyContext = true then
    Context:=arg[1];
elif Length(arg) = 1 and IsInt(arg[1]) then
    R:=SGC_ResolutionForIT(arg[1]);
    Context:=SGC_CohomologyData(arg[1],R);
elif Length(arg) = 2 and IsInt(arg[1]) then
    R:=arg[2];
    Context:=SGC_CohomologyData(arg[1],R);
else
    Error("GroupCohomologyMod2: expected Context, IT, or (IT,R)\n");
fi;

SGC_PrintMod2RingData(Context.ring);
end;

#####################################################################
#####################################################################

WPCohomologyTable:=function(arg)
local Context, R, WPData, entry, equalPosition, key;

if Length(arg) >= 1 and (IsRecord(arg[1]) or IsComponentObjectRep(arg[1]))
   and IsBound(arg[1].isSGCCohomologyContext)
   and arg[1].isSGCCohomologyContext = true then
    Context:=arg[1];
    if Length(arg) = 1 then
        WPData:=SGC_WPCohomologyData(Context);
    elif Length(arg) = 2 and (IsRecord(arg[2]) or IsComponentObjectRep(arg[2]))
         and IsBound(arg[2].isSGCWPCohomologyData)
         and arg[2].isSGCWPCohomologyData = true then
        WPData:=arg[2];
        if not IsIdenticalObj(WPData.context,Context) then
            Error("WPCohomologyTable: WPData belongs to a different Context\n");
        fi;
    else
        Error("WPCohomologyTable: expected Context or (Context,WPData)\n");
    fi;
elif Length(arg) = 1 and IsInt(arg[1]) then
    R:=SGC_ResolutionForIT(arg[1]);
    Context:=SGC_CohomologyData(arg[1],R);
    WPData:=SGC_WPCohomologyData(Context);
elif Length(arg) = 2 and IsInt(arg[1]) then
    R:=arg[2];
    Context:=SGC_CohomologyData(arg[1],R);
    WPData:=SGC_WPCohomologyData(Context);
else
    Error("WPCohomologyTable: expected Context, (Context,WPData), IT, or (IT,R)\n");
fi;

for entry in WPData.rawIWP do
    equalPosition:=Position(entry[1],'=');
    if equalPosition = fail then
        key:=entry[1];
    else
        key:=entry[1]{[1..equalPosition-1]};
    fi;
    Print(key," ");
    SGC_PrintCohomologyClass(Context,WPData.table.(key));
od;
end;

#####################################################################
#####################################################################

WPCohomologyClass:=function(arg)
local R, WPs, Context, WPData, rawLabels, coordinate, wp, equalPosition, key, i, Class;

if Length(arg) >= 2 and (IsRecord(arg[1]) or IsComponentObjectRep(arg[1]))
   and IsBound(arg[1].isSGCCohomologyContext)
   and arg[1].isSGCCohomologyContext = true then
    Context:=arg[1];
    if Length(arg) = 2 then
        WPs:=arg[2];
        WPData:=SGC_WPCohomologyData(Context);
    elif Length(arg) = 3 and (IsRecord(arg[2]) or IsComponentObjectRep(arg[2]))
         and IsBound(arg[2].isSGCWPCohomologyData)
         and arg[2].isSGCWPCohomologyData = true then
        WPData:=arg[2];
        WPs:=arg[3];
        if not IsIdenticalObj(WPData.context,Context) then
            Error("WPCohomologyClass: WPData belongs to a different Context\n");
        fi;
    else
        Error("WPCohomologyClass: expected (Context,WPs) or (Context,WPData,WPs)\n");
    fi;
elif Length(arg) = 2 and IsInt(arg[1]) then
    WPs:=arg[2];
    R:=SGC_ResolutionForIT(arg[1]);
    Context:=SGC_CohomologyData(arg[1],R);
    WPData:=SGC_WPCohomologyData(Context);
elif Length(arg) = 3 and IsInt(arg[1]) then
    R:=arg[2];
    WPs:=arg[3];
    Context:=SGC_CohomologyData(arg[1],R);
    WPData:=SGC_WPCohomologyData(Context);
else
    Error("WPCohomologyClass: expected prepared Context data or an IT signature\n");
fi;
if not IsList(WPs) or IsString(WPs)
   or not ForAll(WPs,IsString) then
    Error("WPCohomologyClass: WPs must be a list of Wyckoff-position labels\n");
fi;
rawLabels:=List(WPData.rawIWP,x->x[1]);
coordinate:=List([1..Length(Context.base3Letters)],x->0);
for wp in WPs do
    equalPosition:=Position(wp,'=');
    if equalPosition = fail then
        key:=wp;
    else
        if not wp in rawLabels then
            Error("WPCohomologyClass: unknown annotated Wyckoff position ",wp,
                  " for IT=",Context.IT,"\n");
        fi;
        key:=wp{[1..equalPosition-1]};
    fi;
    if not IsBound(WPData.coordinates.(key)) then
        Error("WPCohomologyClass: unknown Wyckoff position ",wp,
              " for IT=",Context.IT,"\n");
    fi;
    for i in [1..Length(coordinate)] do
        coordinate[i]:=RemInt(coordinate[i]+WPData.coordinates.(key)[i],2);
    od;
od;
Class:=Context.indToElem(coordinate,Context.base3Letters);
SGC_PrintCohomologyClass(Context,Class);
end;

#####################################################################
#####################################################################

SpaceGroupCohomologyRingGapInterface:=function(arg)
local R, Context, WPData;

if Length(arg) <> 1 or not IsInt(arg[1]) or arg[1] < 1 or arg[1] > 230 then
    Error("SpaceGroupCohomologyRingGapInterface: IT must be 1..230\n");
fi;

R:=SGC_ResolutionForIT(arg[1]);
Context:=SGC_CohomologyData(arg[1],R);
Print("===========================================\n");
Print("Mod-2 Cohomology Ring of Group No. ", Context.IT, ":\n");
SGC_PrintMod2RingData(Context.ring);
Print("===========================================\n");
WPData:=SGC_WPCohomologyData(Context);
SGC_PrintWPCohomologyData(Context,WPData);
return true;
end;


#####################################################################
#####################################################################
