#import Pkg; Pkg.add("AbstractAlgebra");
#import AbstractAlgebra: degree;
###??? import AbstractAlgebra: is_one;  #is_one;

# A function similar to this is in StatsBase.jl but I wish to avoid
# having an external dependency (for just 1 simple fn!)

# Generate random m-subset of {1,2,3,...,n}
# Result is a Int64[]; entries are NOT SORTED!!
function random_subset(n::Int64, m::Int64)
    # assume n >= 1, m >= 0 and m <= n
    if m == 0 #=then=#
        return Int64[];
    end #=if=#
    L = collect(1:n);
    if m == n #=then=#
        return L;
    end #=if=#
    for j in 1:m #=do=#
        k = rand(j:n);
        L[j],L[k] = L[k],L[j]; # just a SWAP
    end #=for=#
    L = first(L,m);
    #?? sort!(L); ??
    return L;
end #=function=#

# There must be a better way...!
# Split a "list" into 2 parts determined by a predicate.
# Returns 2-tuple: list-of-sat-elems, list-of-unsat-elems
function filter2(pred::Function, L::Vector)
    sat = [];
    unsat = [];
    for x in L #=do=#
        if pred(x)
            push!(sat,x);
        else
            push!(unsat,x);
        end #=if=#
    end #=for=#
    return sat,unsat;
end #=function=#

############################################
# Code for representing & manipulating PPs (power products, aka. monomials)
# All this code is "local" to this file, & not exported!

# type alias
const HSNumVar = AbstractAlgebra.Generic.LaurentMPolyWrap{QQFieldElem, QQMPolyRingElem, AbstractAlgebra.Generic.LaurentMPolyWrapRing{QQFieldElem, QQMPolyRing}}


# Each PP is represented as Vector{PP_exponent}

PP_exponent = Int64;  # UInt ???  Strange: Int32 was slower on my machine ?!?


#= mutable =# struct PP
    expv::Vector{PP_exponent};
end # struct

function Base.copy(t::PP)
    return PP(copy(t.expv));
end #=function=#

# RETURN VALUE???  Perhaps index or 0? (or -1?)
function IsSimplePowerPP(t::PP)
    CountNZ = 0;
    for i in 1:length(t.expv) #do
        @inbounds if (t.expv[i] == 0)  continue;  end #if
        if (CountNZ > 0)  return false; end #if
        CountNZ = i;
    end #for
    if (CountNZ != 0)  return true; end #if    MAYBE RETURN index & exp???
    return false; # because t == 1
end #function

# Should be is_one, but julia complained :-(
function isone(t::PP)
    return all(t.expv .== 0);
end #function

function degree(t::PP)
    return sum(t.expv);
end #function

function IsDivisible(t::PP, s::PP)  # is t divisible by s
    n = length(t.expv); # assume equal to length(s.expv);
    for i in 1:n #do
        @inbounds if t.expv[i] < s.expv[i] #then
            return false;
        end #if
    end #for
    return true;
end #function


# modifies first arg
function mult_by_var!(t::PP, j::Int64)
    @inbounds t.expv[j] += 1;
end #function

function mult(t1::PP, t2::PP)
    # ASSUMES: length(t1.expv) == length(t2.expv)
    return PP(t1.expv + t2.expv);
end #function


function divide(t1::PP, t2::PP)
    # ASSUMES: length(t1.expv) == length(t2.expv), also that t1 is mult of t2
    return PP(t1.expv - t2.expv);
end #function

function is_coprime(t1::PP, t2::PP)
    # ASSUMES: length(t1.expv) == length(t2.expv)
    n = length(t1.expv);
    for i in 1:n #do
        @inbounds if t1.expv[i] != 0 && t2.expv[i] != 0  #=then=#
            return false;
        end #=if=#
    end #for
    return true;
end #function

function lcm(t1::PP, t2::PP)
    # ASSUMES: length(t1.expv) == length(t2.expv)
    n = length(t1.expv);
    expv = [0  for _ in 1:n];
    for i in 1:n #do
        @inbounds expv[i] = max(t1.expv[i], t2.expv[i]);
    end #for
    return PP(expv);
end #function

function gcd(t1::PP, t2::PP)
    # ASSUMES: length(t1.expv) == length(t2.expv)
    n = length(t1.expv);
    expv = [0  for _ in 1:n];
    for i in 1:n #do
        @inbounds expv[i] = min(t1.expv[i], t2.expv[i]);
    end #for
    return PP(expv);
end #function

function gcd3(t1::PP, t2::PP, t3::PP)
    # ASSUMES: length(t1.expv) == length(t2.expv) == length(t3.expv)
    n = length(t1.expv);
    expv = [0  for _ in 1:n];
    for i in 1:n #do
        @inbounds expv[i] = min(t1.expv[i], t2.expv[i], t3.expv[i]);
    end #for
    return PP(expv);
end #function

function colon(t1::PP, t2::PP)
    # ASSUMES: length(t1.expv) == length(t2.expv)
    n = length(t1.expv);
    expv = [0  for _ in 1:n];
    for i in 1:n #do
        @inbounds expv[i] = max(0, t1.expv[i]-t2.expv[i]);
    end #for
    return PP(expv);
end #function

function saturatePP(t1::PP, t2::PP)
    # ASSUMES: length(t1.expv) == length(t2.expv)
    n = length(t1.expv);
    expv = [0  for _ in 1:n];
    for i in 1:n #do
        @inbounds if t2.expv[i] == 0 #then
            @inbounds expv[i] = t1.expv[i];
        end #if
    end #for
    return PP(expv);
end #function

function radical(t::PP)
    n = length(t.expv);
    expv = [0  for _ in 1:n];
    for i in 1:n #do
        @inbounds if t.expv[i] > 0 #then
            expv[i] = 1;
        end #if
    end #for
    return PP(expv);
end #function


# if t1 == t2 then returns false.
function DegRevLexLess(t1::PP, t2::PP)
    d1 = sum(t1.expv);
    d2 = sum(t2.expv);
    if d1 != d2 #=then=# return (d1 < d2); end #=if=#
    nvars = length(t1.expv);
    for i in nvars:-1:1 #=do=#
        if t1.expv[i] != t2.expv[i] #=then=#
            return (t1.expv[i] > t2.expv[i]);
        end #=if=#
    end #=for=#
    return false;
end #=function=#

# Test whether PP involves at least 1 indet from the list K (of indexes)
function involves(t::PP, K::Vector{Int64}) # K is index set
    # ASSUMES: indexes in K are all in range
    for i in K #do
        if t.expv[i] != 0 #then
            return true;
        end #if
    end #for
    return  false;
end #function


function Base.show(io::IO, t::PP)
    if (all(t.expv .== 0))   # (isone(PP))  # why doesn't isone work ???
        print(io, "1");
        return;
    end #if
    str = "";
    n = length(t.expv);
    for i in 1:n #do
        if (t.expv[i] != 0) #then
            if (str != "")  str = str * "*";  end #if
            str = str * "x[$(i)]";
            if (t.expv[i] > 1) #then
                str = str * "^$(t.expv[i])";
            end #if
        end #if
    end #for
    print(io, str);
end #function


# interreduce list of PPs; equiv find min set of gens for PP monoideal
function interreduce(L::Vector{PP})  # L is list of PPs
    sort!(L, by=degree);
    MinGens = PP[];##empty Vector{PP}();
    for t in L #do
        discard = false;
        for s in MinGens #do
            if IsDivisible(t,s) #then
                discard = true;
                break;
            end #if
        end #for
        if !discard #then
            push!(MinGens, t);
        end #if
    end # for
    return MinGens;
end # function


# Is t a multiple of at least one element of L?
# Add degree truncation???
function NotMultOf(L::Vector{PP}, t::PP)
    for s in L #=do=#
        if IsDivisible(t,s)
            return false;
        end #=if=#
    end #=for=#
    return true;
end #=function=#



# "project" PP onto sub-monoid of PPs gen by indets in indexes
function ProjectIndets(t::PP, indexes::Vector{Int})
    # # expv = [0  for _ in 1:length(indexes)];
    # # for i in 1:length(indexes) #do
    # #     expv[i] = t.expv[indexes[i]];
    # # end #for
    # # return PP(expv);
    return PP([t.expv[k]  for k in indexes]);
end #function


# NOT SURE THIS IS USEFUL: how many indets are really needed in the list L?
function TrueNumVars(L::Vector{PP})
    # assume L non empty?
    if  isempty(L)  return 0;  end #if
##    MaxNumVars = length(L[1].expv);
##    AllVars = PP([1  for _ in 1:MaxNumVars]);
    t = radical(L[1]);
    for j in 2:length(L)  #do
        t = radical(mult(t,L[j]));
    end #for
    return degree(t);        
end  #function


#-----------------------------------------------------------------------------
# This function is also private/local to this file.

# Each PP in the list L is interpreted as saying that all indets
# involved in that PP are "connected".  The task is to find (minimal)
# connected componnts of indets.

# Connected components of variables

# Input: non-empty list of PPs
# Output: list of lists of var indexes, each sublist is a connected component
function ConnectedComponents(L::Vector{PP})
    ConnCompt::Vector{Vector{Int64}} = [];
    nvars = length(L[1].expv);
    IgnoreVar = [ false   for _ in 1:nvars];
    VarAppears = copy(IgnoreVar);
    for t in L #do
        for j in 1:nvars #do
            @inbounds if t.expv[j] != 0 #then
                @inbounds VarAppears[j] = true;
            end #if
        end #for
    end #for
    CountIgnore = 0;
    for j in 1:nvars #do
        @inbounds if !VarAppears[j] #then
            @inbounds IgnoreVar[j] = true;
            CountIgnore += 1;
        end #if
    end #for
###println("CountIgnore=$(CountIgnore)");
###println("nvars=$(nvars)");
    while CountIgnore < nvars #do
        ## Maybe use findfirst instead of loop below?
        # j=1;
        # while IgnoreVar[j] #do
        #     j += 1;
        # end #while
        j = findfirst(!, IgnoreVar); # j is index of some var which appears in at least one PP
        k = findfirst((t -> t.expv[j] != 0), L); # pick some PP involving j-th var
        lcm = L[k];
###println("lcm=$(lcm)");
        DoAnotherIteration = true;
        while DoAnotherIteration #do
            DoAnotherIteration = false;
            for t in L #do
                if is_coprime(lcm,t)  continue; end
                s = saturatePP(t,lcm);
                if isone(s)  continue;  end
                lcm = mult(lcm,s); ### lcm *= s;
                DoAnotherIteration = true;
            end #for
        end #while
        vars = filter((k -> lcm.expv[k] > 0), 1:nvars);
        # remove conn compt from L???
#seems to be slower with this line ?!?        L = filter((t -> is_coprime(t,lcm)), L);
        push!(ConnCompt, vars);
        for k in vars #do
            IgnoreVar[k] = true;
            CountIgnore += 1;
        end #for
    end #while
    return ConnCompt;
end #function

#############################################################################
#-----------------------------------------------------------------------------

# Base cases for HibertFn -- see Bigatti 1997 (JPAA)

# Data needed:
#  convention: indets are called x[1] to x[n] -- start from 1 because julia does
#  (only at topmost call) Weight matrix:  col k is weight of x[k]
#  PP:  internal repr is expv  Vector{PP_exponent}  equiv to Vector{Int} or similar
#  Gens:  interreduced non-empty list of PPs
#  HS "indets":  seems useful to have power-prods of them by the corr weights.

# Case gens are simple powers
function HSNum_base_SimplePowers(SimplePPs::Vector{PP}, T::Vector{HSNumVar}) # T is list of HSNum PPs, one for each grading dim
    ans = 1; #one(T[1])    ???
    for t in SimplePPs #do
        k = findfirst(entry -> (entry > 0), t.expv);
        ans = ans * (1 - T[k]^t.expv[k]); # ???? ans -= ans*T[k]^t;
    end #for
    return ans;
end #function


function HSNum_base_case1(t::PP, SimplePPs::Vector{PP}, T::Vector{HSNumVar}) # T is list of HSNum PPs, one for each grading dim
#    println("HSNum_base_case1: t = $(t)");
#    println("HSNum_base_case1: SimplePPs = $(SimplePPs)");
    # t is not a "simple power", all the others are
    ans = HSNum_base_SimplePowers(SimplePPs, T);
###    println("HSNum_base_case1: init ans = $(ans)");
    ReducedSimplePPs::Vector{PP} = []; # Vector{PP}
    for j in 1:length(SimplePPs) #do
        @inbounds e = SimplePPs[j].expv;
        k = findfirst((entry -> (entry > 0)), e); # ispositive
        if  t.expv[k] == 0 #=then=# push!(ReducedSimplePPs,SimplePPs[j]); continue; end #=if=#  # no need to make a copy
        tt = copy(SimplePPs[j]);
        tt.expv[k] -= t.expv[k]; # guaranteed > 0
##        println("in loop:  SimplePPs = $(SimplePPs)");
        push!(ReducedSimplePPs, tt);
    end #for
#    println("HSNum_base_case1: input  SimplePPs = $(SimplePPs)");
#    println("HSNum_base_case1: ReducedSimplePPs = $(ReducedSimplePPs)");
    e = t.expv;
    nvars = length(e);
    @inbounds scale = prod([T[k]^e[k]  for k in 1:nvars]);
###    println("HSNum_base_case1: scale = $(scale)");
    ans = ans - scale * HSNum_base_SimplePowers(ReducedSimplePPs, T)
###println("HSNum_base_case1: final ans = $(ans)");
    return ans;
end # function


## CC contains at least 2 connected components (each compt repr as Vector{Int64} of the variable indexes in the compt)
function HSNum_splitting_case(CC::Vector{Vector{Int64}}, SimplePPs::Vector{PP}, NonSimplePPs::Vector{PP},  T::Vector{HSNumVar}, PivotStrategy::Symbol)
    HSNumList = []; ## list of HSNums
    # Now find any simple PPs which are indep of the conn compts found
    nvars = length(NonSimplePPs[1].expv);
    FoundVars = PP([0  for _ in 1:nvars]); # will be prod of all vars appearing in some conn compt
    for IndexSet in CC #=do=#
        for j in IndexSet #=do=#
            mult_by_var!(FoundVars, j);
        end #=for=#
    end #=for=#
    IsolatedSimplePPs = filter((t -> is_coprime(t,FoundVars)), SimplePPs);
    for IndexSet in CC #do
        SubsetNonSimplePPs = filter((t -> involves(t,IndexSet)), NonSimplePPs);
        SubsetSimplePPs = filter((t -> involves(t,IndexSet)), SimplePPs);
##println(" -- SPLIT --  recursive call to LOOP  Simple=$(SubsetSimplePPs)    NonSimple=$(SubsetGens)");
        push!(HSNumList, HSNum_loop(SubsetSimplePPs, SubsetNonSimplePPs, T, PivotStrategy));
#?            SubsetGens = [ProjectIndets(t, IndexSet)  for t in SubsetGens];
#?            SubsetSimplePPs = [ProjectIndets(t, IndexSet)  for t in SubsetSimplePPs];
#?            push!(HSNumList, HSNum_loop(SubsetSimplePPs, SubsetGens, [T[k]  for k in IndexSet], PivotStrategy));
    end #for
    HSNum_combined = prod(HSNumList);
    if !isempty(IsolatedSimplePPs) #then
        HSNum_combined *= HSNum_loop(IsolatedSimplePPs, PP[], T, PivotStrategy);
##            HSNum_combined *= prod([HSNum_loop([t],PP[],T,PivotStrategy)  for t in IsolatedSimplePPs]);
    end #if
    return HSNum_combined;
end #=function=#


function HSNum_total_splitting_case(VarIndexes::Vector{Int64}, SimplePPs::Vector{PP}, NonSimplePPs::Vector{PP},  T::Vector{HSNumVar}, PivotStrategy::Symbol)
    HSNumList = []; ## list of HSNums
    # Now find any simple PPs which are indep of the conn compts found
    nvars = length(NonSimplePPs[1].expv);
    FoundVars = PP([0  for _ in 1:nvars]); # will be prod of all vars appearing in some conn compt
    for i in VarIndexes #=do=#
        mult_by_var!(FoundVars, i);
    end #=for=#
    IsolatedSimplePPs = filter((t -> is_coprime(t,FoundVars)), SimplePPs);
    for t in NonSimplePPs #do
        SubsetSimplePPs = filter((s -> !is_coprime(s,t)), SimplePPs);
##println(" -- SPLIT --  recursive call to LOOP  Simple=$(SubsetSimplePPs)    NonSimple=$(SubsetGens)");
        push!(HSNumList, HSNum_base_case1(t, SubsetSimplePPs, T));
    end #for
    HSNum_combined = prod(HSNumList);
    if !isempty(IsolatedSimplePPs) #then
        HSNum_combined *= HSNum_loop(IsolatedSimplePPs, PP[], T, PivotStrategy);
##            HSNum_combined *= prod([HSNum_loop([t],PP[],T)  for t in IsolatedSimplePPs]);
    end #if
    return HSNum_combined;
end #=function=#

# term-order corr to matrix with -1 on anti-diag: [[0,0,-1],[0,-1,0],...]
function IsRevLexSmaller(t1::PP, t2::PP)
    n = length(t1.expv);
    for j in n:-1:1 #=do=#
        if t1.expv[j] != t2.expv[j]  #=then=# return (t1.expv[j] > t2.expv[j]); end #=if=#
    end #=for=#
    return false; # t1 and t2 were equal (should not happen in this code)
end #=function=#

function RevLexMin(L::Vector{PP})
    # assume length(L) > 0
    if isempty(L) #=then=# return L[1];  end #=if=#
    IndexMin = 1;
    for j in 2:length(L) #=do=#
        if !IsRevLexSmaller(L[j], L[IndexMin])
            IndexMin = j;
        end #=if=#
    end #=for=#
    return L[IndexMin];
end #=function=#

function RevLexMax(L::Vector{PP})
    # assume length(L) > 0
    if length(L) == 1 #=then=# return L[1];  end #=if=#
    IndexMax = 1;
    for j in 2:length(L) #=do=#
        if !IsRevLexSmaller(L[IndexMax], L[j])
            IndexMax = j;
        end #=if=#
    end #=for=#
    return L[IndexMax];
end #=function=#

function HSNum_Bayer_Stillman(SimplePPs::Vector{PP}, NonSimplePPs::Vector{PP},  T::Vector{HSNumVar})
##println("HSNum_BS: Simple:    $(SimplePPs)");
##println("HSNum_BS: NonSimple: $(NonSimplePPs)");
    # Maybe sort the gens???
    if isempty(NonSimplePPs)  #=then=# return HSNum_base_SimplePowers(SimplePPs, T); end #=if=#
    if length(NonSimplePPs) == 1 #=then=# return HSNum_base_case1(NonSimplePPs[1], SimplePPs, T); end #=if=#
    # NonSimplePPs contains at least 2 elements
#    BSPivot = last(NonSimplePPs); # pick one somehow -- this is just simple to impl
#??    BSPivot = RevLexMin(NonSimplePPs); # VERY SLOW on Hilbert-test-rnd6.jl
    BSPivot = RevLexMax(NonSimplePPs);
#VERY SLOW?!?    BSPivot = RevLexMax(vcat(SimplePPs,NonSimplePPs)); # VERY SLOW on Hilbert-test-rnd6.jl
println("BSPivot = $(BSPivot)");
    NonSPP = filter((t -> (t != BSPivot)), NonSimplePPs);
    SPP = SimplePPs;#filter((t -> (t != BSPivot)), SimplePPs);
    part1 = HSNum_loop(SPP, NonSPP, T, :bayer_stillman);
    ReducedPPs = interreduce([colon(t,BSPivot)  for t in vcat(SPP,NonSPP)]);
    NewSimplePPs, NewNonSimplePPs = SeparateSimplePPs(ReducedPPs);
    part2 = HSNum_loop(NewSimplePPs, NewNonSimplePPs, T, :bayer_stillman);
    e = BSPivot.expv;
    return part1 - prod([T[k]^e[k]  for k in 1:length(e)])*part2;
end #=function=#

#--------------------------------------------
# Pivot selection strategies

# [[AUXILIARY FN]]
# Return 2-tuple:
#   (part 1) list of indexes of the indets which appear in
#            the greatest number of PPs in gens.
#   (part 2) corr freq
function HSNum_most_freq_indets(gens::Vector{PP})
    # ASSUMES: gens is non-empty
    nvars = length(gens[1].expv);
    freq = [0  for i in 1:nvars];  # Vector{Int}  or Vector{UInt}  ???
    for t in gens #do
        e = t.expv;
        for i in 1:nvars #do
            @inbounds if (e[i] != 0) #then
                @inbounds freq[i] += 1
            end #if
        end #for
    end #for
    MaxFreq = maximum(freq);
###    if MaxFreq == 1  #=then=# println("MaxFreq = 1"); end #=if=#
    MostFreq = findall((x -> x==MaxFreq), freq);
    return MostFreq, MaxFreq;
end #=function=#


# Returns index of the indet
function HSNum_most_freq_indet1(gens::Vector{PP})
    # ASSUMES: gens is non-empty
    MostFreq,_ = HSNum_most_freq_indets(gens);
    return MostFreq[1];
end #=function=#

# Returns index of the indet
function HSNum_most_freq_indet_rnd(gens::Vector{PP})
    # ASSUMES: gens is non-empty
    MostFreq,_ = HSNum_most_freq_indets(gens);
    return rand(MostFreq);
end #=function=#



function HSNum_choose_pivot_indet(MostFreq::Vector{Int64}, gens::Vector{PP})
    PivotIndet = rand(MostFreq); ##HSNum_most_freq_indet_rnd(gens); # or HSNum_most_freq_indet1
    nvars = length(gens[1].expv);
    PivotExpv = [0  for _ in 1:nvars];
    PivotExpv[PivotIndet] = 1;
    PivotPP = PP(PivotExpv);
end #function

function HSNum_choose_pivot_simple_power_median(MostFreq::Vector{Int64}, gens::Vector{PP})
    # simple-power-pivot from Bigatti JPAA, 1997
    PivotIndet = rand(MostFreq); ##HSNum_most_freq_indet_rnd(gens); # or HSNum_most_freq_indet1
    exps = [t.expv[PivotIndet]  for t in gens];
    exps = filter((e -> e>0), exps);
    sort!(exps);
    exp = exps[div(1+length(exps),2)];  # "median"
    nvars = length(gens[1].expv);
    PivotExpv = [0  for _ in 1:nvars];
    PivotExpv[PivotIndet] = exp;
    PivotPP = PP(PivotExpv);
end #function

function HSNum_choose_pivot_simple_power_max(MostFreq::Vector{Int64}, gens::Vector{PP})
    # simple-power-pivot from Bigatti JPAA, 1997
    PivotIndet = rand(MostFreq); ##HSNum_most_freq_indet_rnd(gens); # or HSNum_most_freq_indet1
    exps = [t.expv[PivotIndet]  for t in gens];
    exp = max(exps...);
    nvars = length(gens[1].expv);
    PivotExpv = [0  for _ in 1:nvars];
    PivotExpv[PivotIndet] = exp;
    PivotPP = PP(PivotExpv);
end #function

function HSNum_choose_pivot_gcd2simple(MostFreq::Vector{Int64}, gens::Vector{PP})
    PivotIndet = rand(MostFreq); ##???HSNum_most_freq_indet_rnd(gens); # or HSNum_most_freq_indet1
    cand = filter((t -> t.expv[PivotIndet]>0), gens);
    if length(cand) == 1 #=then=#
        println(">>>>>> SHOULD NEVER HAPPEN <<<<<<"); error("!OUCH!");
        return cand[1]; # can happen only if there is a "simple splitting case"
    end #=if=#
    nvars = length(gens[1].expv);
    expv = [0  for _ in 1:nvars];
    expv[PivotIndet] = min(cand[1].expv[PivotIndet],  cand[2].expv[PivotIndet]);
    return PP(expv);
end #function

function HSNum_choose_pivot_gcd2max(MostFreq::Vector{Int64}, gens::Vector{PP})
    PivotIndet = rand(MostFreq); ##???HSNum_most_freq_indet_rnd(gens); # or HSNum_most_freq_indet1
    cand = filter((t -> t.expv[PivotIndet]>0), gens);
    if length(cand) == 1 #=then=#
        println(">>>>>> SHOULD NEVER HAPPEN <<<<<<"); error("!OUCH!");
        return cand[1]; # can happen only if there is a "simple splitting case"
    end #=if=#
    pick2 = [cand[k]  for k in random_subset(length(cand),2)];
    t = gcd(pick2[1], pick2[2]);
#    println("BEFORE: t= $(t)");
    d = 0; for i in 1:length(t.expv) #=do=# d = max(d,t.expv[i]); end #=for=#
    for i in 1:length(t.expv)  #=do=# if t.expv[i] < d #=then=#  t.expv[i]=0; end #=if=#  end #=for=#
#    println("AFTER: t= $(t)");
    return t;
end #function


# May produce a non-simple pivot!!!
function HSNum_choose_pivot_gcd3(MostFreq::Vector{Int64}, gens::Vector{PP})
    PivotIndet = rand(MostFreq); ##???HSNum_most_freq_indet_rnd(gens); # or HSNum_most_freq_indet1
    cand = filter((t -> t.expv[PivotIndet]>0), gens);
    if length(cand) == 1 #=then=#
        println(">>>>>> SHOULD NEVER HAPPEN <<<<<<"); error("!OUCH!");
        return cand[1]; # can happen only if there is a "simple splitting case"
    end #=if=#
    if length(cand) == 2 #=then=#
        return gcd(cand[1], cand[2]);
    end #=if=#
    pick3 = [cand[k]  for k in random_subset(length(cand),3)];
    return gcd3(pick3[1], pick3[2], pick3[3]);
end #function


function HSNum_choose_pivot_gcd3simple(MostFreq::Vector{Int64}, gens::Vector{PP})
    PivotIndet = rand(MostFreq); ##???HSNum_most_freq_indet_rnd(gens); # or HSNum_most_freq_indet1
    cand = filter((t -> t.expv[PivotIndet]>0), gens);
    if length(cand) == 1 #=then=#
        println(">>>>>> SHOULD NEVER HAPPEN <<<<<<"); error("!OUCH!");
        return cand[1]; # can happen only if there is a "simple splitting case"
    end #=if=#
    if length(cand) == 2 #=then=#
        t = gcd(cand[1], cand[2]);
    else
    pick3 = [cand[k]  for k in random_subset(length(cand),3)];
    t = gcd3(pick3[1], pick3[2], pick3[3]);
    end #=if=#
#    println("BEFORE: t= $(t)");
    j = 1;
    d = t.expv[1];
    for i in 2:length(t.expv) #=do=#
        if t.expv[i] <= d  #=then=#
            t.expv[i] = 0;
            continue;
        end #=if=#
        t.expv[j] = 0;
        j=i;
        d = t.expv[i];
    end #=for=#
#    println("AFTER: t= $(t)");
    return t;
end #function

function HSNum_choose_pivot_gcd3max(MostFreq::Vector{Int64}, gens::Vector{PP})
    PivotIndet = rand(MostFreq); ##???HSNum_most_freq_indet_rnd(gens); # or HSNum_most_freq_indet1
    cand = filter((t -> t.expv[PivotIndet]>0), gens);
    if length(cand) == 1 #=then=#
        println(">>>>>> SHOULD NEVER HAPPEN <<<<<<"); error("!OUCH!");
        return cand[1]; # can happen only if there is a "simple splitting case"
    end #=if=#
    if length(cand) == 2 #=then=#
        t = gcd(cand[1], cand[2]);
    else
    pick3 = [cand[k]  for k in random_subset(length(cand),3)];
    t = gcd3(pick3[1], pick3[2], pick3[3]);
    end #=if=#
#    println("BEFORE: t= $(t)");
    d = 0; for i in 1:length(t.expv) #=do=# d = max(d,t.expv[i]); end #=for=#
    for i in 1:length(t.expv)  #=do=# if t.expv[i] < d #=then=#  t.expv[i]=0; end #=if=#  end #=for=#
#    println("AFTER: t= $(t)");
    return t;
end #function


# May produce a non-simple pivot!!!
function HSNum_choose_pivot_gcd4(MostFreq::Vector{Int64}, gens::Vector{PP})
    PivotIndet = rand(MostFreq); ##???HSNum_most_freq_indet_rnd(gens); # or HSNum_most_freq_indet1
    cand = filter((t -> t.expv[PivotIndet]>0), gens);
    if length(cand) == 1 #=then=#
        println(">>>>>> SHOULD NEVER HAPPEN <<<<<<"); error("!OUCH!");
        return cand[1]; # can happen only if there is a "simple splitting case"
    end #=if=#
    if length(cand) == 2 #=then=#
        return gcd(cand[1], cand[2]);
    end #=if=#
    if length(cand) ==3 #=then=#
        return gcd3(cand[1], cand[2], cand[3]);
    end #=if=#
    pick4 = [cand[k]  for k in random_subset(length(cand),4)];
    return gcd(gcd(pick4[1], pick4[2]), gcd(pick4[3],pick4[4]));
end #function



# Assume SimplePPs+NonSimplePPs are interreduced
function HSNum_loop(SimplePPs::Vector{PP}, NonSimplePPs::Vector{PP},  T::Vector{HSNumVar}, PivotStrategy::Symbol)
# println("LOOP: SimplePPs=$(SimplePPs)");
# println("LOOP: first <=5 NonSimplePPs=$(first(NonSimplePPs,5))");
# CHECK=vcat(SimplePPs,NonSimplePPs);if length(CHECK) != length(interreduce(CHECK)) #=then=# println("LOOP: !!!BAD INPUT!!!"); error("BANG!"); end #=if=#
    # Check if we have base case 0
    if  isempty(NonSimplePPs)  #then
        # println("LOOP:  END OF CALL --> delegate base case 0");
        return HSNum_base_SimplePowers(SimplePPs, T);
    end #if
    # Check if we have base case 1
    if  length(NonSimplePPs) == 1  #then
        # println("LOOP:  END OF CALL --> delegate base case 1");
        return HSNum_base_case1(NonSimplePPs[1], SimplePPs, T);
    end #if
    # ----------------------
    MostFreq,freq = HSNum_most_freq_indets(NonSimplePPs);
    if  freq == 1 #=then=#
        # total splitting case
#print("+");
        return HSNum_total_splitting_case(MostFreq, SimplePPs, NonSimplePPs, T, PivotStrategy);
    end #=if=#
    
    if PivotStrategy == :bayer_stillman  #=then=#
        return HSNum_Bayer_Stillman(SimplePPs, NonSimplePPs, T);
    end #=if=#
    # Check for "splitting case"
if #=length(NonSimplePPs[1].expv) > 4 &&=# length(NonSimplePPs) <= length(NonSimplePPs[1].expv)#=nvars=# #then
    CC = ConnectedComponents(NonSimplePPs);
else CC = []
    end #if
    if #=false &&=# length(CC) > 1 #then
###        print("*");#println("LOOP: END OF CALL  --> delegate SPLITTING CASE  $(CC)");
        return HSNum_splitting_case(CC, SimplePPs, NonSimplePPs, T, PivotStrategy);
    end #if  (splitting case)
    # ----------------------
    # Pivot case: first do the ideal sum, then do ideal quotient
    # (the ideas are relatively simple, the code is long, tedious, and a bit fiddly)
#     # REDUCE BIG EXPS -- THIS IDEA SEEMS TO WORK DISAPPONTINGLY ON THE RND TEST EXAMPLES
#     nvars = length(NonSimplePPs[1].expv);
#     TopPP = PP([0  for _ in 1:nvars]);
#     for t in NonSimplePPs #=do=# TopPP = lcm(TopPP, t); end #=for=#
# DegBound = 9; BigExp = max(TopPP.expv...);
#     if false && BigExp > DegBound #=then=#
#         expv = [0  for _ in 1:nvars];
#         for i in 1:nvars #=do=#  if TopPP.expv[i] == BigExp #=then=# expv[i] = div(1+TopPP.expv[i],2); break; end #=if=#
# ##NO        for i in 1:nvars #=do=#  if TopPP.expv[i] > 9 #=then=# expv[i] = div(1+TopPP.expv[i],2); #=break;=# end #=if=#
#         end #=for=#
#         PivotPP = PP(expv);
# #println("BIG pivot $(PivotPP)");
#     else
#    PivotPP = HSNum_choose_pivot_gcd2simple(MostFreq, NonSimplePPs)
    if PivotStrategy == :indet
        PivotPP = HSNum_choose_pivot_indet(MostFreq, NonSimplePPs)
end #=if=#
    if PivotStrategy == :simple_power_max
        PivotPP = HSNum_choose_pivot_simple_power_max(MostFreq, NonSimplePPs)
end #=if=#
    if PivotStrategy == :gcd2simple
        PivotPP = HSNum_choose_pivot_gcd2simple(MostFreq, NonSimplePPs)
end #=if=#
    if PivotStrategy == :gcd2max
        PivotPP = HSNum_choose_pivot_gcd2max(MostFreq, NonSimplePPs)
end #=if=#
    if PivotStrategy == :gcd3
        PivotPP = HSNum_choose_pivot_gcd3(MostFreq, NonSimplePPs)
end #=if=#
    if PivotStrategy == :gcd3simple
        PivotPP = HSNum_choose_pivot_gcd3simple(MostFreq, NonSimplePPs)
end #=if=#
    if PivotStrategy == :gcd3max
        PivotPP = HSNum_choose_pivot_gcd3max(MostFreq, NonSimplePPs)
end #=if=#
    if PivotStrategy == :gcd4
        PivotPP = HSNum_choose_pivot_gcd4(MostFreq, NonSimplePPs)
end #=if=#
    if PivotStrategy == :simple_power_median  || PivotStrategy == :auto
        PivotPP = HSNum_choose_pivot_simple_power_median(MostFreq, NonSimplePPs)
    end #=if=#
#     end #=if=#
##=DEVELOP/DEBUG=# if (degree(PivotPP) > 3 && !IsSimplePowerPP(PivotPP))   println(">>>>>PIVOT<<<<<  $(PivotPP)"); end;
    PivotIsSimple = IsSimplePowerPP(PivotPP);
    PivotIndex = findfirst((e -> e>0), PivotPP.expv); # used only if PivotIsSimple == true
    USE_SAFE_VERSION_SUM=false;
    if  USE_SAFE_VERSION_SUM  #=then=#
        # Safe but slow version: just add new gen, then interreduce
        RecurseSum = vcat(SimplePPs, NonSimplePPs);
        push!(RecurseSum, PivotPP);
        RecurseSum = interreduce(RecurseSum);
        RecurseSum_SimplePPs_OLD, RecurseSum_NonSimplePPs_OLD = SeparateSimplePPs(RecurseSum);
    RecurseSum_NonSimplePPs = RecurseSum_NonSimplePPs_OLD;
    RecurseSum_SimplePPs = RecurseSum_SimplePPs_OLD;
    else # use "clever" version:
        # We know that SimplePPs + NonSimplePPs are already interreduced
        if  PivotIsSimple #=then=#
            RecurseSum_SimplePPs_NEW = copy(SimplePPs);
            k = findfirst((t -> t.expv[PivotIndex] > 0), SimplePPs);
            if k === nothing
                push!(RecurseSum_SimplePPs_NEW, PivotPP);
            else
                RecurseSum_SimplePPs_NEW[k] = PivotPP;
            end #if
            RecurseSum_NonSimplePPs_NEW = filter((t -> t.expv[PivotIndex] < PivotPP.expv[PivotIndex]), NonSimplePPs);
        else # PivotPP is not simple -- so this is the "general case"
            RecurseSum_SimplePPs_NEW = copy(SimplePPs); # need to copy?
            RecurseSum_NonSimplePPs_NEW = filter((t -> !IsDivisible(t,PivotPP)), NonSimplePPs);
            push!(RecurseSum_NonSimplePPs_NEW, PivotPP);
        end #=if=#
    RecurseSum_NonSimplePPs = RecurseSum_NonSimplePPs_NEW;
    RecurseSum_SimplePPs = RecurseSum_SimplePPs_NEW;
    end #=if=# # USE_SAFE_VERSION_SUM
    # println("      SUM recurses on:");
    # println("        SimplePPs    = $(RecurseSum_SimplePPs)");
    # println("        NonSimplePPs = $(RecurseSum_NonSimplePPs)");

    # Check that both approaches gave equivalent answers (well, not obviously inequivalent)
    # # if length(RecurseSum_NonSimplePPs_NEW) != length(RecurseSum_NonSimplePPs_OLD) || length(RecurseSum_SimplePPs_NEW) != length(RecurseSum_SimplePPs_OLD) #=then=#
    # #     println("!!!DIFFER!!!");
    # #     println("PivotPP = $(PivotPP)");
    # #     println("Simple_NEW = $(RecurseSum_SimplePPs_NEW)");
    # #     println("Simple_OLD = $(RecurseSum_SimplePPs_OLD)");
    # #     println("NonSimple_NEW = $(RecurseSum_NonSimplePPs_NEW)");
    # #     println("NonSimple_OLD = $(RecurseSum_NonSimplePPs_OLD)");
    # # end #=if=#
 
    # Now do the quotient...
    # Now get SimplePPs & NonSimplePPs for the quotient while limiting amount of interreduction
    USE_SAFE_VERSION_QUOT = false;
    if USE_SAFE_VERSION_QUOT #=then=#
        #=SAFE VERSION: simpler but probably slower=#
        RecurseQuot = [colon(t,PivotPP)  for t in vcat(SimplePPs,NonSimplePPs)];
        RecurseQuot = interreduce(RecurseQuot);
        RecurseQuot_SimplePPs, RecurseQuot_NonSimplePPs = SeparateSimplePPs(RecurseQuot);
    else  # Use "smart version"
        if !PivotIsSimple #then  # GENERAL CASE (e.g. if not PivotIsSimple)
            # Clever approach (for non-simple pivots) taken from Bigatti 1997 paper (p.247 just after Prop 1 to end of sect 5)
            # println("QUOT: NON-SIMPLE PIVOT  $(PivotPP)");
            # println("     init SimplePPs = $(SimplePPs)");
            # println("     init NonSimplePPs = $(NonSimplePPs)");
            BM = PP[]; NotBM = PP[];
            PivotPlus = mult(PivotPP, radical(PivotPP));
            for t in NonSimplePPs #=do=# if IsDivisible(t,PivotPlus) #=then=# push!(BM,t); else push!(NotBM, t); end #=if=# end #=for=#
            ##        println("PivotPP is $(PivotPP)");
            ##        println("PivotPlus is $(PivotPlus)");
            # println("     NotBM is $(NotBM)");
            BM = [divide(t,PivotPP)  for t in BM];  # divide is same as colon here
            NotBM = vcat(NotBM, SimplePPs);
            NotBM_mixed = PP[];  NotBM_coprime = PP[];
            for t in NotBM #=do=# if is_coprime(t,PivotPP) #=then=# push!(NotBM_coprime,t); else push!(NotBM_mixed, colon(t,PivotPP)); end #=if=# end #=for=#
            # At this poiint we have 3 disjoint lists of PPs: BM (big multiples), NotBM_coprime, NotBM_mixed
            # In all cases the PPs have been colon-ed by PivotPP
            # println("     NotBM_mixed = $(NotBM_mixed)");
            # println("     NotBM_coprime = $(NotBM_coprime)");
            NotBM_mixed = interreduce(NotBM_mixed); # cannot easily be "clever" here
            # println("     NotBM_mixed INTERRED = $(NotBM_mixed)");
            filter!((t -> NotMultOf(NotBM_mixed,t)), NotBM_coprime);
            # println("     NotBM_coprime FILTERED is $(NotBM_coprime)");
            ##        if !isempty(BM)  println("BM is $(BM)"); else print("*"); end
            RecurseQuot = vcat(NotBM_coprime, NotBM_mixed); # already interreduced
##DEBUGGING            if length(RecurseQuot) != length(interreduce(RecurseQuot)) #=then=# println(">>>>>>DIFFER<<<<<<  RQ=$(length(RecurseQuot))   interr=$(length(interreduce(RecurseQuot)))"); end #=if=#
            RQ_SimplePPs, RQ_NonSimplePPs = SeparateSimplePPs(RecurseQuot);
            RecurseQuot_SimplePPs = RQ_SimplePPs; ###???interreduce(vcat(RQ_SimplePPs, [colon(t,PivotPP)  for t in SimplePPs]));
            RecurseQuot_NonSimplePPs = vcat(BM, RQ_NonSimplePPs);
            ###        RecurseQuot = vcat(NotBM, [colon(t,PivotPP)  for t in SimplePPs]);  RecurseQuot = interreduce(RecurseQuot); RecurseQuot_SimplePPs, RecurseQuot_NonSimplePPs = SeparateSimplePPs(RecurseQuot);   RecurseQuot_NonSimplePPs = vcat(RecurseQuot_NonSimplePPs, BM);  # Be cleverer about interreduced
            # println("QUOT FINAL:  PivotPP=$(PivotPP)");
            # println("QUOT FINAL:  RecurseQuot_NonSimplePPs=$(RecurseQuot_NonSimplePPs)");
            # println("QUOT FINAL:  RecurseQuot_SimplePPs=$(RecurseQuot_SimplePPs)");
##CHECK=vcat(RecurseQuot_SimplePPs,RecurseQuot_NonSimplePPs);if length(CHECK) != length(interreduce(CHECK)) #=then=# println("LOOP: QUOT(S) FINAL !!!BAD OUTPUT!!!"); error("BANG!"); end #=if=#
        else  # Clever approach when PivotIsSimple
            # println("QUOT: SIMPLE PIVOT  $(PivotPP)");
            # println("      init SimplePPs = $(SimplePPs)");
            # println("      init NonSimplePPs = $(NonSimplePPs)");
            # The idea behind this code is fairly simple; sadly the code itself is not :-(
            RecurseQuot_SimplePPs = copy(SimplePPs);
            k = findfirst((t -> t.expv[PivotIndex] > 0), RecurseQuot_SimplePPs);
            if !(k === nothing)
                RecurseQuot_SimplePPs[k] = copy(RecurseQuot_SimplePPs[k]);  #???copy needed???
                RecurseQuot_SimplePPs[k].expv[PivotIndex] -= PivotPP.expv[PivotIndex];
            end #if
            DegPivot = degree(PivotPP);
            NonSimpleTbl = [PP[]  for _ in 0:DegPivot];  ## WARNING: indexes are offset by 1 -- thanks Julia!
            NonSimple1 = PP[]; # will contain all PPs divisible by PivotPP^(1+epsilon)
            for t in NonSimplePPs  #=do=#
                degt = t.expv[PivotIndex];
                if degt > DegPivot #=then=#
                    push!(NonSimple1, divide(t, PivotPP));
                else
                    push!(NonSimpleTbl[degt+1], colon(t,PivotPP));
                end #=if=#
            end #=for=#
            ##        println("PivotPP = $(PivotPP)");
            # println("(Quot simple) TBL: $(NonSimpleTbl)");
            NonSimple2 = NonSimpleTbl[DegPivot+1];
            for i in DegPivot:-1:1 #=do=#
                NewPPs = filter((t -> NotMultOf(NonSimple2,t)), NonSimpleTbl[i]);
                NonSimple2 = vcat(NonSimple2, NewPPs);
            end #=for=#
            ##        println("After interred: NonSimple1=$(NonSimple1)");
            ##        println("After interred: NonSimple2=$(NonSimple2)");
            NewSimplePPs = filter(IsSimplePowerPP, NonSimple2);
            NonSimple2 = filter(!IsSimplePowerPP, NonSimple2);  ## SeparateSimplePPs???
            if !isempty(NewSimplePPs) #=then=#
                RecurseQuot_SimplePPs = interreduce(vcat(RecurseQuot_SimplePPs, NewSimplePPs));
            end #=if=#
            RecurseQuot_NonSimplePPs = vcat(NonSimple1, NonSimple2);
            # println("QUOT(S) FINAL:  PivotPP=$(PivotPP)");
            # println("QUOT(S) FINAL:  RecurseQuot_NonSimplePPs=$(RecurseQuot_NonSimplePPs)");
            # println("QUOT(S) FINAL:  RecurseQuot_SimplePPs=$(RecurseQuot_SimplePPs)");
##CHECK=vcat(RecurseQuot_SimplePPs,RecurseQuot_NonSimplePPs);if length(CHECK) != length(interreduce(CHECK)) #=then=# println("LOOP: QUOT(S) FINAL !!!BAD OUTPUT!!!"); error("BANG!"); end #=if=#
        end #=if=# #PivotIsSimple
    end #=if=# # end of USE_SAFE_VERSION_QUOT
    # Now put the two pieces together:
    nvars = length(PivotPP.expv);
    scale = prod([T[k]^PivotPP.expv[k]  for k in 1:nvars]);
    # println("RECURSION:");
    # println("   SUM recursion:  simple    $(RecurseSum_SimplePPs)");
    # println("   SUM recursion:  nonsimple $(RecurseSum_NonSimplePPs)");
    # println("   QUOT recursion: simple    $(RecurseQuot_SimplePPs)");
    # println("   QUOT recursion: nonsimple $(RecurseQuot_NonSimplePPs)");
    HSNum_sum = HSNum_loop(RecurseSum_SimplePPs, RecurseSum_NonSimplePPs, T, PivotStrategy);
    # println("RECURSION: after HSNum_sum");
    # println("   SUM recursion:  simple    $(RecurseSum_SimplePPs)");
    # println("   SUM recursion:  nonsimple $(RecurseSum_NonSimplePPs)");
    # println("   QUOT recursion: simple    $(RecurseQuot_SimplePPs)");
    # println("   QUOT recursion: nonsimple $(RecurseQuot_NonSimplePPs)");
    HSNum_quot = HSNum_loop(RecurseQuot_SimplePPs, RecurseQuot_NonSimplePPs, T, PivotStrategy);
    # println("LOOP:  END OF CALL");
    return HSNum_sum + scale*HSNum_quot;
end #function

                  
function SeparateSimplePPs(gens::Vector{PP})
    SimplePPs::Vector{PP} = [];
    NonSimplePPs::Vector{PP} = [];
    for g in gens  #do
        if IsSimplePowerPP(g)  #then
            push!(SimplePPs, g)
        else
            push!(NonSimplePPs, g)
        end #if
    end #for
    return SimplePPs, NonSimplePPs;
end #function


# Check args: either throws or returns nothing.
function HSNum_CheckArgs(gens::Vector{PP}, W::Vector{Vector{Int}})
    if isempty(gens)  throw("HSNum: need at least 1 generator");  end #if
    if isempty(W)  throw("HSNum: weight matrix must have at least 1 row");  end #if
    nvars = length(gens[1].expv);
    if !all((t -> length(t.expv)==nvars), gens)
        throw("HSNum: generators must all have same size exponent vectors");
    end #if
    if !all((row -> length(row)==nvars), W)
        throw("HSNum: weight matrix must have 1 column for each variable")
    end #if
    # Zero weights are allowed???
    # Args are OK, so simply return (without throwing)
end # function


function HSNum(gens::Vector{PP}, W::Vector{Vector{Int}}, PivotStrategy::Symbol)
###    println("HSNum: gens = $(gens)");
###    println("HSNum: W = $(W)");
    HSNum_CheckArgs(gens, W);
    # Grading is over ZZ^m
    m = size(W)[1];  # NumRows
    ncols = size(W[1])[1];  # how brain-damaged is Julia???
    nvars = length(gens[1].expv);
    if  ncols != nvars #then
        throw(ArgumentError("weights matrix has wrong number of columns ($(ncols)); should be same as number of variables ($(nvars))"));
    end #if
    HPRing, t = LaurentPolynomialRing(QQ, ["t$k"  for k in 1:m]);
    T = [one(HPRing)  for k in 1:nvars];
    for k in 1:nvars #do
        s = one(HPRing);
        for j in 1:m #do
            s *= t[j]^W[j][k];
        end #for
        T[k] = s;
    end #for
    # Now have T[1] = t1^W[1,1] * t2^W[2,1] * ...,  etc
###                  println("HSNum: T vector is $(T)");
    SimplePPs,NonSimplePPs = SeparateSimplePPs(interreduce(gens));
    sort!(NonSimplePPs, lt=DegRevLexLess); # recommended by Bayer+Stillman (for their criterion)
    return HSNum_loop(SimplePPs, NonSimplePPs, T, PivotStrategy);
end #function

#-----------------------------------------------------------------------------
# This fn copied from GradedModule.jl (in dir OSCAR/HILBERT/)
function gen_repr(d)
    grading_dim = length(gens(parent(d)));
    return [getindex(d,k)  for k in 1:grading_dim];
end #function

function HSNum_fudge(PmodI::MPolyQuoRing, PivotStrategy::Symbol = :auto)
if PivotStrategy == :indet   return nothing; end
    I = PmodI.I;
    P = base_ring(I);##parent(gens(I)[1]); # there MUST be a better way!!
    nvars = length(gens(P));
    grading_dim = length(gens(parent(degree(gen(P,1))))); # better way???
    weights = [degree(var)  for var in gens(P)];
    W = [[0  for _ in 1:nvars]  for _ in 1:grading_dim];
    for i in 1:nvars #do
        expv = [Int64(exp)  for exp in gen_repr(degree(gen(P,i)))];
        for j in 1:grading_dim #do
            W[j][i] = expv[j];
        end #for
    end #for
    # ##    W = [[Int64(exp)  for exp in gen_repr(d)]  for d in weights];
    # W=[]
    # for d in weights #do
    #     expv = [Int64(exp)  for exp in gen_repr(d)];
    #     if isempty(W) #then
    #         W = expv;
    #     else
    #         W = hcat(W, expv);
    #     end #if
    # end #for
    # W = [W[:,i] for i in 1:size(W,2)]
    # # ?transpose?  hcat
##    println("W is $(W)");
    LTs = gens(leading_ideal(I));
    PPs = [PP(degrees(t))  for t in LTs];
    return HSNum(PPs, W, PivotStrategy);
end #function
