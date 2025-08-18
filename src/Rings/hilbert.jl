##################################################################
# Auxiliary functions
# Perhaps document them and make them publicly available?

# random_subset:
# A function similar to this is in StatsBase.jl but I wish to avoid
# having an external dependency (for just 1 simple fn!)

# Generate random m-subset of {1,2,3,...,n}
# Result is a Int64[]; entries are NOT SORTED!!
function random_subset(n::Int64, m::Int64)
  # assume n >= 1, m >= 0 and m <= n
  if m == 0
    return Int64[]
  end
  L = collect(1:n)
  if m == n
    return L
  end
  for j in 1:m
    k = rand(j:n)
    L[j],L[k] = L[k],L[j] # just a SWAP
  end
  L = first(L,m)
  #?? sort!(L) ??
  return L
end

##################################################################
# Code for representing & manipulating PPs (power products, aka. monomials)
# All this code is "local" to this file, & not exported!


# Each PP is represented as Vector{PP_exponent}

const PP_exponent = Int64 # FIXME  # UInt ???  Strange: Int32 was slower on my machine ?!?


struct PP
  expv::Vector{PP_exponent}
end

function Base.copy(t::PP)
  return PP(copy(t.expv))
end

function length(t::PP)
  return length(t.expv)
end

function getindex(t::PP, i::Int)
  return t.expv[i]
end

function setindex!(t::PP, i::Int, d::Int)
  return t.expv[i] = d
end

# Should be is_one, but julia complained :-(
function isone(t::PP)
  return all(t.expv .== 0)
end

function degree(t::PP)
  return sum(t.expv)
end


# RETURN VALUE???  Perhaps index or 0? (or -1?)
function is_simple_power_pp(t::PP)
  CountNZ = 0
  for i in 1:length(t)
    @inbounds t[i] == 0 && continue
    CountNZ > 0 && return false
    CountNZ = i
  end
  CountNZ != 0 && return true # MAYBE RETURN index & exp???
  return false # because t == 1
end


function is_divisible(t::PP, s::PP)  # is t divisible by s
  nvars = length(t) # assume equal to length(s)
  for i in 1:nvars
    @inbounds t[i] < s[i] && return false
  end
  return true
end


# modifies first arg
function mult_by_var!(t::PP, j::Int64)
  @inbounds t.expv[j] += 1
end

function mult(t1::PP, t2::PP)
  # ASSUMES: length(t1) == length(t2)
  return PP(t1.expv + t2.expv)
end


function divide(t1::PP, t2::PP)
  # ASSUMES: length(t1) == length(t2), also that t1 is mult of t2
  return PP(t1.expv - t2.expv)
end

function is_coprime(t1::PP, t2::PP)
  # ASSUMES: length(t1) == length(t2)
  nvars = length(t1)
  for i in 1:nvars
    @inbounds if t1[i] != 0 && t2[i] != 0
      return false
    end
  end
  return true
end

function lcm(t1::PP, t2::PP)
  # ASSUMES: length(t1.expv) == length(t2.expv)
  nvars = length(t1)
  expv = [0  for _ in 1:nvars]
  for i in 1:nvars
    @inbounds expv[i] = max(t1[i], t2[i])
  end
  return PP(expv)
end

function gcd(t1::PP, t2::PP)
  # ASSUMES: length(t1.expv) == length(t2.expv)
  nvars = length(t1)
  expv = [0  for _ in 1:nvars]
  for i in 1:nnvars
    @inbounds expv[i] = min(t1[i], t2[i])
  end
  return PP(expv)
end

function gcd3(t1::PP, t2::PP, t3::PP)
  # ASSUMES: length(t1) == length(t2) == length(t3)
  nvars = length(t1)
  expv = [0  for _ in 1:nvars]
  for i in 1:nvars
    @inbounds expv[i] = min(t1[i], t2[i], t3[i])
  end
  return PP(expv)
end


# Computes t1/gcd(t1,t2)
function colon(t1::PP, t2::PP)
  # ASSUMES: length(t1) == length(t2)
  nvars = length(t1)
  expv = [0  for _ in 1:nvars]
  for i in 1:nvars
    @inbounds expv[i] = max(0, t1[i]-t2[i])
  end
  return PP(expv)
end


# Computes t1/gcd(t1,t2^infty)
function saturatePP(t1::PP, t2::PP)
  # ASSUMES: length(t1) == length(t2)
  nvars = length(t1)
  expv = [0  for _ in 1:nvars]
  for i in 1:nvars
    @inbounds if t2[i] == 0
      @inbounds expv[i] = t1[i]
    end
  end
  return PP(expv)
end

# Computes radical = product of vars which divide t
function radical(t::PP)
  nvars = length(t)
  expv = [0  for _ in 1:nvars]
  for i in 1:nvars
    @inbounds if t[i] > 0
      expv[i] = 1
    end
  end
  return PP(expv)
end


# Test if t1 < t1 in DegRevLex ordering;
# NB if t1 == t2 then returns false.
function deg_rev_lex_less(t1::PP, t2::PP)
  d1 = degree(t1)
  d2 = degree(t2)
  d1 != d2 && return d1 < d2
  nvars = length(t1)
  for i in nvars:-1:1
    t1[i] != t2[i] && return t1[i] > t2[i]
  end
  return false
end


# Test whether PP involves at least 1 indet with index in IndexList
function involves(t::PP, IndexList::Vector{Int64})
  # ASSUMES: indexes in IndexList are all in range
  for i in IndexList
    @inbounds t[i] != 0 && return true
  end
  return false
end


function Base.show(io::IO, t::PP)
  if all(t.expv .== 0)   # isone(PP)  # why doesn't this work ???
    print(io, "1")
    return
  end
  str = ""
  nvars = length(t)
  for i in 1:nvars
    if (t[i] != 0)
      if (str != "")
        str = str * "*"
      end
      str = str * "x[$(i)]"
      if (t[i] > 1)
        str = str * "^$(t[i])"
      end
    end
  end
  print(io, str)
end


#######################################################
# Interreduction of a list of PPs

# interreduce list of PPs; equiv find min set of gens for monoideal gen by L
function interreduce(L::Vector{PP})
  sort!(L, by=degree)
  MinGens = PP[]
  for t in L
    discard = false
    for s in MinGens
      if is_divisible(t,s)
        discard = true
        break
      end
    end
    if !discard
      push!(MinGens, t)
    end
  end
  return MinGens
end


# Is t a multiple of at least one element of L?
# Add degree truncation???
function not_mult_of_any(L::Vector{PP}, t::PP)
  for s in L
    if is_divisible(t,s)
      return false
    end
  end
  return true
end



# "project" PP onto sub-monoid of PPs gen by indets in indexes
function project_indets(t::PP, indexes::Vector{Int})
  return PP([t[k]  for k in indexes])
end


# NOT SURE THIS IS USEFUL: how many indets are really needed in the list L?
function true_num_vars(L::Vector{PP})
  isempty(L) && return 0
  t = radical(L[1])
  for j in 2:length(L)
    t = lcm(t, radical(L[j]))
  end
  return degree(t)
end


###################################################################
# This function is also private/local to this file.

# Think of a graph where each vertex is labeled by an indeterminate.
# Each PP in the list L is interpreted as saying that all indeterminates
# involved in that PP are "connected".  The task is to find (minimal)
# connected components of the entire graph.

# Input: non-empty list of PPs
# Output: list of lists of var indexes, each sublist is a connected component
function connected_components(L::Vector{PP})
  ConnCompt = Vector{Vector{Int64}}()
  nvars = length(L[1])
  IgnoreVar = [ false   for _ in 1:nvars]
  VarAppears = copy(IgnoreVar)
  for t in L
    for j in 1:nvars
      @inbounds if t[j] != 0
        @inbounds VarAppears[j] = true
      end
    end
  end
  CountIgnore = 0
  for j in 1:nvars
    @inbounds if !VarAppears[j]
      @inbounds IgnoreVar[j] = true
      CountIgnore += 1
    end
  end
  while CountIgnore < nvars
    j = findfirst(!, IgnoreVar) # j is index of some var which appears in at least one PP
    k = findfirst((t -> t[j] != 0), L) # pick some PP involving j-th var
    lcm = L[k]
    DoAnotherIteration = true
    while DoAnotherIteration
      DoAnotherIteration = false
      for t in L
        is_coprime(lcm,t) && continue
        s = saturatePP(t,lcm)
        isone(s) && continue
        lcm = mult(lcm,s) ### lcm *= s
        DoAnotherIteration = true
      end
    end #while
    vars = filter((k -> lcm[k] > 0), 1:nvars)
    # remove conn compt just found from L???
    #seems to be slower with this line ?!?        L = filter((t -> is_coprime(t,lcm)), L)
    push!(ConnCompt, vars)
    for k in vars
      IgnoreVar[k] = true
      CountIgnore += 1
    end
  end #while
  return ConnCompt
end


#############################################################################
# Implementation of "Hilbert series numerator":
# main reference Bigatti 1997 (JPAA) "Computation of Hilbert-Poincare series"
#-----------------------------------------------------------------------------

# Two base cases for HibertFn

# Data needed:
#  convention:      indets are called x[1] to x[n] -- start from 1 because julia does
#  Weight matrix:   col k is weight of x[k]  (only at topmost call)
#  PP:              internal repr is expv  Vector{PP_exponent}  equiv to Vector{Int} or similar
#  Gens:
#     SimplePPs:    PPs which are "simple" power of an indet
#     NonSimplePPs: PPs which are divisible by at least 2 distinct indets
#                   SimplePPs union NonSimplePPs is interreduced & non-empty!!
#  HS "indets":     seems useful to have power-prods of them by the corr weights -- parameter T in the fn defns

# Case gens are simple powers
function HSNum_base_SimplePowers(SimplePPs::Vector{PP}, T::Vector{RingElemType}) where {RingElemType <: RingElem} # T is list of HSNum PPs, one for each grading dim
  ans = one(T[1])
  for t in SimplePPs
    k = findfirst(>(0), t.expv)
    ans = ans * (1 - T[k]^t[k]) # ???? ans -= ans*T[k]^t[k]
  end
  return ans
end


function HSNum_base_case1(t::PP, SimplePPs::Vector{PP}, T::Vector{RingElemType}) where {RingElemType <: RingElem} # T is list of HSNum PPs, one for each grading dim
  # t is not a "simple power", all the others are
  @vprintln  :hilbert 1  "HSNum_base_case1: t = $(t)"
  @vprintln  :hilbert 1  "HSNum_base_case1: SimplePPs = $(SimplePPs)"
  ans = HSNum_base_SimplePowers(SimplePPs, T)
  ReducedSimplePPs = PP[]
  for spp in SimplePPs
    e = spp.expv
    k = findfirst(>(0), e) # ispositive
    if  t[k] == 0
      push!(ReducedSimplePPs, spp)
      continue
    end   # no need to make a copy
    tt = copy(spp)
    tt.expv[k] -= t[k] # guaranteed > 0
    push!(ReducedSimplePPs, tt)
  end
  nvars = length(t)
  @inbounds scale = prod([T[k]^t[k]  for k in 1:nvars])
  ans = ans - scale * HSNum_base_SimplePowers(ReducedSimplePPs, T)
  @vprintln  :hilbert 1  "HSNum_base_case1: returning $(ans)"
  return ans
end # function


## CC contains at least 2 connected components (each compt repr as Vector{Int64} of the variable indexes in the compt)
function HSNum_splitting_case(CC::Vector{Vector{Int64}}, SimplePPs::Vector{PP}, NonSimplePPs::Vector{PP},  T::Vector{RET}, PivotStrategy::Symbol) where {RET <: RingElem}
  @vprintln  :hilbert 1  "Splitting case: CC = $(CC)"
  HSNumList = [] ## list of HSNums
  # Now find any simple PPs which are indep of the conn compts found
  nvars = length(NonSimplePPs[1])
  FoundVars = PP([0  for _ in 1:nvars]) # will be prod of all vars appearing in some conn compt
  for IndexSet in CC
    for j in IndexSet
      mult_by_var!(FoundVars, j)
    end
  end
  IsolatedSimplePPs = filter((t -> is_coprime(t,FoundVars)), SimplePPs)
  for IndexSet in CC
    SubsetNonSimplePPs = filter((t -> involves(t,IndexSet)), NonSimplePPs)
    SubsetSimplePPs = filter((t -> involves(t,IndexSet)), SimplePPs)
    push!(HSNumList, HSNum_loop(SubsetSimplePPs, SubsetNonSimplePPs, T, PivotStrategy))
    # Next 3 lines commented out: seemed to bring no benefit (??but why??)
    #?            SubsetNonSimplePPs = [project_indets(t, IndexSet)  for t in SubsetNonSimplePPs]
    #?            SubsetSimplePPs = [project_indets(t, IndexSet)  for t in SubsetSimplePPs]
    #?            push!(HSNumList, HSNum_loop(SubsetSimplePPs, SubsetGens, [T[k]  for k in IndexSet], PivotStrategy))
  end
  HSNum_combined = prod(HSNumList)
  if !isempty(IsolatedSimplePPs)
    HSNum_combined *= HSNum_base_SimplePowers(IsolatedSimplePPs, T)
    ##OLD        HSNum_combined *= HSNum_loop(IsolatedSimplePPs, PP[], T, PivotStrategy)
  end
  return HSNum_combined
end


function HSNum_total_splitting_case(VarIndexes::Vector{Int64}, SimplePPs::Vector{PP}, NonSimplePPs::Vector{PP},  T::Vector{RET}, PivotStrategy::Symbol) where {RET <: RingElem}
  @vprintln  :hilbert 1 "Total splitting case: VarIndexes = $(VarIndexes)"
  HSNumList = [] ## list of HSNums
  # Now find any simple PPs which are indep of the conn compts found
  nvars = length(NonSimplePPs[1])
  FoundVars = PP([0  for _ in 1:nvars]) # will be prod of all vars appearing in some conn compt
  for i in VarIndexes
    mult_by_var!(FoundVars, i)
  end
  IsolatedSimplePPs = filter((t -> is_coprime(t,FoundVars)), SimplePPs)
  for t in NonSimplePPs
    SubsetSimplePPs = filter((s -> !is_coprime(s,t)), SimplePPs)
    push!(HSNumList, HSNum_base_case1(t, SubsetSimplePPs, T))
  end
  HSNum_combined = prod(HSNumList)
  if !isempty(IsolatedSimplePPs)
    HSNum_combined *= HSNum_base_SimplePowers(IsolatedSimplePPs, T)
    ##OLD        HSNum_combined *= HSNum_loop(IsolatedSimplePPs, PP[], T, PivotStrategy)
  end
  return HSNum_combined
end


#------------------------------------------------------------------
# Pivot selection strategies: (see Bigatti 1997 paper)

# term-order corr to matrix with -1 on anti-diag: [[0,0,-1], [0,-1,0],...]
function is_rev_lex_smaller(t1::PP, t2::PP)
  n = length(t1)
  for j in n:-1:1
    t1[j] != t2[j] && return t1[j] > t2[j]
  end
  return false # t1 and t2 were equal (should not happen in this code)
end

function rev_lex_min(L::Vector{PP})
  # assume length(L) > 0
  if isempty(L)  return L[1]  end
  IndexMin = 1
  for j in 2:length(L)
    if !is_rev_lex_smaller(L[j], L[IndexMin])
      IndexMin = j
    end
  end
  return L[IndexMin]
end

function rev_lex_max(L::Vector{PP})
  # assume length(L) > 0
  if length(L) == 1  return L[1]  end
  IndexMax = 1
  for j in 2:length(L)
    if !is_rev_lex_smaller(L[IndexMax], L[j])
      IndexMax = j
    end
  end
  return L[IndexMax]
end

function HSNum_bayer_stillman(SimplePPs::Vector{PP}, NonSimplePPs::Vector{PP},  T::Vector{RET}) where {RET <: RingElem}
  ##println("HSNum_BS: Simple:    $(SimplePPs)")
  ##println("HSNum_BS: NonSimple: $(NonSimplePPs)")
  # Maybe sort the gens???
  if isempty(NonSimplePPs)
    return HSNum_base_SimplePowers(SimplePPs, T)
  end
  if length(NonSimplePPs) == 1
    return HSNum_base_case1(NonSimplePPs[1], SimplePPs, T)
  end
  # NonSimplePPs contains at least 2 elements
  #    BSPivot = last(NonSimplePPs) # pick one somehow -- this is just simple to impl
  #??    BSPivot = rev_lex_min(NonSimplePPs) # VERY SLOW on Hilbert-test-rnd6.jl
  BSPivot = rev_lex_max(NonSimplePPs)
  #VERY SLOW?!?    BSPivot = rev_lex_max(vcat(SimplePPs,NonSimplePPs)) # VERY SLOW on Hilbert-test-rnd6.jl
  @vprintln :hilbert 2 "BSPivot = $(BSPivot)"
  NonSPP = filter(!=(BSPivot), NonSimplePPs)
  SPP = SimplePPs#filter((t -> (t != BSPivot)), SimplePPs)
  part1 = HSNum_loop(SPP, NonSPP, T, :bayer_stillman)
  ReducedPPs = interreduce([colon(t,BSPivot)  for t in vcat(SPP,NonSPP)])
  NewSimplePPs, NewNonSimplePPs = SeparateSimplePPs(ReducedPPs)
  part2 = HSNum_loop(NewSimplePPs, NewNonSimplePPs, T, :bayer_stillman)
  nvars = length(BSPivot)
  return part1 - prod([T[k]^BSPivot[k]  for k in 1:nvars])*part2
end

#--------------------------------------------
# Pivot selection strategies

# [[AUXILIARY FN]]
# Return 2-tuple:
#   (part 1) list of indexes of the indets which appear in
#            the greatest number of PPs in gens.
#   (part 2) corr freq
function HSNum_most_freq_indets(gens::Vector{PP})
  # ASSUMES: gens is non-empty
  nvars = length(gens[1])
  freq = [0  for i in 1:nvars]  # Vector{Int}  or Vector{UInt}  ???
  for t in gens
    for i in 1:nvars
      @inbounds if (t[i] != 0)
        @inbounds freq[i] += 1
      end
    end
  end
  MaxFreq = maximum(freq)
  MostFreq = findall(==(MaxFreq), freq)
  return MostFreq, MaxFreq
end


# Returns index of the indet
function HSNum_most_freq_indet1(gens::Vector{PP})
  # ASSUMES: gens is non-empty
  MostFreq,_ = HSNum_most_freq_indets(gens)
  return MostFreq[1]
end

# Returns index of the indet
function HSNum_most_freq_indet_rnd(gens::Vector{PP})
  # ASSUMES: gens is non-empty
  MostFreq,_ = HSNum_most_freq_indets(gens)
  return rand(MostFreq)
end



# THIS PIVOT STRATEGY WAS AWFULLY SLOW!
# function HSNum_choose_pivot_indet(MostFreq::Vector{Int64}, gens::Vector{PP})
#   PivotIndet = rand(MostFreq) ##HSNum_most_freq_indet_rnd(gens) # or HSNum_most_freq_indet1
#   nvars = length(gens[1])
#   PivotExpv = [0  for _ in 1:nvars]
#   PivotExpv[PivotIndet] = 1
#   PivotPP = PP(PivotExpv)
# end

function HSNum_choose_pivot_simple_power_median(MostFreq::Vector{Int64}, gens::Vector{PP})
  # variant of simple-power-pivot from Bigatti JPAA, 1997
  PivotIndet = rand(MostFreq)
  exps = [t[PivotIndet]  for t in gens]
  exps = filter(>(0), exps)
  sort!(exps)
  exp = exps[div(1+length(exps),2)]  # "median"
  nvars = length(gens[1])
  PivotExpv = [0  for _ in 1:nvars]
  PivotExpv[PivotIndet] = exp
  PivotPP = PP(PivotExpv)
end

function HSNum_choose_pivot_simple_power_max(MostFreq::Vector{Int64}, gens::Vector{PP})
  # variant of simple-power-pivot from Bigatti JPAA, 1997
  PivotIndet = rand(MostFreq)
  exps = [t[PivotIndet]  for t in gens]
  exp = max(exps...)
  nvars = length(gens[1])
  PivotExpv = [0  for _ in 1:nvars]
  PivotExpv[PivotIndet] = exp
  PivotPP = PP(PivotExpv)
end

function HSNum_choose_pivot_gcd2simple(MostFreq::Vector{Int64}, gens::Vector{PP})
  # simple-power-pivot from Bigatti JPAA, 1997
  PivotIndet = rand(MostFreq)
  cand = filter((t -> t[PivotIndet]>0), gens)
  if length(cand) == 1
    error("Failed to detect total splitting case") ### !!SHOULD NEVER GET HERE!!
  end
  nvars = length(gens[1])
  expv = [0  for _ in 1:nvars]
  expv[PivotIndet] = min(cand[1].expv[PivotIndet],  cand[2].expv[PivotIndet])
  return PP(expv)
end

function HSNum_choose_pivot_gcd2max(MostFreq::Vector{Int64}, gens::Vector{PP})
  PivotIndet = rand(MostFreq)
  cand = filter((t -> t[PivotIndet]>0), gens)
  if length(cand) == 1
    error("Failed to detect total splitting case") ### !!SHOULD NEVER GET HERE!!
  end
  pick2 = [cand[k]  for k in random_subset(length(cand),2)]
  t = gcd(pick2...)
  nvars = length(t)
  d = maximum(t.expv)
  for i in 1:nvars
    if t[i] < d
      t.expv[i]=0
    end
  end
  return t
end


# # May produce a non-simple pivot!!!
# function HSNum_choose_pivot_gcd3(MostFreq::Vector{Int64}, gens::Vector{PP})
#   PivotIndet = rand(MostFreq)
#   cand = filter((t -> t[PivotIndet]>0), gens)
#   if length(cand) == 1
#     error("Failed to detect total splitting case") ### !!SHOULD NEVER GET HERE!!
#   end
#   if length(cand) == 2
#     return gcd(cand[1], cand[2])
#   end
#   pick3 = [cand[k]  for k in random_subset(length(cand),3)]
#   return gcd3(pick3[1], pick3[2], pick3[3])
# end


function HSNum_choose_pivot_gcd3simple(MostFreq::Vector{Int64}, gens::Vector{PP})
  PivotIndet = rand(MostFreq)
  cand = filter((t -> t[PivotIndet]>0), gens)
  if length(cand) == 1
    error("Failed to detect total splitting case") ### !!SHOULD NEVER GET HERE!!
  end
  if length(cand) == 2
    t = gcd(cand[1], cand[2])
  else
    pick3 = [cand[k]  for k in random_subset(length(cand),3)]
    t = gcd3(pick3...) ##    t = gcd3(pick3[1], pick3[2], pick3[3])
  end
  # t is gcd of 3 rnd PPs (or of 2 if there are only 2)
  # Now set all exps to 0 except the first max one.
  j = 1
  d = t[1]
  for i in 2:length(t)
    if t[i] <= d
      t.expv[i] = 0
      continue
    end
    t.expv[j] = 0
    j = i
    d = t[i]
  end
  return t
end


function HSNum_choose_pivot_gcd3max(MostFreq::Vector{Int64}, gens::Vector{PP})
  PivotIndet = rand(MostFreq)
  cand = filter((t -> t[PivotIndet]>0), gens)
  if length(cand) == 1
    error("Failed to detect total splitting case") ### !!SHOULD NEVER GET HERE!!
  end
  if length(cand) == 2
    t = gcd(cand[1], cand[2])
  else
    pick3 = [cand[k]  for k in random_subset(length(cand),3)]
    t = gcd3(pick3...)  ##    t = gcd3(pick3[1], pick3[2], pick3[3])
  end
  d = maximum(t.expv)
  # Now set to 0 all exps which are less than the max
  for i in 1:length(t)
    if t[i] < d
      t.expv[i]=0
    end
  end
  return t
end


# # May produce a non-simple pivot!!!
# function HSNum_choose_pivot_gcd4(MostFreq::Vector{Int64}, gens::Vector{PP})
#   PivotIndet = rand(MostFreq)
#   cand = filter((t -> t[PivotIndet]>0), gens)
#   if length(cand) == 1
#     error("Failed to detect total splitting case") ### !!SHOULD NEVER GET HERE!!
#   end
#   if length(cand) == 2
#     return gcd(cand[1], cand[2])
#   end
#   if length(cand) ==3
#     return gcd3(cand[1], cand[2], cand[3])
#   end
#   pick4 = [cand[k]  for k in random_subset(length(cand),4)]
#   return gcd(gcd(pick4[1], pick4[2]), gcd(pick4[3],pick4[4]))
# end



# Assume SimplePPs+NonSimplePPs are interreduced
function HSNum_loop(SimplePPs::Vector{PP}, NonSimplePPs::Vector{PP},  T::Vector{RET}, PivotStrategy::Symbol) where {RET <: RingElem}
  @vprintln  :hilbert 1 "HSNum_loop: SimplePPs=$(SimplePPs)"
  @vprintln  :hilbert 1 "HSNum_loop: NonSimplePPs=$(NonSimplePPs)"
#  @vprintln  :hilbert 1 "LOOP: first <=5 NonSimplePPs=$(first(NonSimplePPs,5))"
  # Check if we have base case 0
  if  isempty(NonSimplePPs)
    @vprintln :hilbert  1 "HSNum_loop:  --> delegate base case 0"
    return HSNum_base_SimplePowers(SimplePPs, T)
  end
  # Check if we have base case 1
  if  length(NonSimplePPs) == 1
    @vprintln :hilbert  1 "HSNum_loop:  --> delegate base case 1"
    return HSNum_base_case1(NonSimplePPs[1], SimplePPs, T)
  end
  # ----------------------
  MostFreq,freq = HSNum_most_freq_indets(NonSimplePPs)
  if  freq == 1
    return HSNum_total_splitting_case(MostFreq, SimplePPs, NonSimplePPs, T, PivotStrategy)
  end

  if PivotStrategy == :bayer_stillman
    return HSNum_bayer_stillman(SimplePPs, NonSimplePPs, T)
  end
  # Check for "splitting case"
  if length(NonSimplePPs) <= length(NonSimplePPs[1])#=nvars=#
    CC = connected_components(NonSimplePPs)
  else
    CC = []
  end
  if length(CC) > 1
    @vprintln :hilbert  1 "HSNum_loop:  --> delegate Splitting case"
    return HSNum_splitting_case(CC, SimplePPs, NonSimplePPs, T, PivotStrategy)
  end
  # ----------------------
  # Pivot case: first do the ideal sum, then do ideal quotient
  # (the ideas are relatively simple, the code is long, tedious, and a bit fiddly)
  if PivotStrategy == :simple_power_median  || PivotStrategy == :auto
    PivotPP = HSNum_choose_pivot_simple_power_median(MostFreq, NonSimplePPs)
  end
  if PivotStrategy == :simple_power_max
    PivotPP = HSNum_choose_pivot_simple_power_max(MostFreq, NonSimplePPs)
  end
  if PivotStrategy == :gcd2simple
    PivotPP = HSNum_choose_pivot_gcd2simple(MostFreq, NonSimplePPs)
  end
  if PivotStrategy == :gcd2max
    PivotPP = HSNum_choose_pivot_gcd2max(MostFreq, NonSimplePPs)
  end
  # if PivotStrategy == :gcd3
  #   PivotPP = HSNum_choose_pivot_gcd3(MostFreq, NonSimplePPs)
  # end
  if PivotStrategy == :gcd3simple
    PivotPP = HSNum_choose_pivot_gcd3simple(MostFreq, NonSimplePPs)
  end
  if PivotStrategy == :gcd3max
    PivotPP = HSNum_choose_pivot_gcd3max(MostFreq, NonSimplePPs)
  end
  # if PivotStrategy == :gcd4
  #   PivotPP = HSNum_choose_pivot_gcd4(MostFreq, NonSimplePPs)
  # end
  #     end
  @vprintln :hilbert  1 "HSNum_loop:  pivot = $(PivotPP)"
  PivotIsSimple = is_simple_power_pp(PivotPP)
  PivotIndex = findfirst(>(0), PivotPP.expv) # used only if PivotIsSimple == true
  USE_SAFE_VERSION_SUM=false
  if  USE_SAFE_VERSION_SUM
    # Safe but slow version: just add new gen, then interreduce
    RecurseSum = vcat(SimplePPs, NonSimplePPs)
    push!(RecurseSum, PivotPP)
    RecurseSum = interreduce(RecurseSum)
    RecurseSum_SimplePPs_OLD, RecurseSum_NonSimplePPs_OLD = SeparateSimplePPs(RecurseSum)
    RecurseSum_NonSimplePPs = RecurseSum_NonSimplePPs_OLD
    RecurseSum_SimplePPs = RecurseSum_SimplePPs_OLD
  else # use "clever" version:
    # We know that SimplePPs + NonSimplePPs are already interreduced
    if  PivotIsSimple
      RecurseSum_SimplePPs_NEW = copy(SimplePPs)
      k = findfirst((t -> t[PivotIndex] > 0), SimplePPs)
      if k === nothing
        push!(RecurseSum_SimplePPs_NEW, PivotPP)
      else
        RecurseSum_SimplePPs_NEW[k] = PivotPP
      end
      RecurseSum_NonSimplePPs_NEW = filter((t -> t[PivotIndex] < PivotPP[PivotIndex]), NonSimplePPs)
    else # PivotPP is not simple -- so this is the "general case"
      RecurseSum_SimplePPs_NEW = copy(SimplePPs) # need to copy?
      RecurseSum_NonSimplePPs_NEW = filter((t -> !is_divisible(t,PivotPP)), NonSimplePPs)
      push!(RecurseSum_NonSimplePPs_NEW, PivotPP)
    end
    RecurseSum_NonSimplePPs = RecurseSum_NonSimplePPs_NEW
    RecurseSum_SimplePPs = RecurseSum_SimplePPs_NEW
  end  # USE_SAFE_VERSION_SUM

  # Now do the quotient...
  # Now get SimplePPs & NonSimplePPs for the quotient while limiting amount of interreduction
  USE_SAFE_VERSION_QUOT = false
  if USE_SAFE_VERSION_QUOT
    #=SAFE VERSION: simpler but probably slower=#
    RecurseQuot = [colon(t,PivotPP)  for t in vcat(SimplePPs,NonSimplePPs)]
    RecurseQuot = interreduce(RecurseQuot)
    RecurseQuot_SimplePPs, RecurseQuot_NonSimplePPs = SeparateSimplePPs(RecurseQuot)
  else  # Use "smart version"
    if !PivotIsSimple   # GENERAL CASE (e.g. if not PivotIsSimple)
      # Clever approach (for non-simple pivots) taken from Bigatti 1997 paper (p.247 just after Prop 1 to end of sect 5)
      BM = PP[]
      NotBM = PP[]
      PivotPlus = mult(PivotPP, radical(PivotPP))
      for t in NonSimplePPs
        if is_divisible(t,PivotPlus)
          push!(BM,t)
        else
          push!(NotBM,t)
        end
      end
      # BM is short for "big multiple" (see Bigatti 97, p.247 between Prop 1 and Prop 2)
      BM = [divide(t,PivotPP)  for t in BM]  # divide is same as colon here
      NotBM = vcat(NotBM, SimplePPs)
      NotBM_mixed = PP[]
      NotBM_coprime = PP[]
      for t in NotBM
        if is_coprime(t,PivotPP)
          push!(NotBM_coprime,t)
        else
          push!(NotBM_mixed, colon(t,PivotPP))
        end
      end
      # At this point we have 3 disjoint lists of PPs:
      #    BM (big multiples),
      #    NotBM_coprime,
      #    NotBM_mixed
      # In all cases the PPs have been colon-ed by PivotPP
      NotBM_mixed = interreduce(NotBM_mixed) # cannot easily be "clever" here
      filter!((t -> not_mult_of_any(NotBM_mixed,t)), NotBM_coprime)
      RecurseQuot = vcat(NotBM_coprime, NotBM_mixed) # already interreduced
      RQ_SimplePPs, RQ_NonSimplePPs = SeparateSimplePPs(RecurseQuot)
      RecurseQuot_SimplePPs = RQ_SimplePPs
      RecurseQuot_NonSimplePPs = vcat(BM, RQ_NonSimplePPs)
    else  # Clever approach when PivotIsSimple
      # The idea behind this code is fairly simple; sadly the code itself is not :-(
      RecurseQuot_SimplePPs = copy(SimplePPs)
      k = findfirst((t -> t[PivotIndex] > 0), RecurseQuot_SimplePPs)
      if !(k === nothing)
        RecurseQuot_SimplePPs[k] = copy(RecurseQuot_SimplePPs[k])
        RecurseQuot_SimplePPs[k].expv[PivotIndex] -= PivotPP[PivotIndex]
      end
      DegPivot = degree(PivotPP)
      NonSimpleTbl = [PP[]  for _ in 0:DegPivot]  ## WARNING: indexes are offset by 1 -- thanks Julia!
      NonSimple1 = PP[] # will contain all PPs divisible by PivotPP^(1+epsilon)
      for t in NonSimplePPs
        degt = t[PivotIndex]
        if degt > DegPivot
          push!(NonSimple1, divide(t, PivotPP))
        else
          push!(NonSimpleTbl[degt+1], colon(t,PivotPP))
        end
      end
      NonSimple2 = NonSimpleTbl[DegPivot+1]
      for i in DegPivot:-1:1
        NewPPs = filter((t -> not_mult_of_any(NonSimple2,t)), NonSimpleTbl[i])
        NonSimple2 = vcat(NonSimple2, NewPPs)
      end
      NewSimplePPs = filter(is_simple_power_pp, NonSimple2) ## Use instead
      NonSimple2 = filter(!is_simple_power_pp, NonSimple2)  ## SeparateSimplePPs???
      if !isempty(NewSimplePPs)
        RecurseQuot_SimplePPs = interreduce(vcat(RecurseQuot_SimplePPs, NewSimplePPs))
      end
      RecurseQuot_NonSimplePPs = vcat(NonSimple1, NonSimple2)
    end  #PivotIsSimple
  end  # USE_SAFE_VERSION_QUOT
  # Now put the two pieces together:
  nvars = length(PivotPP)
  scale = prod([T[k]^PivotPP[k]  for k in 1:nvars])
  @vprintln  :hilbert 2  "HSNum_loop:   SUM recursion:  simple    $(RecurseSum_SimplePPs)"
  @vprintln  :hilbert 2  "HSNum_loop:   SUM recursion:  nonsimple $(RecurseSum_NonSimplePPs)"
  @vprintln  :hilbert 2  "HSNum_loop:   QUOT recursion: simple    $(RecurseQuot_SimplePPs)"
  @vprintln  :hilbert 2  "HSNum_loop:   QUOT recursion: nonsimple $(RecurseQuot_NonSimplePPs)"

  HSNum_sum = HSNum_loop(RecurseSum_SimplePPs, RecurseSum_NonSimplePPs, T, PivotStrategy)
  HSNum_quot = HSNum_loop(RecurseQuot_SimplePPs, RecurseQuot_NonSimplePPs, T, PivotStrategy)
  @vprintln  :hilbert 1 "HSNum_loop:  END OF CALL"
  return HSNum_sum + scale*HSNum_quot
end


function separate_simple_pps(gens::Vector{PP})
  SimplePPs = PP[]
  NonSimplePPs = PP[]
  for g in gens
    if is_simple_power_pp(g)
      push!(SimplePPs, g)
    else
      push!(NonSimplePPs, g)
    end
  end
  return SimplePPs, NonSimplePPs
end


# !!!OBSOLESCENT!!!   2023-08-17 this fn will no be needed after Wolfram's PR is merged
# Returns nothing; throws if ker(W) contains a non-zero vector >= 0
#function _hilbert_series_check_weights(W::Vector{Vector{Int}})
#  # assumes W is rectangular (and at least 1x1)
#  # Transpose while converting:
#  ncols = length(W)
#  nrows = length(W[1])
#  A = zero_matrix(ZZ, nrows,ncols)
#  for i in 1:nrows, j in 1:ncols  A[i,j] = W[j][i]  end
#  b = zero_matrix(ZZ, nrows,1)
#  try
#    solve_non_negative(A, b) # any non-zero soln gives rise to infinitely many, which triggers an exception
#  catch e
#    if !(e isa ArgumentError && e.msg == "Polyhedron not bounded")
#      rethrow(e)
#    end
#    # solve_non_negative must have thrown because there is a non-zero soln
#    error("given weights permit infinite dimensional homogeneous spaces")
#  end
#  return nothing # otherwise it returns the result of solve_non_negative (Doh!!)
#end


# Check args: either throws or returns nothing.
function HSNum_check_args(gens::Vector{PP}, W::Vector{Vector{Int}})
  if isempty(W)
    throw("HSNum: weight matrix must have at least 1 row")
  end
  nvars = length(W[1])
  if !all((row -> length(row)==nvars), W)
    throw("HSNum: weight matrix must have 1 column for each variable")
  end
  if !all((t -> length(t)==nvars), gens)  # OK also if isempty(gens)
    throw("HSNum: generators must all have same size exponent vectors")
  end
  # Args are OK, so simply return (without throwing)
end



function HSNum_abbott_PPs(PP_gens::Vector{PP}, W::Vector{Vector{Int}}, PivotStrategy::Symbol, HSRing::Ring)
  # ASSUME W is "rectangular"
  @vprintln :hilbert 1 "HSNum: PP_gens = $(PP_gens)"
  @vprintln :hilbert 1 "HSNum: weight matrix W = $(W)"
  HSNum_check_args(PP_gens, W) #throws or does nothing
  # Grading is over ZZ^m
  m = length(W)
  ncols = length(W[1])
  nvars = ncols
  t = gens(HSRing)
  @assert  length(t) >= m  "supplied Hilbert series ring contains too few variables"
  T = [one(HSRing)  for k in 1:nvars]
  for k in 1:nvars
    s = one(HSRing)
    for j in 1:m
      s *= t[j]^W[j][k]
    end
    T[k] = s
  end
  # Now have T[1] = t1^W[1,1] * t2^W[2,1] * ...,  etc
  SimplePPs,NonSimplePPs = separate_simple_pps(interreduce(PP_gens))
  sort!(NonSimplePPs, lt=deg_rev_lex_less) # recommended by Bayer+Stillman (for their criterion)
  return HSNum_loop(SimplePPs, NonSimplePPs, T, PivotStrategy)
end

#-----------------------------------------------------------------------------
# This fn copied from GradedModule.jl (in dir OSCAR/HILBERT/)
function gen_repr(d)
  grading_dim = length(gens(parent(d)))
  return [d[k]  for k in 1:grading_dim]
end

@doc raw"""
    HSNum_abbott(A::MPolyQuoRing, HSRing::Ring; pivot_strategy::Symbol = :auto)

Don't use this internal function: use `hilbert_series` or `multi_hilbert_series` instead.
Compute numerator of Hilbert series of the quotient `A`.
Result is a pair: `N, D` being the numerator `N` (as a laurent polynomial) and the denominator `D` as
a factor list of laurent polynomials.

!!! note
    Applied to an ideal `I`, the function first homogenizes the generators of `I` in the extended ring.
    It then creates the ideal generated by these homogenizations, and saturates this ideal
    with respect to the ideal which is generated by the homogenizing variables.

# Examples
```jldoctest
julia> R, (x,y,z) = graded_polynomial_ring(QQ, [:x, :y,:z])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> I = ideal(R, [x^3+y^2*z,  y^3+x*z^2,  z^3+x^2*y]);

julia> RmodI,_ = quo(R,I);

julia> HSRing1,_ = polynomial_ring(ZZ, :t);

julia> Oscar.HSNum_abbott(RmodI, HSRing1)
-t^9 + 3*t^6 - 3*t^3 + 1

julia> R, (x,y,z) = graded_polynomial_ring(QQ, [:x, :y,:z], [1 2 3; 3 2 1])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> I = ideal(R, [x*z+y^2,  y^6+x^3*z^3,  z^6, x^6]);

julia> RmodI,_ = quo(R,I);

julia> HSRing2,_ = polynomial_ring(ZZ, :t => 1:2);

julia> Oscar.HSNum_abbott(RmodI, HSRing2)
-t[1]^28*t[2]^28 + t[1]^24*t[2]^24 + t[1]^22*t[2]^10 - t[1]^18*t[2]^6 + t[1]^10*t[2]^22 - t[1]^6*t[2]^18 - t[1]^4*t[2]^4 + 1

```
"""
function HSNum_abbott(PmodI::MPolyQuoRing, HSRing::Ring; pivot_strategy::Symbol = :auto)
  I = modulus(PmodI)
  P = base_ring(I)
  nvars = length(gens(P))
  grading_dim = length(gens(Oscar.parent(degree(gen(P,1))))) # better way???
  weights = [degree(var)  for var in gens(P)]
  W = [[0  for _ in 1:nvars]  for _ in 1:grading_dim]
  for i in 1:nvars
    expv = [Int64(exp)  for exp in gen_repr(degree(gen(P,i)))]
    for j in 1:grading_dim
      W[j][i] = expv[j]
    end
  end
  LTs = gens(leading_ideal(I))
  PPs = [PP(degrees(t))  for t in LTs]
  return HSNum_abbott_PPs(PPs, W, pivot_strategy, HSRing)
end



# W[k] is ZZ^m-grading of x[k]
# Expect that HSRing is PolynomialRing or LaurentPolynomialRing
function _hilbert_series_denominator(HSRing::Ring, W::Vector{Vector{Int}})
  n = length(W)
  m = length(W[1])
  @req  all(r -> length(r)==m, W)  "Grading VecVec must be rectangular"
  @req  length(gens(HSRing)) >= m  "Hilbert series ring has too few variables"
  t = gens(HSRing)
  fac_dict = Dict{elem_type(HSRing), Integer}()
  for i = 1:n
    # adjoin factor ( 1 - prod(t_j^W[i,j]) )
    new_fac = 1 - prod([t[k]^W[i][k]  for k in 1:m])
    # ???BUG??? DOES NOT WORK if HSRing is Univariate poly ring over ZZ
    # B = MPolyBuildCtx(HSRing)
    # push_term!(B, 1, W[i])
    # new_fac = 1-finish(B)
    if haskey(fac_dict, new_fac)
      fac_dict[new_fac] += 1  # coalesce identical factors
    else
      fac_dict[new_fac] = 1
    end
  end
  return FacElem(HSRing, fac_dict)
end


function _hilbert_series_ring(parent::Union{Nothing, Ring}, m::Int)
  if parent !== nothing
    # Caller supplied a Ring; check that it is reasonable.
    # The following line might complain in unforeseen ways in case 
    # the argument for `parent` was too unreasonable. However, we would 
    # like to keep the possibilities for that input rather broad.
    @req ngens(parent) >= m "parent ring of the output does not contain sufficiently many variables"
    return parent
  end
  # Caller did not supply a Ring, so we create one.
  if m == 1
    VAR = [:t]
  else
    VAR = [Symbol("t[$i]") for i = 1:m]
  end
  HSRing, _ = laurent_polynomial_ring(ZZ, VAR)
  return HSRing
end

# !!! Note
#   We create a Laurent poly ring even if all weights are non-negative
#   because negative exponents may be needed when dealing with modules
#   with shifts.
