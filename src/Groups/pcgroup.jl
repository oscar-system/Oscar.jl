############################################################################

# Create an Oscar collector object, its type parameter `T` describes
# the type of integers that occur as exponents.
#
# The idea is
# - to formulate the presentation in Julia,
#   and to create the corresponding GAP object on demand
#   when it gets accessed,
# - to create a GAP collector in `GAP.Globals.IsSingleCollectorRep`
#   if all relative orders are primes (in particular finite),
#   and one in `GAP.Globals.IsFromTheLeftCollectorRep` otherwise.
#
# Collectors in Julia and in GAP are mutable,
# we transfer changes in the Julia object to the GAP side
# as soon as the GAP object exists.

abstract type Collector{T} end

mutable struct GAP_Collector{T} <: Collector{T}
  ngens::Int
  relorders::Vector{T}
  powers::Vector{Vector{Pair{Int, T}}}
  conjugates::Matrix{Vector{Pair{Int, T}}}
  X::GapObj   # a collector in GAP, if defined
  F::FPGroup  # the free group in Oscar that belongs to `X`, if defined

  function GAP_Collector{T}(n::Int) where T <: IntegerUnion
    return new(n, zeros(T, n), # relative orders undefined
                  repeat([Pair{Int, T}[]], n), # default powers are identity
                  Matrix{Vector{Pair{Int, T}}}(undef, n, n)) # conjugates undefined
  end
end

"""
    collector(n::Int, ::Type{T} = ZZRingElem) where T <: IntegerUnion

Return an empty collector object for a pc group with `n` generators
and exponents of type `T`.
It describes a free abelian group of rank `n`.

# Examples
```jldoctest
julia> G = pc_group(collector(2))
Pc group of infinite order

julia> is_abelian(G)
true
```
"""
collector(n::Int, ::Type{T} = ZZRingElem) where T <: IntegerUnion = GAP_Collector{T}(n)


# Provide functions for entering data into the collector and accessing data.

# utility:
# Convert a vector of pairs to a generator-exponent vector in GAP.
function _GAP_generator_exponent_vector(l::Vector{Pair{Int, T}}) where T <: IntegerUnion
  res = T[]
  for p in l
    push!(res, p.first)
    push!(res, p.second)
  end
  return GapObj(res; recursive = true)
end

"""
    set_relative_order!(c::Collector{T}, i::Int, relord::T) where T <: IntegerUnion

Set the relative order of the `i`-th generator of `c` to `relord`,
which must be either `0` (meaning infinite order) or a positive integer.

Currently the relative orders of collectors that describe finite groups
must be primes.

# Examples
```jldoctest
julia> c = collector(2, Int);

julia> set_relative_order!(c, 1, 2)
```
"""
function set_relative_order!(c::Collector{T}, i::Int, relord::T) where T <: IntegerUnion
  @req (0 < i && i <= c.ngens) "the collector has only $(c.ngens) generators not $i"
  c.relorders[i] = relord
  all(is_positive, c.relorders) && @req all(is_prime, c.relorders) "relative orders of collectors for finite groups must be primes"

  if relord == 0
    # Remove the `i`-th power relation (as is done in GAP).
    c.powers[i] = Pair{Int, T}[]
  end

  # If the GAP collector has already been created then update it.
  if isdefined(c, :X)
    if GAP.Globals.IsSingleCollectorRep(c.X)
      error("cannot change stored relative orders of $(c.X)")
    elseif GAP.Globals.IsFromTheLeftCollectorRep(c.X)
      GAP.Globals.SetRelativeOrder(c.X, i, GAP.Obj(relord))
    else
      error("unknown GAP collector")
    end
  end
end

"""
    get_relative_order(c::Collector{T}, i::Int) where T <: IntegerUnion

Get the relative order of the `i`-th generator of `c`.

# Examples
```jldoctest
julia> c = collector(2, Int);

julia> get_relative_order(c, 1)
0

julia> set_relative_order!(c, 1, 2)

julia> get_relative_order(c, 1)
2
```
"""
function get_relative_order(c::Collector{T}, i::Int) where T <: IntegerUnion
  @req (0 < i && i <= c.ngens) "the collector has only $(c.ngens) generators not $i"
  return c.relorders[i]
end

"""
    set_relative_orders!(c::Collector{T}, relords::Vector{T})

Set all relative orders of the generators of `c`,
where the length of `relords` must be equal to the number of generators of
`c`, and `relords[i]` denotes the relative order of the `i`-th generator.
which must be either `0` (meaning infinite order) or a positive integer.

Currently the relative orders of collectors that describe finite groups
must be primes.

# Examples
```jldoctest
julia> c = collector(2);

julia> set_relative_orders!(c, ZZRingElem[2, 0])
```
"""
function set_relative_orders!(c::Collector{T}, relords::Vector{T}) where T <: IntegerUnion
  @req length(relords) == c.ngens "the collector has $(c.ngens) generators not $(length(relords))"
  all(is_positive, relords) && @req all(is_prime, relords) "relative orders of collectors for finite groups must be primes"
  copyto!(c.relorders, relords)

  # If the GAP collector has already been created then update it.
  if isdefined(c, :X)
    if GAP.Globals.IsSingleCollectorRep(c.X)
      error("cannot change stored relative orders of $(c.X)")
    elseif GAP.Globals.IsFromTheLeftCollectorRep(c.X)
      for i in 1:c.ngens
        GAP.Globals.SetRelativeOrder(c.X, i, GAP.Obj(relords[i]))
      end
    else
      error("unknown GAP collector")
    end
  end
end

"""
    get_relative_orders(c::Collector{T})

Get the `Vector{T}` of all relative orders of the generators of `c`.

# Examples
```jldoctest
julia> c = collector(2);

julia> get_relative_orders(c)
2-element Vector{ZZRingElem}:
 0
 0

julia> set_relative_orders!(c, ZZRingElem[2, 0])

julia> get_relative_orders(c)
2-element Vector{ZZRingElem}:
 2
 0
```
"""
function get_relative_orders(c::Collector{T}) where T <: IntegerUnion
  return c.relorders
end

"""
    set_power!(c::Collector{T}, i::Int, rhs::Vector{Pair{Int, T}}) where T <: IntegerUnion

Set the `c.relorders[i]`-th power of the `i`-th generator of `c`
to the group element described by `rhs`.

# Examples
```jldoctest
julia> c = collector(2, Int);

julia> set_relative_order!(c, 1, 2)

julia> set_relative_order!(c, 2, 3)

julia> set_power!(c, 1, [2 => 1])
```
"""
function set_power!(c::Collector{T}, i::Int, rhs::Vector{Pair{Int, T}}) where T <: IntegerUnion
  @req 0 < i <= c.ngens "the collector has only $(c.ngens) generators not $i"
  c.powers[i] = copy(rhs)

  # If the GAP collector has already been created then update it.
  if isdefined(c, :X)
    if GAP.Globals.IsSingleCollectorRep(c.X)
      G = c.F
      GAP.Globals.SetPower(c.X, i, G(c.powers[i]).X)
    elseif GAP.Globals.IsFromTheLeftCollectorRep(c.X)
      GAP.Globals.SetPower(c.X, i, _GAP_generator_exponent_vector(rhs))
    else
      error("unknown GAP collector")
    end
  end
end

"""
    get_power(c::Collector{T}, i::Int) where T <: IntegerUnion

Get the `Vector{Pair{Int, T}}` that describes the `c.relorders[i]`-th power
of the `i`-th generator of `c`.

# Examples
```jldoctest
julia> c = collector(2, Int);

julia> set_relative_order!(c, 1, 2)

julia> set_relative_order!(c, 2, 3)

julia> get_power(c, 1)
Pair{Int64, Int64}[]

julia> set_power!(c, 1, [2 => 1])

julia> get_power(c, 1)
1-element Vector{Pair{Int64, Int64}}:
 2 => 1
```
"""
function get_power(c::Collector{T}, i::Int) where T <: IntegerUnion
  @req 0 < i <= c.ngens "the collector has only $(c.ngens) generators not $i"
  return c.powers[i]
end

"""
    set_conjugate!(c::Collector{T}, j::Int, i::Int, rhs::Vector{Pair{Int, T}}) where T <: IntegerUnion

Set the value of the conjugate of the `j`-th generator of `c` by the
`i`-th generator of `c`, for `i < j`, to the group element
described by `rhs`.

# Examples
```jldoctest
julia> c = collector(2, Int);

julia> set_relative_orders!(c, [2, 3])

julia> set_conjugate!(c, 2, 1, [2 => 2])
```
"""
function set_conjugate!(c::Collector{T}, j::Int, i::Int, rhs::Vector{Pair{Int, T}}) where T <: IntegerUnion
  @req 0 < i <= c.ngens "the collector has only $(c.ngens) generators not $i"
  @req i < j "only for i < j, but i = $i, j = $j"
  c.conjugates[i,j] = copy(rhs)

  # If the GAP collector has already been created then update it.
  if isdefined(c, :X)
    if GAP.Globals.IsSingleCollectorRep(c.X)
      G = c.F
      GAP.Globals.SetConjugate(c.X, j, i, G(rhs).X)
    elseif GAP.Globals.IsFromTheLeftCollectorRep(c.X)
      GAP.Globals.SetConjugate(c.X, j, i, _GAP_generator_exponent_vector(rhs))
    else
      error("unknown GAP collector")
    end
  end
end

"""
    get_conjugate(c::Collector{T}, j::Int, i::Int) where T <: IntegerUnion

Get the `Vector{Pair{Int, T}}` that describes the conjugate of the `j`-th
generator of `c` by the `i`-th generator of `c`, for `i < j`.

# Examples
```jldoctest
julia> c = collector(2, Int);

julia> set_relative_orders!(c, [2, 3])

julia> get_conjugate(c, 2, 1)
1-element Vector{Pair{Int64, Int64}}:
 2 => 1

julia> set_conjugate!(c, 2, 1, [2 => 2])

julia> get_conjugate(c, 2, 1)
1-element Vector{Pair{Int64, Int64}}:
 2 => 2
```
"""
function get_conjugate(c::Collector{T}, j::Int, i::Int) where T <: IntegerUnion
  @req 0 < i <= c.ngens "the collector has only $(c.ngens) generators not $i"
  @req i < j "only for i < j, but i = $i, j = $j"
  conj = c.conjugates
  return isassigned(conj, i, j) ? conj[i, j] : [j => T(1)]
end

"""
    set_commutator!(c::Collector{T}, j::Int, i::Int, rhs::Vector{Pair{Int, T}}) where T <: IntegerUnion

Set the value of the commutator of the `i`-th and the `j`-th generator of `c`,
for `i < j`, to the group element described by `rhs`.

# Examples
```jldoctest
julia> c = collector(2, Int);

julia> set_relative_orders!(c, [2, 3])

julia> set_commutator!(c, 2, 1, [2 => 1])
```
"""
function set_commutator!(c::Collector{T}, j::Int, i::Int, rhs::Vector{Pair{Int, T}}) where T <: IntegerUnion
  @req 0 < j <= c.ngens "the collector has only $(c.ngens) generators not $j"
  @req i < j "only for i < j, but i = $i, j = $j"
  if length(rhs) > 0 && rhs[1].first == j
    # freely reduce
    e = rhs[1].second + 1
    if e != 0
      rhs = deepcopy(rhs)
      rhs[1] = j => e
    else
      rhs = rhs[2:end]
    end
  else
    rhs = vcat([j => 1], rhs)
  end
  set_conjugate!(c, j, i, rhs)
end

# Create a collector independent of `c`.
function Base.deepcopy_internal(c::GAP_Collector{T}, dict::IdDict) where T <: IntegerUnion
  cc = GAP_Collector{T}(c.ngens,
         deepcopy_internal(c.relorders, dict::IdDict),
         deepcopy_internal(c.powers, dict::IdDict), 
         deepcopy_internal(c.conjugates, dict::IdDict))

  # If the GAP collector has already been created then copy it.
  if isdefined(c, :X)
    if GAP.Globals.IsSingleCollectorRep(c.X)
      # For this type of collectors, `GAP.Globals.ShallowCopy`
      # returns an independent copy.
      cc.X = GAP.Globals.ShallowCopy(c.X)::GapObj
    elseif GAP.Globals.IsFromTheLeftCollectorRep(c.X)
      # Currently there is no `GAP.Globals.ShallowCopy` method for `c.X`.
      # Create the GAP object anew.
      cc.X = _GAP_collector_from_the_left(c)::GapObj
    else
      error("unknown GAP collector")
    end
  end

  return cc
end

# Create a new GAP collector using `GAP.Globals.SingleCollector`.
function _GAP_single_collector(c::GAP_Collector)
  @req all(is_prime, c.relorders) "the single collector requires prime relative orders"
  G = free_group(c.ngens; eltype = :syllable)
  cGAP = GAP.Globals.SingleCollector(GapObj(G), GapObj(c.relorders; recursive = true))::GapObj
  for i in 1:c.ngens
    GAP.Globals.SetPower(cGAP, i, G(c.powers[i]).X)
  end
  for j in 1:c.ngens
    for i in 1:(j-1)
      if isassigned(c.conjugates, i, j)
        GAP.Globals.SetConjugate(cGAP, j, i, G(c.conjugates[i, j]).X)
      else
        GAP.Globals.SetConjugate(cGAP, j, i, G([j => 1]).X)
      end
    end
  end
  GAP.Globals.UpdatePolycyclicCollector(cGAP)

  return cGAP, G
end

# Create a new GAP collector using `GAP.Globals.FromTheLeftCollector`
# (from GAP's Polycyclic package).
function _GAP_collector_from_the_left(c::GAP_Collector)
  cGAP = GAP.Globals.FromTheLeftCollector(c.ngens)::GapObj
  for i in 1:c.ngens
    if c.relorders[i] != 0
      # We are not allowed to set relative orders to zero.
      GAP.Globals.SetRelativeOrder(cGAP, i, GAP.Obj(c.relorders[i]))
    end
  end
  for i in 1:c.ngens
    if c.relorders[i] != 0
      GAP.Globals.SetPower(cGAP, i, _GAP_generator_exponent_vector(c.powers[i]))
    end
  end
  for j in 1:c.ngens
    for i in 1:(j-1)
      if isassigned(c.conjugates, i, j)
        GAP.Globals.SetConjugate(cGAP, j, i, _GAP_generator_exponent_vector(c.conjugates[i, j]))
      else
        GAP.Globals.SetConjugate(cGAP, j, i, GapObj([j, 1]))
      end
    end
  end
  GAP.Globals.UpdatePolycyclicCollector(cGAP)
  
  return cGAP::GapObj
end

# Create the collector on the GAP side on demand
function underlying_gap_object(c::GAP_Collector)
  if ! isdefined(c, :X)
    # We have to decide which type of GAP collector we create.
#TODO: Here we could specify the desired collector type.
    if 0 in c.relorders
      # create a collector from the Polycyclic package.
      cGAP = _GAP_collector_from_the_left(c)
    else
      # create a GAP library collector for finite polycyclic groups.
      cGAP, F = _GAP_single_collector(c)
      c.F = F
    end
    c.X = cGAP
  end
  return c.X
end

"""
    pc_group(c::GAP_Collector)

Create a `PcGroup` object from `c`.

# Examples
```jldoctest
julia> c = collector(2, Int);

julia> Oscar.set_relative_orders!(c, [2, 3])

julia> Oscar.set_conjugate!(c, 2, 1, [2 => 2])

julia> gg = pc_group(c)
Pc group of order 6

julia> describe(gg)
"S3"
```
"""
function pc_group(c::GAP_Collector)
  # Create the GAP collector if necessary.
  cGAP = underlying_gap_object(c)::GapObj

  if GAP.Globals.IsFromTheLeftCollectorRep(cGAP)
    # Create an independent collector object on the GAP side,
    # such that later changes to `c` do not affect the collector
    # that is stored in the group object on the GAP side.
    cGAP = _GAP_collector_from_the_left(c)
    return PcGroup(GAP.Globals.PcpGroupByCollector(cGAP)::GapObj)
  elseif GAP.Globals.IsSingleCollectorRep(cGAP)
    # `GAP.Globals.GroupByRws` makes the group independent of the
    # collector data.
    return PcGroup(GAP.Globals.GroupByRws(cGAP)::GapObj)
  else
    error("unknown collector type")
  end
end

"""
    exponent_vector(::Type{T} = ZZRingElem, g::Union{PcGroupElem,SubPcGroupElem})

Return the exponent vector of `g` as an instance of the type `Vector{T}`,
each entry corresponding to a group generator.

# Examples

```jldoctest
julia> g = abelian_group(PcGroup, [0, 5])
Pc group of infinite order

julia> x = g[1]^-3 * g[2]^-3
g1^-3*g2^2

julia> exponent_vector(x)
2-element Vector{ZZRingElem}:
 -3
 2
```

```jldoctest
julia> gg = small_group(6, 1)
Pc group of order 6

julia> x = gg[1]^5*gg[2]^-4
f1*f2^2

julia> exponent_vector(x)
2-element Vector{ZZRingElem}:
 1
 2
```
"""
exponent_vector(g::Union{PcGroupElem,SubPcGroupElem}) = exponent_vector(ZZRingElem, g)

function exponent_vector(
  ::Type{T}, g::Union{PcGroupElem,SubPcGroupElem}
) where {T<:IntegerUnion}
  # check if we have a PcpGroup element
  gObj = GapObj(g)
  if GAPWrap.IsPcpElement(gObj)
    return Vector{T}(GAPWrap.Exponents(gObj))
  else # finite PcGroup
    pcgs = GAPWrap.FamilyPcgs(GapObj(parent(g)))
    return Vector{T}(GAPWrap.ExponentsOfPcElement(pcgs, gObj))
  end
end

"""
    relative_order(::Type{T} = ZZRingElem, g::Union{PcGroupElem,SubPcGroupElem})

Return the relative order of `g` as an instance of the type `T`,
with respect to the defining generators. For generators with infinite order, we return 0.

# Examples

```jldoctest
julia> g = abelian_group(PcGroup, [0, 5])
Pc group of infinite order

julia> x = g[1]^-3 * g[2]^-3
g1^-3*g2^2

julia> relative_order(x)
0
```

```jldoctest
julia> gg = small_group(6, 1)
Pc group of order 6

julia> x = gg[1]^5*gg[2]^-4
f1*f2^2

julia> relative_order(x)
2
```
"""
relative_order(g::Union{PcGroupElem,SubPcGroupElem}) = relative_order(ZZRingElem, g)

function relative_order(
  ::Type{T}, g::Union{PcGroupElem,SubPcGroupElem}
) where {T<:IntegerUnion}
  # check if we have a PcpGroup element
  gObj = GapObj(g)
  if GAPWrap.IsPcpElement(gObj)
    return T(GAPWrap.RelativeOrder(gObj))
  else # finite PcGroup
    pcgs = GAPWrap.FamilyPcgs(GapObj(parent(g)))
    return T(GAPWrap.RelativeOrderOfPcElement(pcgs, gObj))
  end
end

"""
    depth(g::Union{PcGroupElem,SubPcGroupElem})

Return the depth of `g` as integer, relative to the defining generators.

# Examples

```jldoctest
julia> g = abelian_group(PcGroup, [0, 5])
Pc group of infinite order

julia> x = g[1]^-3 * g[2]^-3
g1^-3*g2^2

julia> depth(x)
1
```

```jldoctest
julia> gg = small_group(6, 1)
Pc group of order 6

julia> x = gg[1]^5*gg[2]^-4
f1*f2^2

julia> depth(x)
1
```
"""
function depth(g::Union{PcGroupElem,SubPcGroupElem})
  # check if we have a PcpGroup element
  gObj = GapObj(g)
  if GAPWrap.IsPcpElement(gObj)
    return GAPWrap.Depth(gObj)
  else # finite PcGroup
    return GAPWrap.DepthOfPcElement(GAPWrap.FamilyPcgs(GapObj(parent(g))), gObj)
  end
end

"""
    leading_exponent(::Type{T} = ZZRingElem, g::Union{PcGroupElem,SubPcGroupElem})

Return the leading exponent of `g` as an instance of the type `T`,
relative to the defining generators. Throws an error if `g` is the neutral element.

# Examples

```jldoctest
julia> g = abelian_group(PcGroup, [0, 5])
Pc group of infinite order

julia> x = g[1]^-3 * g[2]^-3
g1^-3*g2^2

julia> leading_exponent(x)
-3
```

```jldoctest
julia> gg = small_group(6, 1)
Pc group of order 6

julia> x = gg[1]^5*gg[2]^-4
f1*f2^2

julia> leading_exponent(x)
1
```
"""
leading_exponent(g::Union{PcGroupElem,SubPcGroupElem}) = leading_exponent(ZZRingElem, g)

function leading_exponent(
  ::Type{T}, g::Union{PcGroupElem,SubPcGroupElem}
) where {T<:IntegerUnion}
  # check if we have a PcpGroup element
  gObj = GapObj(g)
  exp = if GAPWrap.IsPcpElement(gObj)
    GAPWrap.LeadingExponent(gObj)
  else # finite PcGroup
    GAPWrap.LeadingExponentOfPcElement(GAPWrap.FamilyPcgs(GapObj(parent(g))), gObj)
  end

  # if GAP returns fail, error otherwise exp
  return exp == GAP.Globals.fail ? error("element has no leading exponent") : T(exp)
end

"""
    hirsch_length(G::PcGroup)

Return the Hirsch length of `G`.

# Examples

```jldoctest
julia> g = abelian_group(PcGroup, [0, 5])
Pc group of infinite order

julia> hirsch_length(g)
1
```

```jldoctest
julia> gg = small_group(6, 1)
Pc group of order 6

julia> hirsch_length(gg)
0
```
"""
function hirsch_length(G::PcGroup)
  GG = GapObj(G)
  if GAPWrap.IsPcpGroup(GG)
    return GAPWrap.HirschLength(GG)
  else # finite PcGroup
    return 0
  end
end

"""
    letters(g::Union{PcGroupElem, SubPcGroupElem})

Return the letters of `g` as a list of integers, each entry corresponding to
a group generator.

This method can produce letters represented by negative numbers. A negative number 
indicates the inverse of the generator at the corresponding positive index.

For example, as shown below, an output of `-1` refers to the "inverse of the first generator".

See also [`syllables(::Union{PcGroupElem, SubPcGroupElem})`](@ref).

# Examples

```jldoctest
julia> g = abelian_group(PcGroup, [0, 5])
Pc group of infinite order

julia> x = g[1]^-3 * g[2]^-3
g1^-3*g2^2

julia> letters(x)
5-element Vector{Int64}:
 -1
 -1
 -1
  2
  2
```

```jldoctest
julia> gg = small_group(6, 1)
Pc group of order 6

julia> x = gg[1]^5*gg[2]^-4
f1*f2^2

julia> letters(x)
3-element Vector{Int64}:
 1
 2
 2
```
"""
function letters(g::Union{PcGroupElem, SubPcGroupElem})
  # check if we have a PcpGroup element
  if GAPWrap.IsPcpElement(GapObj(g))
    exp = GAPWrap.Exponents(GapObj(g))

    # Should we check if the output is not larger than the
    # amount of generators? Requires use of `parent`.
    # @assert length(exp) == length(gens(parent(g)))

    w = [sign(e) * i for (i, e) in enumerate(exp) for _ in 1:abs(e)]
    return Vector{Int}(w)
  else # finite PcGroup
    w = GAPWrap.UnderlyingElement(GapObj(g))
    return Vector{Int}(GAPWrap.LetterRepAssocWord(w))
  end
end

"""
    syllables(g::Union{PcGroupElem, SubPcGroupElem})

Return the syllables of `g` as a list of pairs of integers, each entry corresponding to
a group generator and its exponent.

See also [`letters(::Union{PcGroupElem, SubPcGroupElem})`](@ref).

# Examples

```jldoctest
julia> gg = small_group(6, 1)
Pc group of order 6

julia> x = gg[1]^5*gg[2]^-4
f1*f2^2

julia> s = syllables(x)
2-element Vector{Pair{Int64, ZZRingElem}}:
 1 => 1
 2 => 2

julia> gg(s)
f1*f2^2

julia> gg(s) == x
true
```

```jldoctest
julia> g = abelian_group(PcGroup, [5, 0])
Pc group of infinite order

julia> x = g[1]^-3 * g[2]^-3
g1^2*g2^-3

julia> s = syllables(x)
2-element Vector{Pair{Int64, ZZRingElem}}:
 1 => 2
 2 => -3

julia> g(s)
g1^2*g2^-3

julia> g(s) == x
true
```
"""
function syllables(g::Union{PcGroupElem, SubPcGroupElem})
  # check if we have a PcpGroup element
  if GAPWrap.IsPcpElement(GapObj(g))
    l = GAPWrap.GenExpList(GapObj(g))
  else # finite PcGroup
    l = GAPWrap.ExtRepOfObj(GapObj(g))
  end

  @assert iseven(length(l))
  return Pair{Int, ZZRingElem}[l[i-1] => l[i] for i = 2:2:length(l)]
end

# Convert syllables in canonical form into exponent vector
function _exponent_vector(sylls::Vector{Pair{Int64, ZZRingElem}}, n)
  res = zeros(ZZRingElem, n)
  for pair in sylls
    @assert res[pair.first] == 0 #just to make sure 
    res[pair.first] = pair.second
  end
  return res
end

# Convert syllables in canonical form into group element
function (G::PcGroup)(sylls::Vector{Pair{Int64, ZZRingElem}}; check::Bool=true)
  # check if the syllables are in canonical form
  if check
    indices = map(p -> p.first, sylls)
    @req allunique(indices) "given syllables have repeating generators"
    @req issorted(indices) "given syllables must be in ascending order"
  end

  e = _exponent_vector(sylls, ngens(G))

  # check if G is an underlying PcpGroup
  GG = GapObj(G)
  if GAPWrap.IsPcpGroup(GG)
    coll = GAPWrap.Collector(GG)
    x = GAPWrap.PcpElementByExponentsNC(coll, GapObj(e, true))
  else # finite PcGroup
    pcgs = GAPWrap.FamilyPcgs(GG)
    x = GAPWrap.PcElementByExponentsNC(pcgs, GapObj(e, true))
  end
  
  return Oscar.group_element(G, x)
end

# Create an Oscar collector from a GAP collector.

const SCP_UNDERLYING_FAMILY = GAP.Globals.SCP_UNDERLYING_FAMILY             # = 1 - the family of our free grp elms
const SCP_RWS_GENERATORS = GAP.Globals.SCP_RWS_GENERATORS                   # = 2 - the free grp generators used
const SCP_NUMBER_RWS_GENERATORS = GAP.Globals.SCP_NUMBER_RWS_GENERATORS     # = 3 - number of generators
const SCP_DEFAULT_TYPE = GAP.Globals.SCP_DEFAULT_TYPE                       # = 4 - default type of the result
const SCP_IS_DEFAULT_TYPE = GAP.Globals.SCP_IS_DEFAULT_TYPE                 # = 5 - tester for default type
const SCP_RELATIVE_ORDERS = GAP.Globals.SCP_RELATIVE_ORDERS                 # = 6 - list of relative orders
const SCP_POWERS = GAP.Globals.SCP_POWERS                                   # = 7 - list of power rhs
const SCP_CONJUGATES = GAP.Globals.SCP_CONJUGATES                           # = 8 - list of list of conjugates rhs
const SCP_INVERSES = GAP.Globals.SCP_INVERSES                               # = 9 - list of inverses of the gens
const SCP_COLLECTOR = GAP.Globals.SCP_COLLECTOR                             # = 10 - collector to use
const SCP_AVECTOR = GAP.Globals.SCP_AVECTOR                                 # = 11 - avector

"""
    collector([::Type{T} = ZZRingElem, ]G::PcGroup) where T <: IntegerUnion

Return a collector object for `G`.

# Examples
```jldoctest
julia> g = small_group(12, 3)
Pc group of order 12

julia> c = collector(g);

julia> gc = pc_group(c)
Pc group of order 12

julia> is_isomorphic(g, gc)
true
```
"""
function collector(::Type{T}, G::PcGroup) where T <: IntegerUnion
  Fam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(GapObj(G)))
  GapC = GAP.getbangproperty(Fam, :rewritingSystem)::GapObj

  n = GAP.getbangindex(GapC, SCP_NUMBER_RWS_GENERATORS)::Int
  c = collector(n, T)

  c.relorders = Vector{T}(GAP.getbangindex(GapC, SCP_RELATIVE_ORDERS)::GapObj)

  Gap_powers = Vector{Union{Nothing, GapObj}}(GAP.getbangindex(GapC, SCP_POWERS)::GapObj, recursive = false)
  for i in 1:length(Gap_powers)
    if Gap_powers[i] !== nothing
      l = GAPWrap.ExtRepOfObj(Gap_powers[i])
      c.powers[i] = Pair{Int,T}[l[k-1] => T(l[k]) for k in 2:2:length(l)]
    end
  end

  Gap_conj = Vector{Vector{Union{Nothing, GapObj}}}(GAP.getbangindex(GapC, SCP_CONJUGATES)::GapObj)
  for i in 1:length(Gap_conj)
    for j in 1:length(Gap_conj[i])
      if Gap_conj[i][j] !== nothing
        l = GAPWrap.ExtRepOfObj(Gap_conj[i][j])
        c.conjugates[j,i] = Pair{Int,T}[l[k-1] => T(l[k]) for k in 2:2:length(l)]
      end
    end
  end

#TODO: deal also with the from-the-left-collector from the Polycyclic package

# c.X = GapC
# c.F = FPGroup(GAP.getbangproperty(GAP.getbangindex(GapC, 1)::GapObj, :freeGroup)::GapObj)
#TODO: Set these known data.
#      Currently this does not work because somehow `GroupByRws`
#      requires a *mutable* GAP collector, and `GapC` is immutable.
#      (Change `pc_group` to not call `GroupByRWS` in this case?
#      Forbid `set_power!` etc. in this case?)

  return c
end

collector(G::PcGroup) = collector(ZZRingElem, G)

# GAP wrappers for group encoding / decoding

"""
   encode(G::PcGroup)

Return a `ZZRingElem` representing the polycyclic group `G`,
using the same encoding as GAP's `CodePcGroup` and Magma's `SmallGroupEncoding`.
Currently only defined for `PcGroup`, not `SubPcGroup`.

# Examples
```jldoctest
julia> G = small_group(12, 2)
Pc group of order 12

julia> code = encode(G)
266

julia> H = pc_group(order(G), code)
Pc group of order 12

julia> encode(G) == encode(H)
true
```
"""
function encode(G::PcGroup)
  return ZZ(GAP.Globals.CodePcGroup(GapObj(G))::GapInt)
end

"""
   pc_group(order::IntegerUnion, code::IntegerUnion)

Given an integer `order` and an integer `code`, return the polycyclic group it encodes.
Both `order` and `code` can be of type `Int`, `BigInt`, or `ZZRingElem`.
The accepted codes and resulting groups match those of GAP's `PcGroupCode` and Magma's `SmallGroupDecoding`.

# Examples
```jldoctest
julia> G = small_group(12, 2)
Pc group of order 12

julia> code = encode(G)
266

julia> H = pc_group(order(G), code)
Pc group of order 12

julia> encode(G) == encode(H)
true
```
"""
function pc_group(order::IntegerUnion, code::IntegerUnion)
  return PcGroup(GAP.Globals.PcGroupCode(GAP.GapInt(BigInt(code)),GAP.GapInt(BigInt(order))))
end
