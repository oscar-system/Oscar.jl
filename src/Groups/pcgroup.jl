export collector,
       pc_group,
       set_commutator!,
       set_conjugate!,
       set_power!,
       set_relative_order!,
       set_relative_orders!

############################################################################

# Create GAP filters that describe
# - the union of `IsMultiplicativeElementWithInverseByPolycyclicCollector`
#   and `IsPcpElement` and
# - the union of `IsPcGroup` and `IsPcpGroup`.
# This must be done at runtime.
function __init_PcGroups()
  GAP.evalstr("""
      DeclareFilter("IsPcElementOrPcpElement");
      InstallTrueMethod(IsPcElementOrPcpElement, IsMultiplicativeElementWithInverseByPolycyclicCollector);
      InstallTrueMethod(IsPcElementOrPcpElement, IsPcpElement);
      BindGlobal("IsPcGroupOrPcpGroup", IsGroup and CategoryCollections(IsPcElementOrPcpElement));
      """)
end

############################################################################

# Create an Oscar collector object, its type parameter `T` describes
# the type of integers that occur as exponents.
#
# The idea is
# - to formulate the presentation in Julia,
#   and to create the corresponding GAP object on demand
#   when it gets accessed,
# - to create a GAP collector in `GAP.Globals.IsSingleCollectorRep`
#   if all relative orders are finite,
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
    collector(n::Int, ::Type{T} = fmpz) where T <: IntegerUnion

Return an empty collector object for a pc group with `n` generators
and exponents of type `T`.
It describes a free abelian group of rank `n`.

# Examples
```jldoctest
julia> G = pc_group(collector(2))
Pcp-group with orders [ 0, 0 ]

julia> is_abelian(G)
true
```
"""
collector(n::Int, ::Type{T} = fmpz) where T <: IntegerUnion = GAP_Collector{T}(n)


# Provide functions for entering data into the collector.

# utility:
# Convert a vector of pairs to a generator-exponent vector in GAP.
function _GAP_generator_exponent_vector(l::Vector{Pair{Int, T}}) where T <: IntegerUnion
  res = T[]
  for p in l
    push!(res, p.first)
    push!(res, p.second)
  end
  return GapObj(res, recursive = true)
end

"""
    set_relative_order!(c::Collector{T}, i::Int, relord::T) where T <: IntegerUnion

Set the relative order of the `i`-th generator of `c` to `relord`,
which must be either `0` (meaning infinite order) or a positive integer.

# Examples
```jldoctest
julia> c = collector(2, Int);

julia> set_relative_order!(c, 1, 2)
```
"""
function set_relative_order!(c::Collector{T}, i::Int, relord::T) where T <: IntegerUnion
  (0 < i && i <= c.ngens) || throw(ArgumentError("the collector has only $(c.ngens) generators not $i"))
  c.relorders[i] = relord

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
    set_relative_orders!(c::Collector{T}, relords::Vector{T})

Set all relative orders of the generators of `c`,
where the length of `relords` must be equal to the number of generators of
`c`, and `relords[i]` denotes the relative order of the `i`-th generator.
which must be either `0` (meaning infinite order) or a positive integer..

# Examples
```jldoctest
julia> c = collector(2);

julia> set_relative_orders!(c, fmpz[2, 0])
```
"""
function set_relative_orders!(c::Collector{T}, relords::Vector{T}) where T <: IntegerUnion
  length(relords) == c.ngens || throw(ArgumentError("the collector has $(c.ngens) generators not $(length(relords))"))
  c.relorders = copy(relords)

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
  (0 < i && i <= c.ngens) || throw(ArgumentError("the collector has only $(c.ngens) generators not $i"))
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
  (0 < i && i <= c.ngens) || throw(ArgumentError("the collector has only $(c.ngens) generators not $i"))
  i < j || throw(ArgumentError("only for i < j, but i = $i, j = $j"))
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
  (0 < i && i <= c.ngens) || throw(ArgumentError("the collector has only $(c.ngens) generators not $i"))
  i < j || throw(ArgumentError("only for i < j, but i = $i, j = $j"))
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
  G = free_group(c.ngens; eltype = :syllable)
  cGAP = GAP.Globals.SingleCollector(G.X, GapObj(c.relorders, recursive = true))::GapObj
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
<pc group of size 6 with 2 generators>

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

