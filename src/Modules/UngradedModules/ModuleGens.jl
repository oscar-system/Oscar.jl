###############################################################################
# ModuleGens constructors
###############################################################################

@doc raw"""
    ModuleGens(O::Vector{<:FreeModElem}, F::FreeMod{T}, SF::Singular.FreeMod) where T

Construct `ModuleGens` from an array of free module elements, specifying the free module 
and Singular free module. 
This function is only useful indirectly.
"""
ModuleGens(O::Vector{<:FreeModElem}, F::FreeMod{T}, SF::Singular.FreeMod) where T = ModuleGens{T}(O, F, SF)

@doc raw"""
    ModuleGens(F::FreeMod{S}, s::Singular.smodule) where {S}

Construct `ModuleGens` from a given Singular submodule.
"""
ModuleGens(F::FreeMod{S}, s::Singular.smodule) where {S} = ModuleGens{S}(F, s)

@doc raw"""
    ModuleGens(O::Vector{<:FreeModElem})

Construct `ModuleGens` from an array of free module elements.

!!! note 
    
    The array must not be empty.
"""
function ModuleGens(O::Vector{<:FreeModElem})
  # TODO Empty generating set
  @assert length(O) > 0
  #SF = singular_module(parent(O[1]))
  #return ModuleGens(O, SF)
  return ModuleGens(O, parent(O[1]))
end

@doc raw"""
    ModuleGens(O::Vector{<:FreeModElem}, F::FreeMod{T}) where {T}

Construct `ModuleGens` from an array of free module elements, specifying the free module.

!!! note

    The array might be empty.
"""
function ModuleGens(O::Vector{<:FreeModElem}, F::FreeMod{T}) where {T}
  #SF = singular_module(F)
  return ModuleGens{T}(O, F)
end

@doc raw"""
    ModuleGens(O::Vector{<:FreeModElem}, F::FreeMod{T}, ordering::ModuleOrdering) where {T}

Construct `ModuleGens` from an array of free module elements, specifying the free module.
Moreover, the ordering is defined by `ordering`. 

!!! note

    The array might be empty.
"""
function ModuleGens(O::Vector{<:FreeModElem}, F::FreeMod{T}, ordering::ModuleOrdering) where {T}
  SF = singular_module(F, ordering)
  M = ModuleGens{T}(O, F, SF)
  M.ordering = ordering
  return M
end

@doc raw"""
    ModuleGens(O::Vector{<:FreeModElem}, SF::Singular.FreeMod)

Construct `ModuleGens` from an array of free module elements, specifying the Singular free module.

!!! note 

    The array might be empty.
"""
function ModuleGens(O::Vector{<:FreeModElem}, SF::Singular.FreeMod)
  return ModuleGens{elem_type(base_ring(parent(O[1])))}(O, parent(O[1]), SF)
end

@doc raw"""
    base_ring(M::ModuleGens)

Return the base ring of `M` (that is, if `M` is an `R`-module, return `R`).
"""
function base_ring(M::ModuleGens)
  return base_ring(M.F)
end

base_ring_type(::Type{ModuleGens{T}}) where {T} = base_ring_type(FreeMod{T})

@doc raw"""
    singular_generators(M::ModuleGens)

Return the generators of `M` from Singular side.
"""
function singular_generators(M::ModuleGens)
  singular_assure(M)
  return M.S
end

@doc raw"""
    oscar_generators(M::ModuleGens)  

Return the generators of `M` from the Oscar side.
"""
function oscar_generators(M::ModuleGens)
  oscar_assure(M)
  return M.O
end

@doc raw"""
    iszero(M::ModuleGens)

Check if `M` is zero.
"""
function iszero(M::ModuleGens)
  oscar_assure(M)
  return all(iszero, M.O)
end

function show(io::IO, F::ModuleGens)
  if F.isGB
    show_gb(io, F)
  else
    print(io, "Module generating system of length ", length(F))
    for i=1:length(F)
      if isassigned(F.O, i)
        print(io, "\n", i, " -> ", F.O[i])
      end
    end
  end
end

@doc raw"""
    length(F::ModuleGens)

Return the number of elements of the module generating set.

!!! note

    This is not the length in the mathematical sense!
"""
length(F::ModuleGens) = length(oscar_generators(F))

@doc raw"""
    number_of_generators(F::ModuleGens)

Return the number of elements of the module generating set.
"""
number_of_generators(F::ModuleGens) = length(oscar_generators(F))

# i-th entry of module generating set on Oscar side
# Todo: clean up, convert or assure
function getindex(F::ModuleGens, ::Val{:O}, i::Int)
  if !isassigned(F.O, i)
    F.O[i] = F.F(singular_generators(F)[i])
  end
  return oscar_generators(F)[i]
end

# i-th entry of module generating set on Singular side
# Todo: clean up, convert or assure
function getindex(F::ModuleGens, ::Val{:S}, i::Int)
  singular_assure(F)
  if !isdefined(F, :S)
    F.S = Singular.Module(base_ring(F.SF), [F.SF(x) for x = oscar_generators(F)]...)
  end
  return F.S[i]
end

@doc raw"""
    oscar_assure(F::ModuleGens)

If fields of `F` from the Oscar side are not defined, they
are computed, given the Singular side.
"""
function oscar_assure(F::ModuleGens)
  if !isdefined(F, :O)
    F.O = [F.F(singular_generators(F)[i]) for i=1:Singular.ngens(singular_generators(F))]
  end
end

@doc raw"""
    singular_assure(F::ModuleGens)

If fields of `F` from the Singular side are not defined, they
are computed, given the Oscar side.
"""
function singular_assure(F::ModuleGens)
  if !isdefined(F, :S) || !isdefined(F, :SF)
    if isdefined(F, :ordering)
      SF = singular_module(F.F, F.ordering)
    else
      SF = singular_module(F.F)
    end
    sr = base_ring(SF)
    F.SF = SF
    if length(F) == 0
      F.S = Singular.Module(sr, Singular.vector(sr, sr(0)))
      return 
    end
    F.S = Singular.Module(base_ring(F.SF), [F.SF(x) for x = oscar_generators(F)]...)
    return
  end
  #F[Val(:S), 1]
end

# i-th entry of module generating set (taken from Oscar side)
getindex(F::ModuleGens, i::Int) = getindex(F, Val(:O), i)

@doc raw"""
    union(M::ModuleGens, N::ModuleGens)

Compute the union of `M` and `N`.
"""
function union(M::ModuleGens, N::ModuleGens)
  @assert M.F === N.F
  O = vcat(M.O, N.O)
  return ModuleGens(M.F, O)
end

@doc raw"""
    singular_module(F::FreeMod)

Create a Singular module from a given free module.
"""
function singular_module(F::FreeMod)
  Sx = singular_poly_ring(base_ring(F), keep_ordering=false)
  return Singular.FreeModule(Sx, dim(F))
end

@doc raw"""
    singular_module(F::FreeMod, ordering::ModuleOrdering)

Create a Singular module from a given free module over the given Singular polynomial ring.
"""
function singular_module(F::FreeMod, ordering::ModuleOrdering)
  Sx = singular_poly_ring(base_ring(F), singular(ordering))
  return Singular.FreeModule(Sx, dim(F))
end

@doc raw"""
    (SF::Singular.FreeMod)(m::FreeModElem)

Convert a free module element to the Singular side.
"""
function (SF::Singular.FreeMod)(m::FreeModElem)
  g = Singular.gens(SF)
  e = SF()
  Sx = base_ring(SF)
  for (p,v) in coordinates(m)
    e += Sx(v)*g[p]
  end
  return e
end

@doc raw"""
    (F::FreeMod)(s::Singular.svector)

Convert a Singular vector to a free module element.
"""
function (F::FreeMod)(s::Singular.svector)
  pos = Int[]
  values = []
  Rx = base_ring(F)
  R = base_ring(Rx)
  for (i, e, c) = s
    f = Base.findfirst(x->x==i, pos)
    if f === nothing
      push!(values, MPolyBuildCtx(base_ring(F)))
      f = length(values)
      push!(pos, i)
    end
    push_term!(values[f], R(c), e)
  end
  pv = Tuple{Int, elem_type(Rx)}[(pos[i], base_ring(F)(finish(values[i]))) for i=1:length(pos)]
  return FreeModElem(sparse_row(base_ring(F), pv), F)
end

# After creating the required infrastruture in Singular,
# to facilitate the double book-keeping, the signature
# lift(G1::ModuleGens{T}, G2::ModuleGens{T}) should go to Singular
# and lift(a::FreeModElem{T}, generators::ModuleGens{T}) call it

function lift(G1::ModuleGens{T}, G2::ModuleGens{T}) where {T <: MPolyRingElem}
  results = Vector{SRow{T, Vector{T}}}()
  for i in 1:ngens(G1)
      s_row = lift(G1[i], G2)
      push!(results, s_row)
  end
  return results
end

function lift(a::FreeModElem{T}, generators::ModuleGens{T}) where {T <: MPolyRingElem}
  if iszero(a)
    return sparse_row(base_ring(parent(a)))
  end
  singular_assure(generators)
  S = singular_generators(generators)
  b = ModuleGens([a], generators.SF)
  singular_assure(b)
  s, r = Singular.lift(S, singular_generators(b))
  if Singular.ngens(s) == 0 || iszero(s[1])
    error("The free module element is not liftable to the given generating system.")
  end
  Rx = base_ring(generators)
  return sparse_row(Rx, s[1], 1:ngens(generators))
end

@doc raw"""
    coordinates(a::FreeModElem{T}, generators::ModuleGens{T})

Compute a sparse row `r` such that `a = sum([r[i]*gen(generators,i) for i in 1:ngens(generators)])`.
If no such `r` exists, an exception is thrown.
"""
function coordinates(a::FreeModElem{T}, generators::ModuleGens{T}) where {T <: MPolyRingElem}
  singular_assure(generators)
  if !Singular.has_global_ordering(base_ring(generators.SF))
    error("Ordering must be global")
  end
  return lift(a, generators)
end

@doc raw"""
    coordinates_via_transform(a::FreeModElem{T}, generators::ModuleGens{T}) where T

Let `generators` be a Gröbner basis and let `A*M = generators`.
Compute a sparse row `r` such that `a = sum([r[i]*gen(M,i) for i in 1:ngens(M)])`.
If no such `r` exists, an exception is thrown.
"""
function coordinates_via_transform(a::FreeModElem{T}, generators::ModuleGens{T}) where T
  A = get_attribute(generators, :transformation_matrix)
  A === nothing && error("No transformation matrix in the Gröbner basis.")
  if iszero(a)
    return sparse_row(base_ring(parent(a)))
  end
  @assert generators.isGB
  if base_ring(generators) isa Union{MPolyQuoRing,MPolyRing}
    if !is_global(generators.ordering)
      error("Ordering is not global")
    end
  end

  singular_assure(generators)
  S = singular_generators(generators)
  S.isGB = generators.isGB
  b = ModuleGens([a], generators.SF)
  singular_assure(b)
  s, r = Singular.lift(S, singular_generators(b)) # Possibly use division with remainder
  if Singular.ngens(s) == 0 || iszero(s[1])
    error("The free module element is not liftable to the given generating system.")
  end
  Rx = base_ring(generators)
  coords_wrt_groebner_basis = sparse_row(Rx, s[1], 1:ngens(generators))

  return coords_wrt_groebner_basis * sparse_matrix(A)
end

@doc raw"""
    coordinates(a::FreeModElem, M::SubModuleOfFreeModule, task::Symbol = :auto)

Compute a sparse row `r` such that `a = sum([r[i]*gen(M,i) for i in 1:ngens(M)])`.
If no such `r` exists, an exception is thrown.
For `task` there are the following options:
- `:auto` (default option): Use `:via_transform` if coefficient ring of base ring is a field, 
else use `:via_lift` 
- `:via_transform`: Compute first a Gröbner basis with a transformation matrix representing
the Gröbner basis in terms of the generators and cache the data
- `:via_lift`: Compute the lift and do not cache the auxiliary data
Note: `:via_lift` is typically faster than `:via_transform` for a single vector while the latter
is faster if many vectors are lifted
"""
function coordinates(a::FreeModElem, M::SubModuleOfFreeModule, task::Symbol = :auto)
  if iszero(a)
    return sparse_row(base_ring(parent(a)))
  end
  if task == :auto
    if coefficient_ring(base_ring(parent(a))) isa Field #base_ring(base_ring(...)) does not work for MPolyQuos
      task = :via_transform
    else
      task = :via_lift
    end
  end
  for i in 1:ngens(M)
    g = gen(M,i)
    if a == g
      R = base_ring(M)
      return sparse_row(R, [(i,R(1))])
    end
  end
  if task == :via_transform
    std, _ = lift_std(M)
    return coordinates_via_transform(a, std)
  elseif task == :via_lift
    return coordinates(a, M.gens)
  else
    error("Invalid task given.")
  end
end

@doc raw"""
    normal_form(M::ModuleGens, GB::ModuleGens)

Compute a normal form of `M` (that is of each element of `M`) with respect to the Gröbner basis `GB`.
"""
function normal_form(M::ModuleGens{T}, GB::ModuleGens{T}) where {T <: MPolyRingElem}
  @assert M.F === GB.F
  @assert GB.isGB # TODO When Singular.jl can handle reduce with non-GB remove this

  P = isdefined(GB, :quo_GB) ? union(GB, GB.quo_GB) : GB

  singular_assure(P)
  singular_assure(M)

  red = _reduce(M.S, P.S)
  res = ModuleGens(M.F, red)
  oscar_assure(res)
  return res
end

@doc raw"""
    normal_form_with_unit(M::ModuleGens, GB::ModuleGens)

Compute a normal form of `M` (that is of each element of `M`) with respect to the Gröbner basis `GB`.
Moreover, return a vector `U` of unit elements such that 
`U[i]*M[i]` is the `i`th element of the normal form `ModuleGens`.
"""
function normal_form_with_unit(M::ModuleGens{T}, GB::ModuleGens{T}) where {T <: MPolyRingElem}
  @assert M.F === GB.F
  @assert GB.isGB # TODO When Singular.jl can handle reduce/nf with non-GB remove this
  if !is_global(GB.ordering)
    error("normal_form_with_unit not yet implemented for non-global orderings") # This function doesn't exist yet in Singular.jl
  end
  R = base_ring(M)

  P = isdefined(GB, :quo_GB) ? union(GB, GB.quo_GB) : GB

  singular_assure(P)
  singular_assure(M)

  red = _reduce(M.S, P.S)
  res = ModuleGens(M.F, red)
  oscar_assure(res)
  return res, [R(1) for _ in 1:ngens(M)]
end

function normal_form_with_unit_and_coefficients(M::ModuleGens{T}, GB::ModuleGens{T}) where {T <: MPolyRingElem}
  # TODO requires additional functionality in Singular.jl
  error("Not yet implemented")
end

@doc raw"""
    normal_form(v::AbstractFreeModElem, GB::ModuleGens)

Compute a normal_form of `v` with respect to the Gröbner basis `GB`.
"""
function normal_form(v::AbstractFreeModElem, GB::ModuleGens)
  @assert GB.isGB
  return normal_form(ModuleGens([v], parent(v)), GB).O[1]
end

@doc raw"""
    normal_form_with_unit(v::AbstractFreeModElem, GB::ModuleGens)

Compute a normal form of `v` with respect to the Gröbner basis `GB`.
Moreover, return a unit `u` such that `u*v` is the normal form.
"""
function normal_form_with_unit(v::AbstractFreeModElem, GB::ModuleGens)
  @assert GB.isGB
  if !is_global(GB.ordering)
    error("normal_form_with_unit not yet implemented for non-global orderings") # This function doesn't exist yet in Singular.jl
  end
  red, unit = normal_form_with_unit(ModuleGens([v], parent(v)), GB)
  return red[1], unit[1]
end

@doc raw"""
    reduce(M::ModuleGens, GB::ModuleGens)

Reduce `M` with respect to the Gröbner basis `GB`.
"""
function reduce(M::ModuleGens{T}, GB::ModuleGens{T}) where {T <: MPolyRingElem}
  @assert GB.isGB
  @assert is_global(GB.ordering)
  return normal_form(M, GB)
end

