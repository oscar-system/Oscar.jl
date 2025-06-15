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

Return the generators of `M` from the Singular side.
"""
function singular_generators(M::ModuleGens)
  singular_assure(M)
  return M.S
end

@doc raw"""
    singular_ordering(M::ModuleGens)

Return the ordering of `M` from the Singular side.
"""
function singular_ordering(M::ModuleGens)
    return Singular.ordering(base_ring(singular_freemodule(M)))
end

@doc raw"""
    singular_freemodule(M::ModuleGens)

Return the ambient free module of `M` from the Singular side.
"""
function singular_freemodule(M::ModuleGens)
    singular_assure(M)
    return M.SF
end


@doc raw"""
    has_global_singular_ordering(M::ModuleGens)

Return whether the ordering of `M` from the Singular side is global.
"""
function has_global_singular_ordering(M::ModuleGens)
    return Singular.has_global_ordering(base_ring(singular_freemodule(M)))
end

@doc raw"""
    oscar_generators(M::ModuleGens)  

Return the generators of `M` from the Oscar side.

If fields of `M` from the Oscar side are not defined, they
are computed, given the Singular side.
"""
function oscar_generators(M::ModuleGens)
  if !isdefined(M, :O)
    SI = singular_generators(M)
    if iszero(SI)
      M.O = elem_type(M.F)[zero(M.F) for _ in 1:ngens(SI)]
    else
      M.O = [M.F(SI[i]) for i=1:Singular.ngens(SI)]
    end
  end
  return M.O
end

@doc raw"""
    oscar_free_module(M::ModuleGens)  

Return the ambient free module of `M` on the Oscar side.
"""
function oscar_free_module(M::ModuleGens)
    return M.F
end

@doc raw"""
    iszero(M::ModuleGens)

Check if `M` is zero.
"""
function iszero(M::ModuleGens)
  return all(iszero, oscar_generators(M))
end

function show(io::IO, F::ModuleGens)
  if F.isGB
    show_gb(io, F)
  else
    print(io, "Module generating system of length ", length(F))
    for i=1:length(F)
      if isassigned(F.O, i)
        print(io, "\n", i, ": ", F.O[i])
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
end

# i-th entry of module generating set (taken from Oscar side)
getindex(F::ModuleGens, i::Int) = oscar_generators(F)[i]

@doc raw"""
    union(M::ModuleGens, N::ModuleGens)

Compute the union of `M` and `N`.
"""
function union(M::ModuleGens, N::ModuleGens)
  @assert oscar_free_module(M) === oscar_free_module(M)
  O = vcat(oscar_generators(M), oscar_generators(N))
  return ModuleGens(oscar_free_module(M), O)
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
function singular_module(F::FreeMod{<:MPolyRingElem}, ordering::ModuleOrdering)
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
function (F::FreeMod{<:MPolyRingElem})(s::Singular.svector)
  Rx = base_ring(F)
  row = _build_sparse_row(Rx, s)
  return FreeModElem(row, F)
end

function (F::FreeMod{<:MPolyQuoRingElem})(s::Singular.svector)
  Qx = base_ring(F)::MPolyQuoRing
  Rx = base_ring(Qx)::MPolyRing
  row = _build_sparse_row(Rx, s; cast=Qx)
  return FreeModElem(row, F)
end

function _build_sparse_row(
    Rx::MPolyRing, s::Singular.svector; 
    cast::Ring=Rx
  )
  is_zero(length(s)) && return sparse_row(cast)
  R = coefficient_ring(Rx)

  # shortcuts in order not to allocate the dictionary
  ctx = MPolyBuildCtx(Rx)
  if isone(length(s))
    (i, e, c) = first(s)
    push_term!(ctx, R(c), e)
    cast !== Rx && return sparse_row(cast, [(i, cast(finish(ctx)))])
    return sparse_row(Rx, [(i, finish(ctx))])
  end

  cache = IdDict{Int, typeof(ctx)}()
  last_index = 0
  for (i, e, c) in s
    if i == last_index
      push_term!(ctx, R(c), e)
      continue
    end
    last_index = i
    ctx = get!(cache, i) do
      MPolyBuildCtx(Rx)
    end
    push_term!(ctx, R(c), e)
  end
  cast !== Rx && return sparse_row(cast, [(i, cast(finish(ctx))) for (i, ctx) in cache])
  return sparse_row(Rx, [(i, finish(ctx)) for (i, ctx) in cache])
end

# After creating the required infrastructure in Singular,
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
  S = singular_generators(generators)
  b = ModuleGens([a], singular_freemodule(generators))
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
function coordinates(a::FreeModElem{T}, generators::ModuleGens{T}) where {T<:MPolyRingElem}
    if !has_global_singular_ordering(generators)
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

  S = singular_generators(generators)
  S.isGB = generators.isGB
  b = ModuleGens([a], singular_freemodule(generators))
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
  @assert oscar_free_module(M) === oscar_free_module(GB)
  @assert GB.isGB # TODO When Singular.jl can handle reduce with non-GB remove this

  P = isdefined(GB, :quo_GB) ? union(GB, GB.quo_GB) : GB

  red = _reduce(singular_generators(M), singular_generators(P))
  res = ModuleGens(oscar_free_module(M), red)
  return res
end

@doc raw"""
    normal_form_with_unit(M::ModuleGens, GB::ModuleGens)

Compute a normal form of `M` (that is of each element of `M`) with respect to the Gröbner basis `GB`.
Moreover, return a vector `U` of unit elements such that 
`U[i]*M[i]` is the `i`th element of the normal form `ModuleGens`.
"""
function normal_form_with_unit(M::ModuleGens{T}, GB::ModuleGens{T}) where {T <: MPolyRingElem}
  @assert oscar_free_module(M) === oscar_free_module(GB)
  @assert GB.isGB # TODO When Singular.jl can handle reduce/nf with non-GB remove this
  if !is_global(GB.ordering)
    error("normal_form_with_unit not yet implemented for non-global orderings") # This function doesn't exist yet in Singular.jl
  end
  R = base_ring(M)

  P = isdefined(GB, :quo_GB) ? union(GB, GB.quo_GB) : GB

  red = _reduce(singular_generators(M), singular_generators(P))
  res = ModuleGens(oscar_free_module(M), red)
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
  nf = normal_form(ModuleGens([v], parent(v)), GB)
  return oscar_generators(nf)[1]
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

