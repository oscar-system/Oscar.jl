const AdmissibleSingularQuoRingElem = Union{<:MPolyQuoRingElem{<:MPolyRingElem{<:FieldElem}}}

@attr Any function kernel(
    f::FreeModuleHom{DomainType, CodomainType}
  ) where {
           T<:AdmissibleSingularQuoRingElem,
           DomainType<:FreeMod{T},
           CodomainType<:FreeMod{T}
          }
  is_graded(f) && return _graded_kernel(f)
  return _simple_kernel(f)
end

function in(v::FreeModElem{T}, I::SubModuleOfFreeModule{T}) where {T <: AdmissibleSingularQuoRingElem}
  F = ambient_free_module(I)
  @assert parent(v) === F
  return iszero(reduce(v, standard_basis(I, ordering=_default_ordering(F))))
end

### missing functions which needed to be overwritten
@attr ModuleOrdering function _default_ordering(F::FreeMod{<:AdmissibleSingularQuoRingElem})
  R = base_ring(F)
  # TODO: This ignores the `default_ordering` on the quotient ring. 
  # The problem is that the latter returns a monomial ordering on the 
  # polynomial `base_ring` and that is incompatible with the ordering 
  # on the generators of `F`. The workaround here is to invent some 
  # dummy ordering on the quotient ring (which maybe does not even make 
  # sense!), just to satisfy the "generic" code for `ModuleGens`. 
  return lex(gens(F))*Orderings.MonomialOrdering(R, Orderings.SymbOrdering(:degrevlex, collect(1:ngens(R))))
end

function singular_module(
    F::FreeMod{<:AdmissibleSingularQuoRingElem}, 
    ordering::ModuleOrdering
  )
  Sx = singular_poly_ring(base_ring(F), induced_ring_ordering(ordering))
  return Singular.FreeModule(Sx, dim(F))
end

# Overwritten only to call `_default_ordering` instead of `default_ordering` 
# for the reasons explained above.
function _simple_kernel(h::FreeModuleHom{<:FreeMod{T}, <:FreeMod{T}}) where {T<:MPolyQuoRingElem}
  F = domain(h)
  G = codomain(h)
  g = images_of_generators(h)
  b = ModuleGens(g, G, _default_ordering(G))
  M = syzygy_module(b)
  v = elem_type(F)[F(coordinates(repres(w))) for w in gens(M) if !is_zero(w)]
  return sub(F, v)
end

# Overwritten because the original method has the restriction to be for 
# polynomial rings only. 
function syzygy_module(F::ModuleGens{T}; sub = FreeMod(base_ring(F.F), length(oscar_generators(F)))) where {T <: AdmissibleSingularQuoRingElem}
  singular_assure(F)
  # TODO Obtain the Gröbner basis and cache it
  s = Singular.syz(singular_generators(F))
  return SubquoModule(sub, s)
end

# Called somewhere in the internals
function default_ordering(I::SubModuleOfFreeModule{T}) where {T <: AdmissibleSingularQuoRingElem}
  return _default_ordering(ambient_free_module(I))
end

# Overwritten only because original signature was for polynomial rings only.
function standard_basis(F::ModuleGens{T}, reduced::Bool=false) where {T <: AdmissibleSingularQuoRingElem}
  R = base_ring(F)
  @req is_exact_type(elem_type(coefficient_ring(R))) "This functionality is only supported over exact fields."
  singular_assure(F)
  if reduced
    @assert Singular.has_global_ordering(base_ring(F.SF))
  end
  if singular_generators(F).isGB && !reduced
    return F
  end
  return ModuleGens(F.F, Singular.std(singular_generators(F), complete_reduction=reduced))
end

# A hacky overwrite to just make things run.
function Orderings.is_global(ord::MonomialOrdering{<:MPolyQuoRing})
  return true # At the moment, it is not possible to define non-global orderings via user-facing methods.
end

# Overwritten only because the original signature was exclusively 
# for polynomial rings.
function normal_form(M::ModuleGens{T}, GB::ModuleGens{T}) where {T <: AdmissibleSingularQuoRingElem}
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

function lift_std(M::ModuleGens{T}) where {T <: AdmissibleSingularQuoRingElem}
  singular_assure(M)
  R = base_ring(M)
  G,Trans_mat = Singular.lift_std(singular_generators(M)) # When Singular supports reduction add it also here
  mg = ModuleGens(M.F, G)
  mg.isGB = true
  mg.S.isGB = true
  mg.ordering = _default_ordering(M.F)
  mat = map_entries(R, transpose(Trans_mat))
  set_attribute!(mg, :transformation_matrix => mat)
  return mg, mat
end

function lift_std(M::ModuleGens{T}, ordering::ModuleOrdering) where {T <: AdmissibleSingularQuoRingElem}
  M = ModuleGens(M.O, M.F, ordering)
  mg, mat = lift_std(M)
  mg.ordering = ordering
  return mg, mat
end

# TODO: This is a hotfix which does not adapt automaticall 
# with changes to `AdmissibleSingularQuoRingElem`.
function sparse_row(
    A::MPolyQuoRing{<:MPolyRingElem{<:FieldElem}},
    svec::Singular.svector, rng::AbstractUnitRange
  )
  pre_res = sparse_row(base_ring(A), svec, rng)
  return map_entries(A, pre_res)
end

# Overwritten because the caching and comparison of monomial orderings 
# is not functional as of yet. Since there is only one monomial ordering 
# available at the moment, we disabled this overhead completely.
function standard_basis(
    submod::SubModuleOfFreeModule{T}; 
    ordering::Union{ModuleOrdering, Nothing} = nothing
  ) where {T <: AdmissibleSingularQuoRingElem}
  if !isempty(submod.groebner_basis)
    for (ord, gb) in submod.groebner_basis
      return gb
    end
  end
    
  @req is_exact_type(elem_type(base_ring(submod))) "This functionality is only supported over exact fields."
  gb = get!(submod.groebner_basis, ordering) do
    compute_standard_basis(submod, ordering)
  end::ModuleGens{T}
  return gb
end

# Should eventually go elsewhere; needed for hashing of 
# `ModuleOrderings` for lookups of standard bases in dictionaries.
number_of_variables(A::MPolyQuoRing) = number_of_variables(base_ring(A))

