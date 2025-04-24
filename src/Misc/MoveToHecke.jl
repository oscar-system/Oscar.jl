###############################################################################
# A place to accumulate code that should eventually be moved to Hecke.jl
###############################################################################

function getindex(r::Hecke.SRow, u::AbstractUnitRange)
  s = sparse_row(base_ring(r))
  shift = 1-first(u)
  for (p,v) = r
    if p in u
      push!(s.pos, p+shift)
      push!(s.values, v)
    end
  end
  return s
end

Oscar.canonical_unit(x::AbsSimpleNumFieldOrderQuoRingElem) = one(parent(x))

###############################################################################
# Part of https://github.com/thofma/Hecke.jl/pull/1800, but breaking for OSCAR 1.3.1.
# To allow progress, we include it here until a future breaking Hecke release.
#
function (::Type{T})(G::FinGenAbGroup) where T <: Group
  return codomain(isomorphism(T, G))
end

function (::Type{FinGenAbGroup})(G::Group)
  return codomain(isomorphism(FinGenAbGroup, G))
end

function isomorphism(::Type{FinGenAbGroup}, G::FinGenAbGroup; on_gens::Bool=false)
  # Known isomorphisms are cached in the attribute `:isomorphisms`.
  on_gens = true # we ignore the on_gens flag, the identity will *always* map gens onto gens
  isos = get_attribute!(Dict{Tuple{Type, Bool}, Any}, G, :isomorphisms)::Dict{Tuple{Type, Bool}, Any}
  return get!(isos, (FinGenAbGroup, on_gens)) do
    return id_hom(G)
  end::FinGenAbGroupHom
end

function isomorphism(::Type{T}, G::FinGenAbGroup; on_gens::Bool=false) where T <: Group
  throw(NotImplementedError(:isomorphism, T, G))
end

function isomorphism(::Type{FinGenAbGroup}, G::Group; on_gens::Bool=false)
  throw(NotImplementedError(:isomorphism, FinGenAbGroup, G))
end
#
###############################################################################

########################################################################
# Part of PR #4706
function is_equal_as_morphism(f::MapFromFunc, g::MapFromFunc)
  f === g && return true
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  f.f === g.f && return true
  error("comparison of maps $f and $g not possible")
end

function is_equal_as_morphism(f::Map, g::Map)
  f === g && return true
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  error("comparison of maps $f and $g not possible")
end
# end of changes in PR #4706
########################################################################
