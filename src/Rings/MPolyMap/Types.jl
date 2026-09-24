### Types for maps from polynomial rings
const _DomainTypes = Union{MPolyRing, MPolyQuoRing}

# Property tags for homomorphisms

abstract type MPolyHomKind end

struct HomUngraded           <: MPolyHomKind end
struct HomGraded             <: MPolyHomKind end
struct HomGradedToUngraded   <: MPolyHomKind end
struct HomUngradedToGraded   <: MPolyHomKind end

@attributes mutable struct MPolyAnyMap{
    D <: _DomainTypes,
    C <: NCRing,
    U,
    V,
    K <: MPolyHomKind} <: Map{D, C, Map, MPolyAnyMap}

  domain::D
  codomain::C
  coeff_map::U
  img_gens::Vector{V}
  temp_ring
  variable_indices::Vector{Int}

  function MPolyAnyMap{D, C, U, V, K}(domain::D,
                                      codomain::C,
                                      coeff_map::U,
                                      img_gens::Vector{V};
                                      check_for_mapping_of_vars::Bool=true
                                      ) where {D, C, U, V, K}
    @assert V === elem_type(C)
    for g in img_gens
      @assert parent(g) === codomain "elements does not have the correct parent"
    end
    result = new{D, C, U, V, K}(domain, codomain, coeff_map, img_gens)
    if check_for_mapping_of_vars
      result.variable_indices = __maps_variables_to_variables(img_gens, codomain)
    end
    return result
  end
end

function MPolyAnyMap(d::D, c::C, cm::U, ig::Vector{V}) where {D, C, U, V}
  return MPolyAnyMap{D, C, U, V, HomUngraded}(d, c, cm, ig)
end

function MPolyAnyMap(::Type{K}, d::D, c::C, cm::U, ig::Vector{V}) where {K <: MPolyHomKind, D, C, U, V}
  return MPolyAnyMap{D, C, U, V, K}(d, c, cm, ig)
end
