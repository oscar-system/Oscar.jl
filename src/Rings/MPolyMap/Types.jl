### Types for maps from polynomial rings
const _DomainTypes = Union{MPolyRing, MPolyQuoRing}

@attributes mutable struct MPolyAnyMap{
    D <: _DomainTypes,
    C <: NCRing,
    U,
    V} <: Map{D, C, Map, MPolyAnyMap}

  domain::D
  codomain::C
  coeff_map::U
  img_gens::Vector{V}
  temp_ring           # temporary ring used when evaluating maps

  function MPolyAnyMap{D, C, U, V}(domain::D,
                                codomain::C,
                                coeff_map::U,
                                img_gens::Vector{V}) where {D, C, U, V}
      @assert V === elem_type(C)
      for g in img_gens
        @assert parent(g) === codomain "elements does not have the correct parent"
      end
    return new{D, C, U, V}(domain, codomain, coeff_map, img_gens)
  end
end

function MPolyAnyMap(d::D, c::C, cm::U, ig::Vector{V}) where {D, C, U, V}
  return MPolyAnyMap{D, C, U, V}(d, c, cm, ig)
end
