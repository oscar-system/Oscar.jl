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
  variable_indices::Vector{Int} # a table where the i-th entry contains the 
                                # index of the variable where it is mapped to
                                # in case the mapping takes such a particularly
                                # simple form.

  function MPolyAnyMap{D, C, U, V}(domain::D,
                                   codomain::C,
                                   coeff_map::U,
                                   img_gens::Vector{V};
                                   check_for_mapping_of_vars::Bool=true
                                  ) where {D, C, U, V}
    @assert V === elem_type(C)
    for g in img_gens
      @assert parent(g) === codomain "elements does not have the correct parent"
    end
    result = new{D, C, U, V}(domain, codomain, coeff_map, img_gens)
    # If it ever turned out that doing the checks within the following if-block
    # is a bottleneck, consider passing on the `check_for_mapping_of_vars` kw 
    # argument to the outer constructors or make the outer constructors 
    # call the inner one with this argument set to `false`. This way the check 
    # can safely be disabled.
    if check_for_mapping_of_vars
      result.variable_indices = __maps_variables_to_variables(img_gens,
                                                              codomain)
    end
    return result
  end
end

function MPolyAnyMap(d::D, c::C, cm::U, ig::Vector{V}) where {D, C, U, V}
  return MPolyAnyMap{D, C, U, V}(d, c, cm, ig)
end

