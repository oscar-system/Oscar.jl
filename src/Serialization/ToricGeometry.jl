################################################################################
# Toric varieties
@register_serialization_type AffineNormalToricVariety uses_id

@register_serialization_type NormalToricVariety uses_id [:cox_ring, :class_group, :cohomology_ring]

# function type_params(ntv::T) where T <: NormalToricVarietyType
#   attrs = attrs_list(T)
#   if !isempty(attrs) && any([has_attribute(ntv, attr) for attr in attrs])
#     dict = Dict{Symbol, Any}()
# 
#     for attr in filter(x -> has_attribute(ntv, x), attrs)
#       params = type_params(get_attribute(ntv, attr))
#       if !isnothing(params)
#         dict[attr] = params
#       end
#     end
#     return Dict(:attrs => dict)
#   else
#     return nothing
#   end
# end

function save_object(s::SerializerState, ntv::T) where T <: NormalToricVarietyType
  save_object(s, ntv.polymakeNTV)
end

function load_object(s::DeserializerState, ::Type{T}) where {T <: Union{NormalToricVariety, AffineNormalToricVariety}}
  return T(load_object(s, Polymake.BigObject))
end

################################################################################
# Torus invariant divisors on toric varieties
@register_serialization_type ToricDivisor uses_params

type_params(obj::ToricDivisor) = toric_variety(obj)

function save_object(s::SerializerState, td::ToricDivisor)
  save_object(s, td.coeffs)
end

function load_object(s::DeserializerState, ::Type{ToricDivisor}, tv::NormalToricVarietyType)
  coeffs = load_object(s, Vector{ZZRingElem})
  all = Polymake._lookup_multi(pm_object(tv), "DIVISOR")
  index = 0
  for i in 1:length(all)
    if Vector{ZZRingElem}(all[i].COEFFICIENTS) == coeffs
      index = i
      break
    end
  end
  pmdiv = Polymake._lookup_multi(pm_object(tv), "DIVISOR", index-1)
  return ToricDivisor(pmdiv, tv, coeffs)
end

################################################################################
# Torus invariant divisor classes on toric varieties
@register_serialization_type ToricDivisorClass uses_params

type_params(obj::ToricDivisorClass) = toric_variety(obj)

function save_object(s::SerializerState, tdc::ToricDivisorClass)
  save_object(s, toric_divisor(tdc).coeffs)
end

function load_object(s::DeserializerState, ::Type{ToricDivisorClass}, tv::NormalToricVarietyType)
  coeffs = load_object(s, Vector{ZZRingElem})
  all = Polymake._lookup_multi(pm_object(tv), "DIVISOR")
  index = 0
  for i in 1:length(all)
    if Vector{ZZRingElem}(all[i].COEFFICIENTS) == coeffs
      index = i
      break
    end
  end
  pmdiv = Polymake._lookup_multi(pm_object(tv), "DIVISOR", index-1)
  return toric_divisor_class(ToricDivisor(pmdiv, tv, coeffs))
end

################################################################################
# Cohomology classes on toric varieties
@register_serialization_type CohomologyClass

type_params(obj::CohomologyClass) = toric_variety(obj)

function save_object(s::SerializerState, cc::CohomologyClass)
  save_data_dict(s) do
    save_object(s, lift(polynomial(cc)), :polynomial)
  end
end

function load_object(s::DeserializerState, ::Type{CohomologyClass}, tv::NormalToricVarietyType)
  poly = load_object(s, MPolyDecRingElem, base_ring(cohomology_ring(tv)), :polynomial)
  return cohomology_class(tv, MPolyQuoRingElem(poly, cohomology_ring(tv)))
end
