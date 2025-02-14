################################################################################
# Toric varieties
@register_serialization_type AffineNormalToricVariety uses_id

@register_serialization_type NormalToricVariety uses_id [:cox_ring, :class_group, :cohomology_ring]


function save_object(s::SerializerState, ntv::T) where T <: NormalToricVarietyType
  attrs = attrs_list(s, T)

  if !isempty(attrs) && any([has_attribute(ntv, attr) for attr in attrs])
    save_data_dict(s) do
      save_attrs(s, ntv)
      save_object(s, ntv.polymakeNTV, :pm_data)
    end
  else
    save_object(s, ntv.polymakeNTV)
  end
end

function load_object(s::DeserializerState, ::Type{T}) where {T <: Union{NormalToricVariety, AffineNormalToricVariety}}
  if haskey(s, :pm_data)
    ntv = T(load_object(s, Polymake.BigObject, :pm_data))
    load_attrs(s, ntv)
    
    return ntv
  end
  return T(load_object(s, Polymake.BigObject))
end

################################################################################
# Torus invariant divisors on toric varieties
@register_serialization_type ToricDivisor uses_params

function save_type_params(s::SerializerState, obj::ToricDivisor)
  save_data_dict(s) do
    save_object(s, encode_type(ToricDivisor), :name)
    save_typed_object(s, obj.toric_variety, :params)
  end
end

function load_type_params(s::DeserializerState, ::Type{<:ToricDivisor})
  return load_typed_object(s)
end

function save_object(s::SerializerState, td::ToricDivisor)
  save_object(s, td.coeffs)
end

function load_object(s::DeserializerState, ::Type{ToricDivisor}, tv::NormalToricVarietyType)
  coeffs = load_object(s, Vector, ZZRingElem)
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

function save_type_params(s::SerializerState, obj::ToricDivisorClass)
  save_data_dict(s) do
    save_object(s, encode_type(ToricDivisorClass), :name)
    save_typed_object(s, obj.toric_variety, :params)
  end
end

function load_type_params(s::DeserializerState, ::Type{<:ToricDivisorClass})
  return load_typed_object(s)
end

function save_object(s::SerializerState, tdc::ToricDivisorClass)
  save_object(s, toric_divisor(tdc).coeffs)
end

function load_object(s::DeserializerState, ::Type{ToricDivisorClass}, tv::NormalToricVarietyType)
  coeffs = load_object(s, Vector, ZZRingElem)
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
@register_serialization_type CohomologyClass uses_params

function save_type_params(s::SerializerState, obj::CohomologyClass)
  save_data_dict(s) do
    save_object(s, encode_type(CohomologyClass), :name)
    save_typed_object(s, obj.v, :params)
  end
end

function load_type_params(s::DeserializerState, ::Type{<:CohomologyClass})
  return load_typed_object(s)
end

function save_object(s::SerializerState, cc::CohomologyClass)
  save_data_dict(s) do
    save_object(s, lift(polynomial(cc)), :polynomial)
  end
end

function load_object(s::DeserializerState, ::Type{CohomologyClass}, tv::NormalToricVarietyType)
  poly = load_object(s, MPolyDecRingElem, base_ring(cohomology_ring(tv)), :polynomial)
  return cohomology_class(tv, MPolyQuoRingElem(poly, cohomology_ring(tv)))
end
