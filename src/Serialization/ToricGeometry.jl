################################################################################
# Toric varieties
@registerSerializationType(AffineNormalToricVariety, true)
@registerSerializationType(NormalToricVariety, true)

function save_object(s::SerializerState, ntv::NormalToricVarietyType)
  save_object(s, ntv.polymakeNTV)
end

function load_object(s::DeserializerState, ::Type{T}, dict::Dict, field::Field) where {T <: Union{NormalToricVariety, AffineNormalToricVariety}}
  return T(load_object(s, Polymake.BigObject, dict))
end

################################################################################
# Torus invariant divisors on toric varieties
@registerSerializationType(ToricDivisor)
type_needs_params(::Type{<:ToricDivisor}) = true

function save_type_params(s::SerializerState, obj::ToricDivisor)
  data_dict(s) do
    save_object(s, encode_type(ToricDivisor), :name)
    save_typed_object(s, obj.toric_variety, :params)
  end
end

function load_type_params(s::DeserializerState, ::Type{<:ToricDivisor}, str::String)
  return load_ref(s, str)
end

function save_object(s::SerializerState, td::ToricDivisor)
  save_object(s, td.coeffs)
end

function load_object(s::DeserializerState, ::Type{ToricDivisor},
                     entries::Vector, tv::NormalToricVarietyType)
  coeffs = load_object(s, Vector, entries, ZZRingElem)
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
