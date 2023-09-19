# This type should not be exported
const BasicTypeUnion = Union{String, QQFieldElem, Symbol,
                       Number, ZZRingElem, TropicalSemiringElem}
function save_object(s::SerializerState, x::T) where T <: Union{BasicTypeUnion, VersionNumber}
  save_data_basic(s, x)
end

################################################################################
# Bool
@register_serialization_type Bool 

function load_object(s::DeserializerState, ::Type{Bool}, str::String)
  if str == "true"
    return true
  end

  if str == "false"
    return false
  end

  error("Error parsing boolean string: $str")
end

################################################################################
# ZZRingElem
@register_serialization_type ZZRingElem

function load_object(s::DeserializerState, ::Type{ZZRingElem}, str::String)
  return ZZRingElem(str)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{ZZRingElem},
                                   str::String,
                                   parent::ZZRing)
  return parent(ZZRingElem(str))
end

################################################################################
# QQFieldElem
@register_serialization_type QQFieldElem

function load_object(s::DeserializerState, ::Type{QQFieldElem}, q::String)
  # TODO: simplify the code below once https://github.com/Nemocas/Nemo.jl/pull/1375
  # is merged and in a Nemo release
  fraction_parts = collect(map(String, split(q, "//")))
  fraction_parts = [ZZRingElem(s) for s in fraction_parts]

  return QQFieldElem(fraction_parts...)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{QQFieldElem},
                                   str::String,
                                   parent::QQField)
  return parent(load_internal(s, QQFieldElem, str))
end

################################################################################
# Number
@register_serialization_type Int8
@register_serialization_type Int16
@register_serialization_type Int32
@register_serialization_type Int64 "Base.Int"
@register_serialization_type Int128

@register_serialization_type UInt8
@register_serialization_type UInt16
@register_serialization_type UInt32
@register_serialization_type UInt64
@register_serialization_type UInt128

@register_serialization_type BigInt

@register_serialization_type Float16
@register_serialization_type Float32
@register_serialization_type Float64

function load_object(s::DeserializerState, ::Type{T}, str::String) where {T<:Number}
  return parse(T, str)
end

################################################################################
# Strings
@register_serialization_type String

function load_object(s::DeserializerState, ::Type{String}, str::String)
  return str
end

################################################################################
# Symbol
@register_serialization_type Symbol

function load_object(s::DeserializerState, ::Type{Symbol}, str::String)
  return Symbol(str)
end

