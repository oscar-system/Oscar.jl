################################################################################
# UUID
@register_serialization_type UUID "UUID"
save_object(s::SerializerState, x::UUID) = save_data_basic(s, x)
load_object(s::DeserializerState, ::Type{UUID}) = load_ref(s, UUID(s.obj))

################################################################################

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

load_object(s::DeserializerState, T::Type{ZZRingElem}, ::ZZRing) = load_object(s, T)

function load_object(s::DeserializerState, ::Type{ZZRingElem})
  load_node(s) do str
    return ZZRingElem(str)
  end
end

################################################################################
# QQFieldElem
@register_serialization_type QQFieldElem

load_object(s::DeserializerState, T::Type{QQFieldElem}, ::QQField) = load_object(s, T)

function load_object(s::DeserializerState, ::Type{QQFieldElem})
  # TODO: simplify the code below once https://github.com/Nemocas/Nemo.jl/pull/1375
  # is merged and in a Nemo release
  load_node(s) do q
    fraction_parts = String.(split(q, "//"))
    fraction_parts = parse.(ZZRingElem, fraction_parts)
    
    return QQFieldElem(fraction_parts...)
  end
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

function load_object(s::DeserializerState, ::Type{T}) where {T<:Number}
  load_node(s) do str
    parse(T, str)
  end
end

@register_serialization_type PosInf

################################################################################
# Strings
@register_serialization_type String

function load_object(s::DeserializerState, ::Type{String})
  return s.obj
end

################################################################################
# Symbol
@register_serialization_type Symbol

function load_object(s::DeserializerState, ::Type{Symbol})
  load_node(s) do str
    Symbol(str)
  end
end

