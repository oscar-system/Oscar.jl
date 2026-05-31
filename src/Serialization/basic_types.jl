################################################################################
# Nothing
@register_serialization_type Nothing 


function save_object(s::SerializerState, ::Nothing)
  save_data_json(s, JSON.json(nothing))
end

function load_object(s::DeserializerState, ::Type{Nothing})
  return nothing
end
################################################################################
  
function save_object(s::SerializerState, x::T) where T <: Union{BasicTypeUnion, VersionNumber}
  save_data_basic(s, x)
end

################################################################################
# Bool
@register_serialization_type Bool 

function load_object(s::DeserializerState, ::Type{Bool})
  load_node(s) do
    return load_json(s, Bool)
  end::Bool
end


################################################################################
# ZZRingElem
@register_serialization_type ZZRingElem

load_object(s::DeserializerState, ::TypeParams{ZZRingElem, ZZRing}) = load_object(s, ZZRingElem)

function load_object(s::DeserializerState, ::Type{ZZRingElem})
  load_node(s) do
    return ZZRingElem(load_json(s, String))
  end::ZZRingElem
end

################################################################################
# QQFieldElem
@register_serialization_type QQFieldElem

load_object(s::DeserializerState, ::TypeParams{QQFieldElem, QQField}) = load_object(s, QQFieldElem)

function load_object(s::DeserializerState, ::Type{QQFieldElem})
  load_node(s) do
    fraction_parts = split(load_json(s, String), "//")
    fraction_parts = parse.(ZZRingElem, fraction_parts)

    return QQFieldElem(fraction_parts...)
  end::QQFieldElem
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
  load_node(s) do
    parse(T, load_json(s, String))
  end
end

@register_serialization_type PosInf

################################################################################
# Strings
@register_serialization_type String

function load_object(s::DeserializerState, ::Type{String})
  return load_json(s, String)
end

################################################################################
# Symbol
@register_serialization_type Symbol

function load_object(s::DeserializerState, ::Type{Symbol})
  load_node(s) do
    Symbol(load_json(s, String))
  end
end

