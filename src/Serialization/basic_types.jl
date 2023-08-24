# This type should not be exported
const BasicTypeUnion = Union{String, QQFieldElem, Symbol,
                       Number, ZZRingElem, TropicalSemiringElem}
function save_object(s::SerializerState, x::T) where T <: Union{BasicTypeUnion, VersionNumber}
  data_basic(s, string(x))
end

################################################################################
# Bool
@registerSerializationType(Bool)

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
@registerSerializationType(ZZRingElem)

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
@registerSerializationType(QQFieldElem)

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
@registerSerializationType(Int8)
@registerSerializationType(Int16)
@registerSerializationType(Int32)
@registerSerializationType(Int64, false, "Base.Int")
@registerSerializationType(Int128)

@registerSerializationType(UInt8)
@registerSerializationType(UInt16)
@registerSerializationType(UInt32)
@registerSerializationType(UInt64)
@registerSerializationType(UInt128)

@registerSerializationType(BigInt)

@registerSerializationType(Float16)
@registerSerializationType(Float32)
@registerSerializationType(Float64)

function load_object(s::DeserializerState, ::Type{T}, str::String) where {T<:Number}
  return parse(T, str)
end

################################################################################
# Strings
@registerSerializationType(String)

function load_object(s::DeserializerState, ::Type{String}, str::String)
  return str
end

################################################################################
# Symbol
@registerSerializationType(Symbol)

function load_object(s::DeserializerState, ::Type{Symbol}, str::String)
  return Symbol(str)
end
