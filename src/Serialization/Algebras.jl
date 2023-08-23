################################################################################
# FreeAssAlgebra

# Free associative algebra serialization
@registerSerializationType(AbstractAlgebra.Generic.FreeAssAlgebra, true)

function save_internal(s::SerializerState, A::FreeAssAlgebra)
  d = Dict(
    :base_ring => save_as_ref(s, base_ring(A)),
    :symbols => save_internal(s, symbols(A)),
  )
  return d
end

function load_internal(s::DeserializerState, ::Type{<:FreeAssAlgebra}, dict::Dict)
  R = load_unknown_type(s, dict[:base_ring])
  gens = load_internal(s, Vector{Symbol}, dict[:symbols])
  return free_associative_algebra(R, gens)[1]
end

# Free associative algebra element serialization
@registerSerializationType(FreeAssAlgElem)

function save_internal(s::SerializerState, f::FreeAssAlgElem)
  d = Dict(
    :parent => save_as_ref(s, parent(f)),
    :coeffs => save_internal(s, collect(coefficients(f))),
    :exps => save_internal(s, collect(exponent_words(f))),
    :length => save_internal(s, length(f)),
  )
  return d
end

function load_internal(s::DeserializerState, ::Type{<:FreeAssAlgElem}, dict::Dict)
  parent = load_ref(s, dict[:parent])
  coeff_type=elem_type(coefficient_ring(parent))
  coeffs = load_internal(s, Vector{coeff_type}, dict[:coeffs]) 
  exps = load_internal(s, Vector{Vector{Int}},  dict[:exps]) 
  length = load_internal(s, Int, dict[:length])
  element = parent(coeffs, exps)
  return element
end

# Free associative algebra element serialization
@registerSerializationType(FreeAssAlgElem)

function save_internal(s::SerializerState, f::FreeAssAlgElem)
  d = Dict(
    :parent => save_as_ref(s, parent(f)),
    :coeffs => save_internal(s, collect(coefficients(f))),
    :exps => save_internal(s, collect(exponent_words(f))),
    :length => save_internal(s, length(f)),
  )
  return d
end

function load_internal(s::DeserializerState, ::Type{<:FreeAssAlgElem}, dict::Dict)
  parent = load_ref(s, dict[:parent])
  coeff_type=elem_type(coefficient_ring(parent))
  coeffs = load_internal(s, Vector{coeff_type}, dict[:coeffs]) 
  exps = load_internal(s, Vector{Vector{Int}},  dict[:exps]) 
  length = load_internal(s, Int, dict[:length])
  element = parent(coeffs, exps)
  return element
end
