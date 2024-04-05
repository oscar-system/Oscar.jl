################################################################################
# FreeAssAlgebra

# Free associative algebra serialization
@register_serialization_type FreeAssAlgebra uses_id

function save_object(s::SerializerState, A::FreeAssAlgebra)
  save_data_dict(s) do
    save_typed_object(s, base_ring(A), :base_ring),
    save_object(s, symbols(A), :symbols)
  end
end

function load_object(s::DeserializerState, ::Type{<:FreeAssAlgebra})
  R = load_typed_object(s, :base_ring)
  gens = load_object(s, Vector{Symbol}, :symbols)
  return free_associative_algebra(R, gens)[1]
end

# Free associative algebra element serialization
@register_serialization_type FreeAssAlgElem uses_params

# see save_type_params in Rings

function save_object(s::SerializerState, f::FreeAssAlgElem)
  save_data_array(s) do
    for term in terms(f)
      save_data_array(s) do
        save_object(s, collect(exponent_words(term))[1])
        save_object(s, collect(coefficients(term))[1])
      end
    end
  end
end

function load_object(s::DeserializerState, ::Type{<:FreeAssAlgElem}, parents::Vector)
  parent_algebra = parents[end]
  coeff_type = elem_type(base_ring(parent_algebra))
  elem = MPolyBuildCtx(parent_algebra)

  load_array_node(s) do _
    loaded_coeff = load_object(s, coeff_type, 2)
    loaded_term = parent_algebra(loaded_coeff)
    e = load_array_node(s, 1) do _
      load_object(s, Int)
    end
    # guarantees e is a Int[]
    e = convert(Vector{Int}, e)
    push_term!(elem, loaded_coeff, e)
  end

  return finish(elem)
end

# Ideals
@register_serialization_type FreeAssAlgIdeal uses_params
