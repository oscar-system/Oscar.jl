################################################################################
# FreeAssociativeAlgebra

# Free associative algebra serialization
@register_serialization_type FreeAssociativeAlgebra uses_id

type_params(R::T) where T <: FreeAssociativeAlgebra = TypeParams(T, base_ring(R))

function save_object(s::SerializerState, A::FreeAssociativeAlgebra)
  save_data_dict(s) do
    save_object(s, symbols(A), :symbols)
  end
end

function load_object(s::DeserializerState, ::Type{<:FreeAssociativeAlgebra}, R::Ring)
  gens = load_object(s, Vector{Symbol}, :symbols)
  return free_associative_algebra(R, gens)[1]
end

# Free associative algebra element serialization
@register_serialization_type FreeAssociativeAlgebraElem

# see type_params in Rings

function save_object(s::SerializerState, f::FreeAssociativeAlgebraElem)
  save_data_array(s) do
    for term in terms(f)
      save_data_array(s) do
        save_object(s, collect(exponent_words(term))[1])
        save_object(s, collect(coefficients(term))[1])
      end
    end
  end
end

function load_object(s::DeserializerState, ::Type{<:FreeAssociativeAlgebraElem},
                     parent_algebra::FreeAssociativeAlgebra)
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
@register_serialization_type FreeAssociativeAlgebraIdeal
