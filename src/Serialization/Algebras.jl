################################################################################
# FreeAssociativeAlgebra

# Free associative algebra serialization
@register_serialization_type FreeAssociativeAlgebra uses_id

type_and_params(R::T) where T <: FreeAssociativeAlgebra = TypeAndParams(T, coefficient_ring(R))

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

# see type_and_params in Rings

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
    loaded_coeff = load_object(s, coeff_type, base_ring(parent_algebra), 2)
    loaded_term = parent_algebra(loaded_coeff)
    e = load_array_node(s, 1; entry_type=Int) do _
      load_object(s, Int)
    end
    push_term!(elem, loaded_coeff, e)
  end

  return finish(elem)
end

# Ideals
@register_serialization_type FreeAssociativeAlgebraIdeal

################################################################################
# QuaternionAlgebra

# uses_id, so that the elements of one file share their parent: elements of
# distinct quaternion algebras never compare equal, and quaternion_algebra does
# not cache
@register_serialization_type Hecke.QuaternionAlgebra "QuaternionAlgebra" uses_id

type_and_params(A::Hecke.QuaternionAlgebra) = TypeAndParams(typeof(A), base_ring(A))

# the pair (a, b) of its standard form, both elements of the base field
function save_object(s::SerializerState, A::Hecke.QuaternionAlgebra)
  save_object(s, Hecke.standard_form(A))
end

function load_object(s::DeserializerState, ::Type{<:Hecke.QuaternionAlgebra}, K::Field)
  E = elem_type(K)
  a, b = load_object(s, Tuple{E, E}, (K, K))
  return quaternion_algebra(K, a, b)
end

# Quaternion algebra element serialization
@register_serialization_type Hecke.AssociativeAlgebraElem{
  <:RingElem, <:Hecke.QuaternionAlgebra} "QuaternionAlgebraElem"

type_and_params(x::Hecke.AssociativeAlgebraElem{<:RingElem, <:Hecke.QuaternionAlgebra}) =
  TypeAndParams(typeof(x), parent(x))

function save_object(
  s::SerializerState, x::Hecke.AssociativeAlgebraElem{<:RingElem, <:Hecke.QuaternionAlgebra}
)
  save_object(s, coefficients(x))
end

function load_object(
  s::DeserializerState,
  ::Type{<:Hecke.AssociativeAlgebraElem{<:RingElem, <:Hecke.QuaternionAlgebra}},
  A::Hecke.QuaternionAlgebra
)
  K = base_ring(A)
  return A(load_object(s, Vector{elem_type(K)}, K))
end
