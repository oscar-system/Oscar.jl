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
    e = load_object(s, Vector{Int}, 1)
    push_term!(elem, loaded_coeff, e)
  end

  return finish(elem)
end

# Ideals
@register_serialization_type FreeAssociativeAlgebraIdeal

################################################################################
# PBW Algebras
@register_serialization_type PBWAlgRing [:is_weyl_algebra]

type_params(x::PBWAlgRing) = TypeParams(
  PBWAlgRing,
  :base_ring => base_ring(x),
  :relations => type_params(relations(x))
)

function save_object(s::SerializerState, D::PBWAlgRing)
  save_data_dict(s::SerializerState) do
    save_object(s, relations(D), :relations)
    save_object(s, ordering(D), :ordering)
  end
end

function load_object(s::DeserializerState, ::Type{PBWAlgRing}, params::Dict)
  R = params[:base_ring]
  rels = load_object(s, MatElem{elem_type(R)}, params[:relations], :relations)
  ord =  load_object(s, MonomialOrdering, R, :ordering)
  return pbw_algebra(R, rels, ord; check=false)[1]
end

@register_serialization_type PBWAlgElem

type_params(x::PBWAlgElem) = TypeParams(PBWAlgElem, parent(x))

function save_object(s::SerializerState, f::PBWAlgElem)
  save_data_array(s) do
    for term in terms(f)
      save_data_array(s) do
        save_object(s, collect(exponents(term))[1])
        save_object(s, collect(coefficients(term))[1])
      end
    end
  end
end

function load_object(s::DeserializerState, ::Type{<:PBWAlgElem}, R::T) where T <: PBWAlgRing
  coeff_type = elem_type(coefficient_ring(R))
  elem = build_ctx(R)
  
  load_array_node(s) do _
    loaded_coeff = load_object(s, coeff_type, coefficient_ring(R), 2)
    loaded_term = R(loaded_coeff)
    e = load_object(s, Vector{Int}, 1)
    push_term!(elem, loaded_coeff, e)
  end
  return finish(elem)
end
