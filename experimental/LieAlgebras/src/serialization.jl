@register_serialization_type AbstractLieAlgebra uses_id

function save_object(s::SerializerState, L::AbstractLieAlgebra)
  save_data_dict(s) do
    save_typed_object(s, coefficient_ring(L), :base_ring)
    save_object(s, _struct_consts(L), :struct_consts)
    save_object(s, symbols(L), :symbols)
    save_attrs(s, L)
  end
end

function load_object(s::DeserializerState, ::Type{AbstractLieAlgebra})
  R = load_typed_object(s, :base_ring)
  struct_consts = load_object(s, Matrix, (sparse_row_type(R), R), :struct_consts)
  symbs = load_object(s, Vector, Symbol, :symbols)
  L = lie_algebra(R, struct_consts, symbs; check=false)
  load_attrs!(s, L)
  return L
end

@register_serialization_type LinearLieAlgebra uses_id

function save_object(s::SerializerState, L::LinearLieAlgebra)
  save_data_dict(s) do
    save_typed_object(s, coefficient_ring(L), :base_ring)
    save_object(s, L.n, :n)
    save_object(s, matrix_repr_basis(L), :basis)
    save_object(s, symbols(L), :symbols)
    save_attrs(s, L)
  end
end

function load_object(s::DeserializerState, ::Type{LinearLieAlgebra})
  R = load_typed_object(s, :base_ring)
  n = load_object(s, Int, :n)
  basis = load_object(s, Vector, (dense_matrix_type(R), matrix_space(R, n, n)), :basis)
  symbs = load_object(s, Vector, Symbol, :symbols)
  L = lie_algebra(R, n, basis, symbs; check=false)
  load_attrs!(s, L)
  return L
end

@register_serialization_type DirectSumLieAlgebra uses_id

function save_object(s::SerializerState, L::DirectSumLieAlgebra)
  save_data_dict(s) do
    save_typed_object(s, coefficient_ring(L), :base_ring)
    save_data_array(s, :summands) do
      for summand in L.summands
        ref = Oscar.save_as_ref(s, summand)
        save_object(s, ref)
      end
    end
    save_attrs(s, L)
  end
end

function load_object(s::DeserializerState, ::Type{DirectSumLieAlgebra})
  R = load_typed_object(s, :base_ring)
  summands = Vector{LieAlgebra{elem_type(R)}}(
    Oscar.load_array_node(s, :summands) do _
      Oscar.load_ref(s)
    end,
  )

  L = direct_sum(R, summands)
  load_attrs!(s, L)
  return L
end

function save_attrs(s::SerializerState, L::LieAlgebra)
  save_data_dict(s, :attrs) do
    for bool_prop in (:is_abelian, :is_nilpotent, :is_perfect, :is_simple, :is_solvable)
      if has_attribute(L, bool_prop)
        save_object(s, get_attribute(L, bool_prop), bool_prop)
      end
    end
    for symbol_prop in (:type,)
      if has_attribute(L, symbol_prop)
        save_object(s, get_attribute(L, symbol_prop), symbol_prop)
      end
    end
    # TODO: handle root_system
  end
end

function load_attrs!(s::DeserializerState, L::LieAlgebra)
  Oscar.load_node(s, :attrs) do _
    for bool_prop in (:is_abelian, :is_nilpotent, :is_perfect, :is_simple, :is_solvable)
      if haskey(s, bool_prop)
        set_attribute!(L, bool_prop, load_object(s, Bool, bool_prop))
      end
    end
    for symbol_prop in (:type,)
      if haskey(s, symbol_prop)
        set_attribute!(L, symbol_prop, load_object(s, Symbol, symbol_prop))
      end
    end
  end
end

@register_serialization_type AbstractLieAlgebraElem uses_params
@register_serialization_type LinearLieAlgebraElem uses_params
@register_serialization_type DirectSumLieAlgebraElem uses_params

function save_object(s::SerializerState, x::LieAlgebraElem)
  save_object(s, coefficients(x))
end

function load_object(s::DeserializerState, ::Type{<:LieAlgebraElem}, L::LieAlgebra)
  return L(load_object(s, Vector, coefficient_ring(L)))
end

function save_type_params(s::SerializerState, x::LieAlgebraElem)
  save_data_dict(s) do
    save_object(s, Oscar.encode_type(typeof(x)), :name)
    parent_x = parent(x)
    parent_ref = Oscar.save_as_ref(s, parent_x)
    save_object(s, parent_ref, :params)
  end
end

function load_type_params(s::DeserializerState, ::Type{<:LieAlgebraElem})
  return load_typed_object(s)
end
