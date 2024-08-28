###############################################################################
#
#   Lie algebras
#
###############################################################################

const lie_algebra_serialization_attributes = [
  :is_abelian, :is_nilpotent, :is_perfect, :is_simple, :is_solvable
]

@register_serialization_type AbstractLieAlgebra uses_id lie_algebra_serialization_attributes

function save_object(s::SerializerState, L::AbstractLieAlgebra)
  save_data_dict(s) do
    save_typed_object(s, coefficient_ring(L), :base_ring)
    save_object(s, _struct_consts(L), :struct_consts)
    save_object(s, symbols(L), :symbols)
    save_root_system_data(s, L)
    save_attrs(s, L)
  end
end

function load_object(s::DeserializerState, ::Type{<:AbstractLieAlgebra})
  R = load_typed_object(s, :base_ring)
  struct_consts = load_object(s, Matrix, (sparse_row_type(R), R), :struct_consts)
  symbs = load_object(s, Vector, Symbol, :symbols)
  L = lie_algebra(R, struct_consts, symbs; check=false)
  load_root_system_data(s, L)
  load_attrs(s, L)
  return L
end

@register_serialization_type LinearLieAlgebra uses_id [
  lie_algebra_serialization_attributes;
  [:type, :form]
]

function save_object(s::SerializerState, L::LinearLieAlgebra)
  save_data_dict(s) do
    save_typed_object(s, coefficient_ring(L), :base_ring)
    save_object(s, L.n, :n)
    save_object(s, matrix_repr_basis(L), :basis)
    save_object(s, symbols(L), :symbols)
    save_root_system_data(s, L)
    save_attrs(s, L)
  end
end

function load_object(s::DeserializerState, ::Type{<:LinearLieAlgebra})
  R = load_typed_object(s, :base_ring)
  n = load_object(s, Int, :n)
  basis = load_object(s, Vector, (dense_matrix_type(R), matrix_space(R, n, n)), :basis)
  symbs = load_object(s, Vector, Symbol, :symbols)
  L = lie_algebra(R, n, basis, symbs; check=false)
  load_root_system_data(s, L)
  load_attrs(s, L)
  return L
end

@register_serialization_type DirectSumLieAlgebra uses_id lie_algebra_serialization_attributes

function save_object(s::SerializerState, L::DirectSumLieAlgebra)
  save_data_dict(s) do
    save_typed_object(s, coefficient_ring(L), :base_ring)
    save_data_array(s, :summands) do
      for summand in L.summands
        ref = save_as_ref(s, summand)
        save_object(s, ref)
      end
    end
    save_root_system_data(s, L)
    save_attrs(s, L)
  end
end

function load_object(s::DeserializerState, ::Type{<:DirectSumLieAlgebra})
  R = load_typed_object(s, :base_ring)
  summands = Vector{LieAlgebra{elem_type(R)}}(
    load_array_node(s, :summands) do _
      load_ref(s)
    end,
  )

  L = direct_sum(R, summands)
  load_root_system_data(s, L)
  load_attrs(s, L)
  return L
end

function save_root_system_data(s::SerializerState, L::LieAlgebra)
  if has_root_system(L)
    save_typed_object(s, root_system(L), :root_system)
    save_object(s, chevalley_basis(L), :chevalley_basis)
  end
end

function load_root_system_data(s::DeserializerState, L::LieAlgebra)
  if haskey(s, :root_system)
    @assert L isa AbstractLieAlgebra # TODO: adapt once we have a proper interface for this
    L.root_system = load_typed_object(s, :root_system)
    chevalley_basis = load_object(
      s, Tuple, [(Vector, (AbstractLieAlgebraElem, L)) for _ in 1:3], :chevalley_basis
    )
    # chevalley basis will become an attribute in the near future
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
    save_object(s, encode_type(typeof(x)), :name)
    parent_x = parent(x)
    parent_ref = save_as_ref(s, parent_x)
    save_object(s, parent_ref, :params)
  end
end

function load_type_params(s::DeserializerState, ::Type{<:LieAlgebraElem})
  return load_typed_object(s)
end

###############################################################################
#
#   Root systems
#
###############################################################################

@register_serialization_type RootSystem uses_id

function save_object(s::SerializerState, R::RootSystem)
  save_data_dict(s) do
    save_object(s, cartan_matrix(R), :cartan_matrix)
    if has_root_system_type(R)
      type, type_ordering = root_system_type_with_ordering(R)
      save_object(s, type, :type)
      if type_ordering != 1:length(type_ordering) # don't save if it's the default
        save_object(s, type_ordering, :type_ordering)
      end
    end
  end
end

function load_object(s::DeserializerState, ::Type{RootSystem})
  cm = load_object(s, Matrix, Int, :cartan_matrix)
  R = root_system(cm; check=false, detect_type=false)
  if haskey(s, :type)
    type = Vector{Tuple{Symbol,Int}}(load_object(s, Vector, (Tuple, [Symbol, Int]), :type)) # coercion needed due to https://github.com/oscar-system/Oscar.jl/issues/3983
    if haskey(s, :type_ordering)
      type_ordering = load_object(s, Vector, Int, :type_ordering)
      set_root_system_type!(R, type, type_ordering)
    else
      set_root_system_type!(R, type)
    end
  end
  return R
end

@register_serialization_type RootSpaceElem uses_params
@register_serialization_type DualRootSpaceElem uses_params

function save_object(s::SerializerState, r::Union{RootSpaceElem,DualRootSpaceElem})
  save_object(s, _vec(coefficients(r)))
end

function load_object(
  s::DeserializerState, T::Type{<:Union{RootSpaceElem,DualRootSpaceElem}}, R::RootSystem
)
  return T(R, load_object(s, Vector, QQ))
end

function save_type_params(s::SerializerState, r::Union{RootSpaceElem,DualRootSpaceElem})
  save_data_dict(s) do
    save_object(s, encode_type(typeof(r)), :name)
    rs_x = root_system(r)
    rs_ref = save_as_ref(s, rs_x)
    save_object(s, rs_ref, :params)
  end
end

function load_type_params(
  s::DeserializerState, ::Type{<:Union{RootSpaceElem,DualRootSpaceElem}}
)
  return load_typed_object(s)
end

###############################################################################
#
#   Weyl groups
#
###############################################################################

@register_serialization_type WeylGroup uses_id

function save_object(s::SerializerState, W::WeylGroup)
  save_data_dict(s) do
    save_typed_object(s, root_system(W), :root_system)
  end
end

function load_object(s::DeserializerState, ::Type{WeylGroup})
  R = load_typed_object(s, :root_system)
  return weyl_group(R)
end

@register_serialization_type WeylGroupElem uses_params

function save_object(s::SerializerState, x::WeylGroupElem)
  save_object(s, word(x))
end

function load_object(s::DeserializerState, ::Type{WeylGroupElem}, W::WeylGroup)
  return W(load_object(s, Vector, UInt8); normalize=false)
end

function save_type_params(s::SerializerState, x::WeylGroupElem)
  save_data_dict(s) do
    save_object(s, encode_type(typeof(x)), :name)
    parent_x = parent(x)
    parent_ref = save_as_ref(s, parent_x)
    save_object(s, parent_ref, :params)
  end
end

function load_type_params(s::DeserializerState, ::Type{WeylGroupElem})
  return load_typed_object(s)
end
