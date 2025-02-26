###############################################################################
#
#   Lie algebras
#
###############################################################################

const lie_algebra_serialization_attributes = [
  :is_abelian, :is_nilpotent, :is_perfect, :is_semisimple, :is_simple, :is_solvable
]

@register_serialization_type AbstractLieAlgebra uses_id lie_algebra_serialization_attributes
@register_serialization_type LinearLieAlgebra uses_id [
  lie_algebra_serialization_attributes;
  [:type, :form]
]
@register_serialization_type DirectSumLieAlgebra uses_id lie_algebra_serialization_attributes

function type_params(L::T) where {T<:AbstractLieAlgebra}
  TypeParams(
    T,
    :base_ring => coefficient_ring(L),
    type_params_for_root_system(L)...,
  )
end

function save_object(s::SerializerState, L::AbstractLieAlgebra)
  save_data_dict(s) do
    save_object(s, _struct_consts(L), :struct_consts)
    save_object(s, symbols(L), :symbols)
    save_root_system_data(s, L)
  end
end

function load_object(s::DeserializerState, ::Type{<:AbstractLieAlgebra}, d::Dict)
  R = d[:base_ring]
  struct_consts = load_object(s, Matrix{sparse_row_type(R)}, R, :struct_consts)
  symbs = load_object(s, Vector{Symbol}, :symbols)
  L = lie_algebra(R, struct_consts, symbs; check=false)
  load_root_system_data(s, L, d)
  return L
end

function type_params(L::T) where {T<:LinearLieAlgebra}
  TypeParams(
    T,
    :base_ring => coefficient_ring(L),
    type_params_for_root_system(L)...,
  )
end

function save_object(s::SerializerState, L::LinearLieAlgebra)
  save_data_dict(s) do
    save_object(s, L.n, :n)
    save_object(s, matrix_repr_basis(L), :basis)
    save_object(s, symbols(L), :symbols)
    save_root_system_data(s, L)
  end
end

function load_object(s::DeserializerState, ::Type{<:LinearLieAlgebra}, d::Dict)
  R = d[:base_ring]
  n = load_object(s, Int, :n)
  basis = Vector{dense_matrix_type(R)}(
    load_object(s, Vector{dense_matrix_type(R)}, matrix_space(R, n, n), :basis)
  ) # coercion needed due to https://github.com/oscar-system/Oscar.jl/issues/3983
  symbs = load_object(s, Vector{Symbol}, :symbols)
  L = lie_algebra(R, n, basis, symbs; check=false)
  load_root_system_data(s, L, d)
  return L
end

function type_params(L::T) where {T<:DirectSumLieAlgebra}
  TypeParams(
    T,
    :base_ring => coefficient_ring(L),
    :summands => Tuple(L.summands),
    type_params_for_root_system(L)...,
  )
end

function save_object(s::SerializerState, L::DirectSumLieAlgebra)
  save_data_dict(s) do
    # summands get saved in type params
    save_root_system_data(s, L)
  end
end

function load_object(s::DeserializerState, ::Type{<:DirectSumLieAlgebra}, d::Dict)
  R = d[:base_ring]
  summands = d[:summands]

  L = direct_sum(R, collect(summands))
  load_root_system_data(s, L, d)
  return L
end

function save_root_system_data(s::SerializerState, L::LieAlgebra)
  if has_root_system(L)
    save_object(s, chevalley_basis(L), :chevalley_basis)
  end
end

function type_params_for_root_system(L::LieAlgebra)
  if has_root_system(L)
    return (:root_system => root_system(L),)
  else
    return ()
  end
end

function load_root_system_data(s::DeserializerState, L::LieAlgebra, d::Dict)
  if haskey(d, :root_system)
    rs = d[:root_system]
    chev = NTuple{3,Vector{elem_type(L)}}(
      load_object(
        s, NTuple{3,Vector{elem_type(L)}}, (L, L, L), :chevalley_basis
      ),
    ) # coercion needed due to https://github.com/oscar-system/Oscar.jl/issues/3983
    set_root_system_and_chevalley_basis!(L, rs, chev)
  end
end

@register_serialization_type AbstractLieAlgebraElem
@register_serialization_type LinearLieAlgebraElem
@register_serialization_type DirectSumLieAlgebraElem

type_params(x::T) where {T<:LieAlgebraElem} = TypeParams(T, parent(x))

function save_object(s::SerializerState, x::LieAlgebraElem)
  save_object(s, coefficients(x))
end

function load_object(s::DeserializerState, ::Type{<:LieAlgebraElem}, L::LieAlgebra)
  R = coefficient_ring(L)
  return L(load_object(s, Vector{elem_type(R)}, R))
end

###############################################################################
#
#   Lie algebra modules
#
###############################################################################

@register_serialization_type LieAlgebraModule uses_id

type_params(V::LieAlgebraModule) = TypeParams(
  LieAlgebraModule,
  :lie_algebra => base_lie_algebra(V),
  type_params_for_construction_data(V)...,
)

function save_object(s::SerializerState, V::LieAlgebraModule)
  save_data_dict(s) do
    save_object(s, dim(V), :dim)
    save_object(s, transformation_matrices(V), :transformation_matrices)
    save_object(s, symbols(V), :symbols)
  end
end

function load_object(s::DeserializerState, T::Type{<:LieAlgebraModule}, d::Dict)
  L = d[:lie_algebra]
  R = coefficient_ring(L)
  dim = load_object(s, Int, :dim)
  transformation_matrices = load_object(
    s, Vector{dense_matrix_type(R)}, matrix_space(R, dim, dim), :transformation_matrices
  )
  symbs = load_object(s, Vector{Symbol}, :symbols)
  V = load_construction_data(T, d)
  if isnothing(V)
    V = abstract_module(L, dim, transformation_matrices, symbs; check=false)
  end
  return V
end

function type_params_for_construction_data(V::LieAlgebraModule)
  if _is_standard_module(V)
    return (:construction_data => (:is_standard_module, base_lie_algebra(V)),)
  elseif ((fl, W) = _is_dual(V); fl)
    return (:construction_data => (:is_dual, W),)
  elseif ((fl, Vs) = _is_direct_sum(V); fl)
    return (:construction_data => (:is_direct_sum, Vs),)
  elseif ((fl, Vs) = _is_tensor_product(V); fl)
    return (:construction_data => (:is_tensor_product, Vs),)
  elseif ((fl, W, k) = _is_exterior_power(V); fl)
    return (:construction_data => (:is_exterior_power, W, k),)
  elseif ((fl, W, k) = _is_symmetric_power(V); fl)
    return (:construction_data => (:is_symmetric_power, W, k),)
  elseif ((fl, W, k) = _is_tensor_power(V); fl)
    return (:construction_data => (:is_tensor_power, W, k),)
  end
  return ()
end

function load_construction_data(T::Type{<:LieAlgebraModule}, d::Dict)
  V = nothing
  if haskey(d, :construction_data)
    data = d[:construction_data]
    if data[1] == :is_standard_module
      _, L = data
      V = standard_module(L)
    elseif data[1] == :is_dual
      _, W = data
      V = dual(W)
    elseif data[1] == :is_direct_sum
      _, Vs = data
      V = direct_sum(Vs...)
    elseif data[1] == :is_tensor_product
      _, Vs = data
      V = tensor_product(Vs...)
    elseif data[1] == :is_exterior_power
      _, W, k = data
      V = exterior_power(W, k)[1]
    elseif data[1] == :is_symmetric_power
      _, W, k = data
      V = symmetric_power(W, k)[1]
    elseif data[1] == :is_tensor_power
      _, W, k = data
      V = tensor_power(W, k)[1]
    end
  end
  return V::Union{T,Nothing}
end

@register_serialization_type LieAlgebraModuleElem

type_params(x::T) where {T<:LieAlgebraModuleElem} = TypeParams(T, parent(x))

function save_object(s::SerializerState, x::LieAlgebraModuleElem)
  save_object(s, coefficients(x))
end

function load_object(
  s::DeserializerState, ::Type{<:LieAlgebraModuleElem}, V::LieAlgebraModule
)
  R = coefficient_ring(V)
  return V(load_object(s, Vector{elem_type(R)}, R))
end
