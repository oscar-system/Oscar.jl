###############################################################################
#
#   Lie algebras
#
###############################################################################

using Oscar.Serialization
import Oscar.Serialization: load_object, save_object, type_and_params, parameters

const lie_algebra_serialization_attributes = [
  :is_abelian, :is_nilpotent, :is_perfect, :is_semisimple, :is_simple, :is_solvable
]

@register_serialization_type AbstractLieAlgebra uses_id lie_algebra_serialization_attributes
@register_serialization_type LinearLieAlgebra uses_id [
  lie_algebra_serialization_attributes;
  [:type, :form]
]
@register_serialization_type DirectSumLieAlgebra uses_id lie_algebra_serialization_attributes

function type_and_params(L::T) where {T<:AbstractLieAlgebra}
  TypeAndParams(
    T,
    :base_ring => coefficient_ring(L),
    type_and_params_for_root_system(L)...,
  )
end

function save_object(s::SerializerState, L::AbstractLieAlgebra)
  save_data_dict(s) do
    save_object(s, _struct_consts(L), :struct_consts)
    save_object(s, symbols(L), :symbols)
    save_root_system_data(s, L)
  end
end

function load_object(s::DeserializerState, d::TypeAndParams{<:AbstractLieAlgebra, <:Tuple{Vararg{Pair}}})
  R = d[:base_ring]
  struct_consts = load_object(s, TypeAndParams(Matrix{sparse_row_type(R)}, R), :struct_consts)
  symbs = load_object(s, Vector{Symbol}, :symbols)
  L = lie_algebra(R, struct_consts, symbs; check=false)
  load_root_system_data(s, L, d)
  return L
end

function type_and_params(L::T) where {T<:LinearLieAlgebra}
  TypeAndParams(
    T,
    :base_ring => coefficient_ring(L),
    type_and_params_for_root_system(L)...,
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

function load_object(s::DeserializerState, d::TypeAndParams{<:LinearLieAlgebra, <:Tuple{Vararg{Pair}}})
  R = d[:base_ring]
  n = load_object(s, Int, :n)
  basis = Vector{dense_matrix_type(R)}(
    load_object(s, TypeAndParams(Vector{dense_matrix_type(R)}, matrix_space(R, n, n)), :basis)
  ) # coercion needed due to https://github.com/oscar-system/Oscar.jl/issues/3983
  symbs = load_object(s, Vector{Symbol}, :symbols)
  L = lie_algebra(R, n, basis, symbs; check=false)
  load_root_system_data(s, L, d)
  return L
end

function type_and_params(L::T) where {T<:DirectSumLieAlgebra}
  TypeAndParams(
    T,
    :base_ring => coefficient_ring(L),
    :summands => Tuple(L.summands),
    type_and_params_for_root_system(L)...,
  )
end

function save_object(s::SerializerState, L::DirectSumLieAlgebra)
  save_data_dict(s) do
    # summands get saved in type params
    save_root_system_data(s, L)
  end
end

function load_object(s::DeserializerState, d::TypeAndParams{<:DirectSumLieAlgebra, <:Tuple{Vararg{Pair}}})
  R = d[:base_ring]
  summands = d[:summands]

  L = direct_sum(R, collect(LieAlgebra{elem_type(R)}, summands))
  load_root_system_data(s, L, d)
  return L
end

function save_root_system_data(s::SerializerState, L::LieAlgebra)
  if has_root_system(L)
    save_object(s, chevalley_basis(L), :chevalley_basis)
  end
end

function type_and_params_for_root_system(L::LieAlgebra)
  if has_root_system(L)
    return (:root_system => root_system(L),)
  else
    return ()
  end
end

function load_root_system_data(s::DeserializerState, L::LieAlgebra, d::TypeAndParams)
  if haskey(d, :root_system)
    rs = d[:root_system]
    chev = NTuple{3,Vector{elem_type(L)}}(
      load_object(
        s, TypeAndParams(NTuple{3,Vector{elem_type(L)}}, (L, L, L)), :chevalley_basis
      ),
    ) # coercion needed due to https://github.com/oscar-system/Oscar.jl/issues/3983
    set_root_system_and_chevalley_basis!(L, rs, chev)
  end
end

@register_serialization_type AbstractLieAlgebraElem
@register_serialization_type LinearLieAlgebraElem
@register_serialization_type DirectSumLieAlgebraElem

function save_object(s::SerializerState, x::LieAlgebraElem)
  save_object(s, coefficients(x))
end

function load_object(s::DeserializerState, tp::TypeAndParams{<:LieAlgebraElem, <:LieAlgebra})
  L = params(tp)
  R = coefficient_ring(L)
  return L(load_object(s, TypeAndParams(Vector{elem_type(R)}, R)))
end

###############################################################################
#
#   Lie algebra modules
#
###############################################################################

@register_serialization_type LieAlgebraModule uses_id

type_and_params(V::LieAlgebraModule) = TypeAndParams(
  LieAlgebraModule,
  :lie_algebra => base_lie_algebra(V),
  type_and_params_for_construction_data(V)...,
)

function save_object(s::SerializerState, V::LieAlgebraModule)
  save_data_dict(s) do
    save_object(s, dim(V), :dim)
    save_object(s, transformation_matrices(V), :transformation_matrices)
    save_object(s, symbols(V), :symbols)
  end
end

function load_object(s::DeserializerState, d::TypeAndParams{<:LieAlgebraModule, <:Tuple{Vararg{Pair}}})
  T = Serialization.type(d)
  L = d[:lie_algebra]
  R = coefficient_ring(L)
  dim = load_object(s, Int, :dim)
  transformation_matrices = load_object(
    s, TypeAndParams(Vector{dense_matrix_type(R)}, matrix_space(R, dim, dim)), :transformation_matrices
  )
  symbs = load_object(s, Vector{Symbol}, :symbols)
  V = load_construction_data(L, T, d)
  if isnothing(V)
    V = abstract_module(L, dim, transformation_matrices, symbs; check=false)
  end
  return V
end

function type_and_params_for_construction_data(V::LieAlgebraModule)
  if _is_standard_module(V)
    return (:_is_standard_module => true,)
  elseif ((fl, W) = _is_dual(V); fl)
    return (:_is_dual => W,)
  elseif ((fl, Vs) = _is_direct_sum(V); fl)
    return (:_is_direct_sum => Tuple(Vs),)
  elseif ((fl, Vs) = _is_tensor_product(V); fl)
    return (:_is_tensor_product => Tuple(Vs),)
  elseif ((fl, W, k) = _is_exterior_power(V); fl)
    return (:_is_exterior_power => (W, k),)
  elseif ((fl, W, k) = _is_symmetric_power(V); fl)
    return (:_is_symmetric_power => (W, k),)
  elseif ((fl, W, k) = _is_tensor_power(V); fl)
    return (:_is_tensor_power => (W, k),)
  end
  return ()
end

function load_construction_data(L::LieAlgebra, T::Type{<:LieAlgebraModule}, d::TypeAndParams)
  V = nothing
  if haskey(d, :_is_standard_module) && d[:_is_standard_module]
    V = standard_module(L)
  elseif haskey(d, :_is_dual)
    W = d[:_is_dual]
    V = dual(W)
  elseif haskey(d, :_is_direct_sum)
    Vs = d[:_is_direct_sum]
    V = direct_sum(Vs...)
  elseif haskey(d, :_is_tensor_product)
    Vs = d[:_is_tensor_product]
    V = tensor_product(Vs...)
  elseif haskey(d, :_is_exterior_power)
    W, k = d[:_is_exterior_power]
    V = exterior_power(W, k)[1]
  elseif haskey(d, :_is_symmetric_power)
    W, k = d[:_is_symmetric_power]
    V = symmetric_power(W, k)[1]
  elseif haskey(d, :_is_tensor_power)
    W, k = d[:_is_tensor_power]
    V = tensor_power(W, k)[1]
  end
  return V::Union{T,Nothing}
end

@register_serialization_type LieAlgebraModuleElem

function save_object(s::SerializerState, x::LieAlgebraModuleElem)
  save_object(s, coefficients(x))
end

function load_object(s::DeserializerState, tp::TypeAndParams{<:LieAlgebraModuleElem, <:LieAlgebraModule})
  V = params(tp)
  R = coefficient_ring(V)
  return V(load_object(s, TypeAndParams(Vector{elem_type(R)}, R)))
end
