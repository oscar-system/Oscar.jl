###############################################################################
#
#   Lie algebras
#
###############################################################################

const lie_algebra_serialization_attributes = [
  :is_abelian, :is_nilpotent, :is_perfect, :is_simple, :is_solvable
]

@register_serialization_type AbstractLieAlgebra uses_id lie_algebra_serialization_attributes

type_params(L::AbstractLieAlgebra) = TypeParams(
  AbstractLieAlgebra,
  :base_ring => coefficient_ring(L),
)

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
  load_root_system_data(s, L)
  return L
end

@register_serialization_type LinearLieAlgebra uses_id [
  lie_algebra_serialization_attributes;
  [:type, :form]
]

type_params(L::LinearLieAlgebra) = TypeParams(
  LinearLieAlgebra,
  :base_ring => coefficient_ring(L),
)

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
  load_root_system_data(s, L)
  return L
end

@register_serialization_type DirectSumLieAlgebra uses_id lie_algebra_serialization_attributes

type_params(L::DirectSumLieAlgebra) = TypeParams(
  DirectSumLieAlgebra,
  :base_ring => coefficient_ring(L),
)

function save_object(s::SerializerState, L::DirectSumLieAlgebra)
  save_data_dict(s) do
    save_data_array(s, :summands) do
      for summand in L.summands
        ref = save_as_ref(s, summand)
        save_object(s, ref)
      end
    end
    save_root_system_data(s, L)
  end
end

function load_object(s::DeserializerState, ::Type{<:DirectSumLieAlgebra}, d::Dict)
  R = d[:base_ring]
  summands = Vector{LieAlgebra{elem_type(R)}}(
    load_array_node(s, :summands) do _
      load_ref(s)
    end,
  )

  L = direct_sum(R, summands)
  load_root_system_data(s, L)
  return L
end

function save_root_system_data(s::SerializerState, L::LieAlgebra)
  if has_root_system(L)
    save_object(s, save_as_ref(s, root_system(L)), :root_system)
    save_object(s, chevalley_basis(L), :chevalley_basis)
  end
end

function load_root_system_data(s::DeserializerState, L::LieAlgebra)
  if haskey(s, :root_system)
    rs = load_node(_ -> load_ref(s), s, :root_system)
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
)

function save_object(s::SerializerState, V::LieAlgebraModule)
  save_data_dict(s) do
    save_object(s, dim(V), :dim)
    save_object(s, transformation_matrices(V), :transformation_matrices)
    save_object(s, symbols(V), :symbols)
    save_construction_data(s, V)
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
  V = load_construction_data(s, T)
  if isnothing(V)
    V = abstract_module(L, dim, transformation_matrices, symbs; check=false)
  end
  return V
end

function save_construction_data(s::SerializerState, V::LieAlgebraModule)
  with_attrs(s) || return nothing # attribute saving disabled
  save_data_dict(s, :construction_data) do
    if _is_standard_module(V)
      save_object(s, save_as_ref(s, base_lie_algebra(V)), :is_standard_module)
    elseif ((fl, W) = _is_dual(V); fl)
      save_object(s, save_as_ref(s, W), :is_dual)
    elseif ((fl, Vs) = _is_direct_sum(V); fl)
      save_object(s, Vs, :is_direct_sum)
    elseif ((fl, Vs) = _is_tensor_product(V); fl)
      save_object(s, Vs, :is_tensor_product)
    elseif ((fl, W, k) = _is_exterior_power(V); fl)
      save_object(s, (W, k), :is_exterior_power)
    elseif ((fl, W, k) = _is_symmetric_power(V); fl)
      save_object(s, (W, k), :is_symmetric_power)
    elseif ((fl, W, k) = _is_tensor_power(V); fl)
      save_object(s, (W, k), :is_tensor_power)
    end
  end
end

function load_construction_data(s::DeserializerState, T::Type{<:LieAlgebraModule})
  V = nothing
  with_attrs(s) && haskey(s, :construction_data) &&
    load_node(s, :construction_data) do _
      if haskey(s, :is_standard_module)
        L = load_node(s, :is_standard_module) do _
          load_ref(s)
        end
        V = standard_module(L)
      elseif haskey(s, :is_dual)
        W = load_node(s, :is_dual) do _
          load_ref(s)
        end
        V = dual(W)
      elseif haskey(s, :is_direct_sum)
        Vs = load_array_node(s, :is_direct_sum) do _
          load_ref(s)
        end
        V = direct_sum(Vs...)
      elseif haskey(s, :is_tensor_product)
        Vs = load_array_node(s, :is_tensor_product) do _
          load_ref(s)
        end
        V = tensor_product(Vs...)
      elseif haskey(s, :is_exterior_power)
        W, k = load_node(s, :is_exterior_power) do _
          (load_node(_ -> load_ref(s), s, 1), load_object(s, Int, 2))
        end
        V = exterior_power(W, k)[1]
      elseif haskey(s, :is_symmetric_power)
        W, k = load_node(s, :is_symmetric_power) do _
          (load_node(_ -> load_ref(s), s, 1), load_object(s, Int, 2))
        end
        V = symmetric_power(W, k)[1]
      elseif haskey(s, :is_tensor_power)
        W, k = load_node(s, :is_tensor_power) do _
          (load_node(_ -> load_ref(s), s, 1), load_object(s, Int, 2))
        end
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
