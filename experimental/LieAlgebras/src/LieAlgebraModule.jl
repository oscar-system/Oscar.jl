@attributes mutable struct LieAlgebraModule{C<:RingElement}
  L::LieAlgebra{C}
  dim::Int
  transformation_matrices::Vector{MatElem{C}}
  s::Vector{Symbol}

  function LieAlgebraModule{C}(
    L::LieAlgebra{C},
    dimV::Int,
    transformation_matrices::Vector{<:MatElem{C}},
    s::Vector{Symbol};
    cached::Bool=true,
    check::Bool=true,
  ) where {C<:RingElement}
    return get_cached!(
      LieAlgebraModuleDict, (L, dimV, transformation_matrices, s), cached
    ) do
      @req dimV == length(s) "Invalid number of basis element names."
      @req dim(L) == length(transformation_matrices) "Invalid number of transformation matrices."
      @req all(m -> size(m) == (dimV, dimV), transformation_matrices) "Invalid transformation matrix dimensions."

      V = new{C}(L, dimV, transformation_matrices, s)
      if check
        for xi in basis(L), xj in basis(L), v in basis(V)
          @req (xi * xj) * v == xi * (xj * v) - xj * (xi * v) "Transformation matrices do not define a module."
        end
      end
      V
    end::LieAlgebraModule{C}
  end
end

const LieAlgebraModuleDict = CacheDictType{
  Tuple{LieAlgebra,Int,Vector{MatElem},Vector{Symbol}},LieAlgebraModule
}()

struct LieAlgebraModuleElem{C<:RingElement}
  parent::LieAlgebraModule{C}
  mat::MatElem{C}
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{LieAlgebraModuleElem{C}}) where {C<:RingElement} = LieAlgebraModule{C}

elem_type(::Type{LieAlgebraModule{C}}) where {C<:RingElement} = LieAlgebraModuleElem{C}

parent(v::LieAlgebraModuleElem{C}) where {C<:RingElement} = v.parent

base_ring(V::LieAlgebraModule{C}) where {C<:RingElement} = base_ring(base_lie_algebra(V))

base_ring(v::LieAlgebraModuleElem{C}) where {C<:RingElement} = base_ring(parent(v))

base_lie_algebra(V::LieAlgebraModule{C}) where {C<:RingElement} = V.L

ngens(L::LieAlgebraModule{C}) where {C<:RingElement} = dim(L)

gens(L::LieAlgebraModule{C}) where {C<:RingElement} = basis(L)

gen(L::LieAlgebraModule{C}, i::Int) where {C<:RingElement} = basis(L, i)

dim(V::LieAlgebraModule{C}) where {C<:RingElement} = V.dim

basis(L::LieAlgebraModule{C}) where {C<:RingElement} =
  [basis(L, i)::elem_type(L) for i in 1:dim(L)]

function basis(L::LieAlgebraModule{C}, i::Int) where {C<:RingElement}
  R = base_ring(L)
  return L([(j == i ? one(R) : zero(R)) for j in 1:dim(L)])
end

function zero(V::LieAlgebraModule{C}) where {C<:RingElement}
  mat = zero_matrix(base_ring(V), 1, dim(V))
  return elem_type(V)(V, mat)
end

function iszero(v::LieAlgebraModuleElem{C}) where {C<:RingElement}
  return iszero(Generic._matrix(v))
end

@inline function Generic._matrix(v::LieAlgebraModuleElem{C}) where {C<:RingElement}
  return (v.mat)::dense_matrix_type(C)
end

@doc raw"""
    getindex(v::LieAlgebraElem{C}, i::Int) where C <: RingElement

Return the $i$-th coefficient of the module element $x$.
"""
function getindex(v::LieAlgebraModuleElem{C}, i::Int) where {C<:RingElement}
  return Generic._matrix(v)[1, i]
end

function Base.deepcopy_internal(
  v::LieAlgebraModuleElem{C}, dict::IdDict
) where {C<:RingElement}
  return parent(v)(deepcopy_internal(Generic._matrix(v), dict))
end

function check_parent(
  v1::LieAlgebraModuleElem{C}, v2::LieAlgebraModuleElem{C}
) where {C<:RingElement}
  parent(v1) !== parent(v2) && error("Incompatible modules.")
end

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::LieAlgebraModule{C}) where {C<:RingElement}
  if has_attribute(V, :show) && get_attribute(V, :show) isa Function
    get_attribute(V, :show)(io, V)
  else
    print(io, "AbstractModule of ")
    print(IOContext(io, :compact => true), base_lie_algebra(V))
  end
end

function show_standard_module(io::IO, V::LieAlgebraModule{C}) where {C<:RingElement}
  print(io, "StdModule of ")
  print(IOContext(io, :compact => true), base_lie_algebra(V))
end

function show_exterior_power(io::IO, V::LieAlgebraModule{C}) where {C<:RingElement}
  print(io, "$(get_attribute(V, :power))-th exterior power of ")
  print(IOContext(io, :compact => true), get_attribute(V, :inner_module))
end

function show_symmetric_power(io::IO, V::LieAlgebraModule{C}) where {C<:RingElement}
  print(io, "$(get_attribute(V, :power))-th symmetric power of ")
  print(IOContext(io, :compact => true), get_attribute(V, :inner_module))
end

function show_tensor_power(io::IO, V::LieAlgebraModule{C}) where {C<:RingElement}
  print(io, "$(get_attribute(V, :power))-th tensor power of ")
  print(IOContext(io, :compact => true), get_attribute(V, :inner_module))
end

function symbols(V::LieAlgebraModule{C}) where {C<:RingElement}
  return V.s
end

function expressify(
  v::LieAlgebraModuleElem{C}, s=symbols(parent(v)); context=nothing
) where {C<:RingElement}
  sum = Expr(:call, :+)
  for (i, c) in enumerate(Generic._matrix(v))
    push!(sum.args, Expr(:call, :*, expressify(c; context=context), s[i]))
  end
  return sum
end

@enable_all_show_via_expressify LieAlgebraModuleElem

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (V::LieAlgebraModule{C})() where {C<:RingElement}
  return zero(V)
end

function (V::LieAlgebraModule{C})(v::Vector{Int}) where {C<:RingElement}
  return V(base_ring(V).(v))
end

function (V::LieAlgebraModule{C})(v::Vector{C}) where {C<:RingElement}
  @req length(v) == dim(V) "Length of vector does not match dimension."
  mat = matrix(base_ring(V), 1, length(v), v)
  return elem_type(V)(V, mat)
end

function (V::LieAlgebraModule{C})(v::MatElem{C}) where {C<:RingElement}
  @req ncols(v) == dim(V) "Length of vector does not match dimension"
  @req nrows(v) == 1 "Not a vector in module constructor"
  return elem_type(V)(V, v)
end

function (V::LieAlgebraModule{C})(v::SRow{C}) where {C<:RingElement}
  mat = dense_row(v, dim(V))
  return elem_type(V)(V, mat)
end

function (V::LieAlgebraModule{C})(v::LieAlgebraModuleElem{C}) where {C<:RingElement}
  @req V == parent(v) "Incompatible modules."
  return v
end

function (V::LieAlgebraModule{C})(
  a::Vector{T}
) where {T<:LieAlgebraModuleElem{C}} where {C<:RingElement}
  @req is_exterior_power(V) || is_symmetric_power(V) || is_tensor_power(V) "Only implemented for power modules."
  @req length(a) == get_attribute(V, :power) "Length of vector does not match power."
  @req all(x -> parent(x) == get_attribute(V, :inner_module), a) "Incompatible modules."
  mat = zero_matrix(base_ring(V), 1, dim(V))
  if is_exterior_power(V)
    for (i, _inds) in enumerate(get_attribute(V, :ind_map)),
      (inds, sgn) in permutations_with_sign(_inds)

      mat[1, i] += sgn * prod(a[j].mat[k] for (j, k) in enumerate(inds))
    end
  elseif is_symmetric_power(V)
    for (i, _inds) in enumerate(get_attribute(V, :ind_map)),
      inds in unique(permutations(_inds))

      mat[1, i] += prod(a[j].mat[k] for (j, k) in enumerate(inds))
    end
  elseif is_tensor_power(V)
    for (i, inds) in enumerate(get_attribute(V, :ind_map))
      mat[1, i] += prod(a[j].mat[k] for (j, k) in enumerate(inds))
    end
  end
  return LieAlgebraModuleElem{C}(V, mat)
end

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

function Base.:-(v::LieAlgebraModuleElem{C}) where {C<:RingElement}
  return parent(v)(-Generic._matrix(v))
end

function Base.:+(
  v1::LieAlgebraModuleElem{C}, v2::LieAlgebraModuleElem{C}
) where {C<:RingElement}
  check_parent(v1, v2)
  return parent(v1)(Generic._matrix(v1) + Generic._matrix(v2))
end

function Base.:-(
  v1::LieAlgebraModuleElem{C}, v2::LieAlgebraModuleElem{C}
) where {C<:RingElement}
  check_parent(v1, v2)
  return parent(v1)(Generic._matrix(v1) - Generic._matrix(v2))
end

function Base.:*(v::LieAlgebraModuleElem{C}, c::C) where {C<:RingElem}
  base_ring(v) != parent(c) && error("Incompatible rings.")
  return parent(v)(Generic._matrix(v) * c)
end

function Base.:*(
  v::LieAlgebraModuleElem{C}, c::U
) where {C<:RingElement,U<:Union{Rational,Integer}}
  return parent(v)(Generic._matrix(v) * c)
end

function Base.:*(c::C, v::LieAlgebraModuleElem{C}) where {C<:RingElem}
  base_ring(v) != parent(c) && error("Incompatible rings.")
  return parent(v)(c * Generic._matrix(v))
end

function Base.:*(
  c::U, v::LieAlgebraModuleElem{C}
) where {C<:RingElement,U<:Union{Rational,Integer}}
  return parent(v)(c * Generic._matrix(v))
end

###############################################################################
#
#   Comparison functions
#
###############################################################################

function Base.:(==)(
  v1::LieAlgebraModuleElem{C}, v2::LieAlgebraModuleElem{C}
) where {C<:RingElement}
  check_parent(v1, v2)
  return Generic._matrix(v1) == Generic._matrix(v2)
end

function Base.hash(v::LieAlgebraModuleElem{C}, h::UInt) where {C<:RingElement}
  b = 0x723913014484513a % UInt
  return xor(hash(Generic._matrix(v), hash(parent(v), h)), b)
end

###############################################################################
#
#   Module action
#
###############################################################################

function Base.:*(x::LieAlgebraElem{C}, v::LieAlgebraModuleElem{C}) where {C<:RingElement}
  return action(x, v)
end

function action(
  x::LieAlgebraElem{C}, v::ElemT
) where {ElemT<:LieAlgebraModuleElem{C}} where {C<:RingElement}
  @req parent(x) == base_lie_algebra(parent(v)) "Incompatible Lie algebras."

  cx = Generic._matrix(x)

  return parent(v)(
    sum(
      cx[i] * Generic._matrix(v) * transpose(transformation_matrix(parent(v), i)) for
      i in 1:dim(parent(x)) if !iszero(cx[i]);
      init=zero_matrix(base_ring(parent(v)), 1, dim(parent(v)))::dense_matrix_type(C),
    ), # equivalent to (x * v^T)^T, since we work with row vectors
  )::ElemT
end

function transformation_matrix(V::LieAlgebraModule{C}, i::Int) where {C<:RingElement}
  return (V.transformation_matrices[i])::dense_matrix_type(C)
end

###############################################################################
#
#   Attribute getters
#
###############################################################################

function is_standard_module(V::LieAlgebraModule{C}) where {C<:RingElement}
  return get_attribute(V, :type, :fallback) == :standard_module
end

function is_exterior_power(V::LieAlgebraModule{C}) where {C<:RingElement}
  return get_attribute(V, :type, :fallback) == :exterior_power
end

function is_symmetric_power(V::LieAlgebraModule{C}) where {C<:RingElement}
  return get_attribute(V, :type, :fallback) == :symmetric_power
end

function is_tensor_power(V::LieAlgebraModule{C}) where {C<:RingElement}
  return get_attribute(V, :type, :fallback) == :tensor_power
end

###############################################################################
#
#   Constructor
#
###############################################################################

function abstract_module(
  L::LieAlgebra{C},
  dimV::Int,
  transformation_matrices::Vector{<:MatElem{C}},
  s::Vector{<:VarName}=[Symbol("v_$i") for i in 1:dimV];
  cached::Bool=true,
  check::Bool=true,
) where {C<:RingElement}
  return LieAlgebraModule{C}(
    L, dimV, transformation_matrices, Symbol.(s); cached, check=check
  )
end

function abstract_module(
  L::LieAlgebra{C},
  dimV::Int,
  struct_consts::Matrix{SRow{C}},
  s::Vector{<:VarName}=[Symbol("v_$i") for i in 1:dimV];
  cached::Bool=true,
  check::Bool=true,
) where {C<:RingElement}
  @req dim(L) == size(struct_consts, 1) "Invalid structure constants dimensions."
  @req dimV == size(struct_consts, 2) "Invalid structure constants dimensions."
  @req dimV == length(s) "Invalid number of basis element names."

  transformation_matrices = [zero_matrix(base_ring(L), dimV, dimV) for _ in 1:dim(L)]
  for i in 1:dim(L), j in 1:dimV
    transformation_matrices[i][:, j] = transpose(dense_row(struct_consts[i, j], dimV))
  end

  return LieAlgebraModule{C}(
    L, dimV, transformation_matrices, Symbol.(s); cached, check=check
  )
end

function highest_weight_module(
  L::LieAlgebra{C}, weight::Vector{Int}; cached::Bool=true
) where {C<:RingElement}
  struct_consts = lie_algebra_highest_weight_module_struct_consts_gap(L, weight)
  dimV = size(struct_consts, 2)
  V = abstract_module(L, dimV, struct_consts; cached, check=false)
  set_attribute!(V, :highest_weight => weight)
  return V
end

function standard_module(L::LinearLieAlgebra{C}; cached::Bool=true) where {C<:RingElement}
  dim_std_V = L.n
  transformation_matrices = matrix_repr_basis(L)
  s = [Symbol("v_$(i)") for i in 1:dim_std_V]
  std_V = LieAlgebraModule{elem_type(base_ring(L))}(
    L, dim_std_V, transformation_matrices, s; cached, check=false
  )
  set_attribute!(std_V, :type => :standard_module, :show => show_standard_module)
  return std_V
end

function exterior_power(
  V::LieAlgebraModule{C}, k::Int; cached::Bool=true
) where {C<:RingElement}
  L = base_lie_algebra(V)
  dim_pow_V = binomial(dim(V), k)
  ind_map = collect(combinations(1:dim(V), k))

  T = tensor_power(V, k)
  basis_change_E2T = zero_matrix(base_ring(V), dim(T), dim_pow_V)
  basis_change_T2E = zero_matrix(base_ring(V), dim_pow_V, dim(T))
  T_ind_map = get_attribute(T, :ind_map)
  for (i, _inds) in enumerate(ind_map), (inds, sgn) in permutations_with_sign(_inds)
    j = findfirst(==(inds), T_ind_map)
    basis_change_E2T[j, i] = sgn//factorial(k)
    basis_change_T2E[i, j] = sgn
  end
  transformation_matrices = map(1:dim(L)) do i
    basis_change_T2E * transformation_matrix(T, i) * basis_change_E2T
  end

  if k == 1
    s = symbols(V)
  else
    if is_standard_module(V)
      parentheses = identity
    else
      parentheses = x -> "($x)"
    end
    s = [Symbol(join(parentheses.(s), " ∧ ")) for s in combinations(symbols(V), k)]
  end

  pow_V = LieAlgebraModule{C}(L, dim_pow_V, transformation_matrices, s; cached, check=false)
  set_attribute!(
    pow_V,
    :type => :exterior_power,
    :power => k,
    :inner_module => V,
    :ind_map => ind_map,
    :show => show_exterior_power,
  )
  return pow_V
end

function symmetric_power(
  V::LieAlgebraModule{C}, k::Int; cached::Bool=true
) where {C<:RingElement}
  L = base_lie_algebra(V)
  dim_pow_V = binomial(dim(V) + k - 1, k)
  ind_map = collect(multicombinations(1:dim(V), k))

  T = tensor_power(V, k)
  basis_change_S2T = zero_matrix(base_ring(V), dim(T), dim_pow_V)
  basis_change_T2S = zero_matrix(base_ring(V), dim_pow_V, dim(T))
  T_ind_map = get_attribute(T, :ind_map)
  for (i, _inds) in enumerate(ind_map), inds in permutations(_inds)
    j = findfirst(==(inds), T_ind_map)
    basis_change_S2T[j, i] += 1//factorial(k)
    basis_change_T2S[i, j] = 1
  end
  transformation_matrices = map(1:dim(L)) do i
    basis_change_T2S * transformation_matrix(T, i) * basis_change_S2T
  end

  if k == 1
    s = symbols(V)
  else
    if is_standard_module(V)
      parentheses = identity
    else
      parentheses = x -> "($x)"
    end
    s = [
      Symbol(
        join(
          (
            begin
              e = count(==(i), inds)
              if e == 1
                parentheses(s)
              else
                "$(parentheses(s))^$e"
              end
            end for (i, s) in enumerate(symbols(V)) if in(i, inds)
          ),
          "*",
        ),
      ) for inds in ind_map
    ]
  end

  pow_V = LieAlgebraModule{C}(L, dim_pow_V, transformation_matrices, s; cached, check=false)
  set_attribute!(
    pow_V,
    :type => :symmetric_power,
    :power => k,
    :inner_module => V,
    :ind_map => ind_map,
    :show => show_symmetric_power,
  )
  return pow_V
end

function tensor_power(
  V::LieAlgebraModule{C}, k::Int; cached::Bool=true
) where {C<:RingElement}
  L = base_lie_algebra(V)
  dim_pow_V = dim(V)^k
  ind_map = reverse.(collect(ProductIterator(1:dim(V), k)))

  transformation_matrices = map(1:dim(L)) do i
    y = transformation_matrix(V, i)
    sum(reduce(kronecker_product, (j == i ? y : one(y) for j in 1:k)) for i in 1:k)
  end

  if k == 1
    s = symbols(V)
  else
    if is_standard_module(V)
      parentheses = identity
    else
      parentheses = x -> "($x)"
    end
    s = [
      Symbol(join(parentheses.(s), " ⊗ ")) for s in reverse.(ProductIterator(symbols(V), k))
    ]
  end

  pow_V = LieAlgebraModule{C}(L, dim_pow_V, transformation_matrices, s; cached, check=false)
  set_attribute!(
    pow_V,
    :type => :tensor_power,
    :power => k,
    :inner_module => V,
    :ind_map => ind_map,
    :show => show_tensor_power,
  )
  return pow_V
end
