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
    check::Bool=true,
  ) where {C<:RingElement}
    @req dimV == length(s) "Invalid number of basis element names."
    @req dim(L) == length(transformation_matrices) "Invalid number of transformation matrices."
    @req all(m -> size(m) == (dimV, dimV), transformation_matrices) "Invalid transformation matrix dimensions."

    V = new{C}(L, dimV, transformation_matrices, s)
    if check
      for xi in basis(L), xj in basis(L), v in basis(V)
        @req (xi * xj) * v == xi * (xj * v) - xj * (xi * v) "Transformation matrices do not define a module."
      end
    end
    return V
  end
end

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
  return iszero(coefficients(v))
end

@inline function Generic._matrix(v::LieAlgebraModuleElem{C}) where {C<:RingElement}
  return (v.mat)::dense_matrix_type(C)
end

function coefficients(v::LieAlgebraModuleElem{C}) where {C<:RingElement}
  return collect(Generic._matrix(v))[1, :]
end

function coeff(v::LieAlgebraModuleElem{C}, i::Int) where {C<:RingElement}
  return Generic._matrix(v)[1, i]
end

@doc raw"""
    getindex(v::LieAlgebraElem{C}, i::Int) where C <: RingElement

Return the $i$-th coefficient of the module element $x$.
"""
function getindex(v::LieAlgebraModuleElem{C}, i::Int) where {C<:RingElement}
  return coeff(v, i)
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

function show_dual(io::IO, V::LieAlgebraModule{C}) where {C<:RingElement}
  print(io, "Dual of ")
  print(IOContext(io, :compact => true), base_module(V))
end

function show_direct_sum(io::IO, V::LieAlgebraModule{C}) where {C<:RingElement}
  print(io, "Direct sum of ")
  print(IOContext(io, :compact => true), base_modules(V))
end

function show_tensor_product(io::IO, V::LieAlgebraModule{C}) where {C<:RingElement}
  print(io, "Tensor product of ")
  print(IOContext(io, :compact => true), base_modules(V))
end

function show_exterior_power(io::IO, V::LieAlgebraModule{C}) where {C<:RingElement}
  print(io, "$(get_attribute(V, :power))-th exterior power of ")
  print(IOContext(io, :compact => true), base_module(V))
end

function show_symmetric_power(io::IO, V::LieAlgebraModule{C}) where {C<:RingElement}
  print(io, "$(get_attribute(V, :power))-th symmetric power of ")
  print(IOContext(io, :compact => true), base_module(V))
end

function show_tensor_power(io::IO, V::LieAlgebraModule{C}) where {C<:RingElement}
  print(io, "$(get_attribute(V, :power))-th tensor power of ")
  print(IOContext(io, :compact => true), base_module(V))
end

function symbols(V::LieAlgebraModule{C}) where {C<:RingElement}
  return V.s
end

function expressify(
  v::LieAlgebraModuleElem{C}, s=symbols(parent(v)); context=nothing
) where {C<:RingElement}
  sum = Expr(:call, :+)
  for (i, c) in enumerate(coefficients(v))
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
  if is_dual(V) && base_module(V) == parent(v)
    return V(coefficients(v))
  end
  @req V == parent(v) "Incompatible modules."
  return v
end

function (V::LieAlgebraModule{C})(
  a::Vector{T}
) where {T<:LieAlgebraModuleElem{C}} where {C<:RingElement}
  if is_direct_sum(V) || is_tensor_product(V)
    @req length(a) == length(base_modules(V)) "Length of vector does not match."
    @req all(i -> parent(a[i]) == base_modules(V)[i], 1:length(a)) "Incompatible modules."
    if is_direct_sum(V)
      return V(vcat([coefficients(x) for x in a]...))
    elseif is_tensor_product(V)
      mat = zero_matrix(base_ring(V), 1, dim(V))
      for (i, inds) in enumerate(get_attribute(V, :ind_map))
        mat[1, i] += prod(a[j].mat[k] for (j, k) in enumerate(inds))
      end
      return LieAlgebraModuleElem{C}(V, mat)
    end
  elseif is_exterior_power(V) || is_symmetric_power(V) || is_tensor_power(V)
    @req length(a) == get_attribute(V, :power) "Length of vector does not match power."
    @req all(x -> parent(x) == base_module(V), a) "Incompatible modules."
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
  else
    throw(ArgumentError("Invalid input."))
  end
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

function Base.:(==)(V1::LieAlgebraModule{C}, V2::LieAlgebraModule{C}) where {C<:RingElement}
  return V1.dim == V2.dim &&
         V1.s == V2.s &&
         V1.L == V2.L &&
         V1.transformation_matrices == V2.transformation_matrices
end

function Base.hash(V::LieAlgebraModule{C}, h::UInt) where {C<:RingElement}
  b = 0x28b0c111e3ff8526 % UInt
  return xor(hash(V.L, hash(V.dim, hash(V.transformation_matrices, hash(V.s, h)))), b)
end

function Base.:(==)(
  v1::LieAlgebraModuleElem{C}, v2::LieAlgebraModuleElem{C}
) where {C<:RingElement}
  check_parent(v1, v2)
  return coefficients(v1) == coefficients(v2)
end

function Base.hash(v::LieAlgebraModuleElem{C}, h::UInt) where {C<:RingElement}
  b = 0x723913014484513a % UInt
  return xor(hash(coefficients(v), hash(parent(v), h)), b)
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

  cx = coefficients(x)

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
  return get_attribute(V, :type, :fallback)::Symbol == :standard_module
end

function is_dual(V::LieAlgebraModule{C}) where {C<:RingElement}
  return get_attribute(V, :type, :fallback)::Symbol == :dual
end

function is_tensor_product(V::LieAlgebraModule{C}) where {C<:RingElement}
  return get_attribute(V, :type, :fallback)::Symbol == :tensor_product
end

function is_direct_sum(V::LieAlgebraModule{C}) where {C<:RingElement}
  return get_attribute(V, :type, :fallback)::Symbol == :direct_sum
end

function is_exterior_power(V::LieAlgebraModule{C}) where {C<:RingElement}
  return get_attribute(V, :type, :fallback)::Symbol == :exterior_power
end

function is_symmetric_power(V::LieAlgebraModule{C}) where {C<:RingElement}
  return get_attribute(V, :type, :fallback)::Symbol == :symmetric_power
end

function is_tensor_power(V::LieAlgebraModule{C}) where {C<:RingElement}
  return get_attribute(V, :type, :fallback)::Symbol == :tensor_power
end

function base_module(V::LieAlgebraModule{C}) where {C<:RingElement}
  @req is_dual(V) || is_exterior_power(V) || is_symmetric_power(V) || is_tensor_power(V) "Not a power module."
  return get_attribute(V, :base_module)::LieAlgebraModule{C}
end

function base_modules(V::LieAlgebraModule{C}) where {C<:RingElement}
  @req is_direct_sum(V) || is_tensor_product(V) "Not a direct sum or tensor product module."
  return get_attribute(V, :base_modules)::Vector{LieAlgebraModule{C}}
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
  check::Bool=true,
) where {C<:RingElement}
  return LieAlgebraModule{C}(L, dimV, transformation_matrices, Symbol.(s); check)
end

function abstract_module(
  L::LieAlgebra{C},
  dimV::Int,
  struct_consts::Matrix{SRow{C}},
  s::Vector{<:VarName}=[Symbol("v_$i") for i in 1:dimV];
  check::Bool=true,
) where {C<:RingElement}
  @req dim(L) == size(struct_consts, 1) "Invalid structure constants dimensions."
  @req dimV == size(struct_consts, 2) "Invalid structure constants dimensions."
  @req dimV == length(s) "Invalid number of basis element names."

  transformation_matrices = [zero_matrix(base_ring(L), dimV, dimV) for _ in 1:dim(L)]
  for i in 1:dim(L), j in 1:dimV
    transformation_matrices[i][:, j] = transpose(dense_row(struct_consts[i, j], dimV))
  end

  return LieAlgebraModule{C}(L, dimV, transformation_matrices, Symbol.(s); check)
end

function highest_weight_module(L::LieAlgebra{C}, weight::Vector{Int}) where {C<:RingElement}
  struct_consts = lie_algebra_highest_weight_module_struct_consts_gap(L, weight)
  dimV = size(struct_consts, 2)
  V = abstract_module(L, dimV, struct_consts; check=false)
  set_attribute!(V, :highest_weight => weight)
  return V
end

function standard_module(L::LinearLieAlgebra{C}) where {C<:RingElement}
  dim_std_V = L.n
  transformation_matrices = matrix_repr_basis(L)
  s = [Symbol("v_$(i)") for i in 1:dim_std_V]
  std_V = LieAlgebraModule{elem_type(base_ring(L))}(
    L, dim_std_V, transformation_matrices, s; check=false
  )
  set_attribute!(std_V, :type => :standard_module, :show => show_standard_module)
  return std_V
end

function dual(V::LieAlgebraModule{C}) where {C<:RingElement}
  L = base_lie_algebra(V)
  dim_dual_V = dim(V)

  transformation_matrices = map(1:dim(L)) do i
    -transpose(transformation_matrix(V, i))
  end

  if is_standard_module(V)
    parentheses = identity
  else
    parentheses = x -> "($x)"
  end
  s = [Symbol("$(parentheses(s))*") for s in symbols(V)]

  pow_V = LieAlgebraModule{C}(L, dim_dual_V, transformation_matrices, s; check=false)
  set_attribute!(pow_V, :type => :dual, :base_module => V, :show => show_dual)
  return pow_V
end

⊕(V::LieAlgebraModule{C}...) where {C<:RingElement} = direct_sum(V...)

function direct_sum(V::LieAlgebraModule{C}...) where {C<:RingElement}
  @req length(V) >= 1 "At least one module must be given."
  L = base_lie_algebra(V[1])
  @req all(x -> base_lie_algebra(x) == L, V) "All modules must have the same base Lie algebra."

  dim_direct_sum_V = sum(dim, V)
  transformation_matrices = map(1:dim(L)) do i
    block_diagonal_matrix([transformation_matrix(Vj, i) for Vj in V])
  end
  parentheses = x -> "($x)"

  if length(V) == 1
    s = symbols(V[1])
  else
    s = [
      Symbol("$s^($j)") for (j, Vj) in enumerate(V) for
      s in (is_standard_module(Vj) ? symbols(Vj) : parentheses.(symbols(Vj)))
    ]
  end

  direct_sum_V = LieAlgebraModule{C}(
    L, dim_direct_sum_V, transformation_matrices, s; check=false
  )
  set_attribute!(
    direct_sum_V,
    :type => :direct_sum,
    :base_modules => collect(V),
    :show => show_direct_sum,
  )
  return direct_sum_V
end

⊗(V::LieAlgebraModule{C}...) where {C<:RingElement} = tensor_product(V...)

function tensor_product(V::LieAlgebraModule{C}...) where {C<:RingElement}
  @req length(V) >= 1 "At least one module must be given."
  L = base_lie_algebra(V[1])
  @req all(x -> base_lie_algebra(x) == L, V) "All modules must have the same base Lie algebra."

  dim_tensor_product_V = prod(dim, V)
  ind_map = collect(reverse.(ProductIterator([1:dim(Vi) for Vi in reverse(V)])))

  transformation_matrices = map(1:dim(L)) do i
    ys = [transformation_matrix(Vj, i) for Vj in V]
    sum(
      reduce(kronecker_product, (j == i ? ys[j] : one(ys[j]) for j in 1:length(V))) for
      i in 1:length(V)
    )
  end

  if length(V) == 1
    s = symbols(V[1])
  else
    parentheses = x -> "($x)"
    s = [
      Symbol(join(s, " ⊗ ")) for s in
      reverse.(
        ProductIterator([
          is_standard_module(Vi) ? symbols(Vi) : parentheses.(symbols(Vi)) for
          Vi in reverse(V)
        ])
      )
    ]
  end

  tensor_product_V = LieAlgebraModule{C}(
    L, dim_tensor_product_V, transformation_matrices, s; check=false
  )
  set_attribute!(
    tensor_product_V,
    :type => :tensor_product,
    :base_modules => collect(V),
    :ind_map => ind_map,
    :show => show_tensor_product,
  )
  return tensor_product_V
end

function exterior_power(V::LieAlgebraModule{C}, k::Int) where {C<:RingElement}
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

  pow_V = LieAlgebraModule{C}(L, dim_pow_V, transformation_matrices, s; check=false)
  set_attribute!(
    pow_V,
    :type => :exterior_power,
    :power => k,
    :base_module => V,
    :ind_map => ind_map,
    :show => show_exterior_power,
  )
  return pow_V
end

function symmetric_power(V::LieAlgebraModule{C}, k::Int) where {C<:RingElement}
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

  pow_V = LieAlgebraModule{C}(L, dim_pow_V, transformation_matrices, s; check=false)
  set_attribute!(
    pow_V,
    :type => :symmetric_power,
    :power => k,
    :base_module => V,
    :ind_map => ind_map,
    :show => show_symmetric_power,
  )
  return pow_V
end

function tensor_power(V::LieAlgebraModule{C}, k::Int) where {C<:RingElement}
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

  pow_V = LieAlgebraModule{C}(L, dim_pow_V, transformation_matrices, s; check=false)
  set_attribute!(
    pow_V,
    :type => :tensor_power,
    :power => k,
    :base_module => V,
    :ind_map => ind_map,
    :show => show_tensor_power,
  )
  return pow_V
end
