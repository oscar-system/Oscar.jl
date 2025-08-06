###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(
  ::Type{LieAlgebraModuleElem{C,LieT}}
) where {C<:FieldElem,LieT<:LieAlgebraElem{C}} = LieAlgebraModule{C,LieT}

elem_type(::Type{LieAlgebraModule{C,LieT}}) where {C<:FieldElem,LieT<:LieAlgebraElem{C}} =
  LieAlgebraModuleElem{C,LieT}

parent(v::LieAlgebraModuleElem) = v.parent

coefficient_ring(V::LieAlgebraModule) = coefficient_ring(base_lie_algebra(V))

coefficient_ring(v::LieAlgebraModuleElem) = coefficient_ring(parent(v))

@doc raw"""
    base_lie_algebra(V::LieAlgebraModule{C}) -> LieAlgebra{C}

Return the Lie algebra `V` is a module over.
"""
base_lie_algebra(V::LieAlgebraModule{C,LieT}) where {C<:FieldElem,LieT<:LieAlgebraElem{C}} =
  V.L::parent_type(LieT)

number_of_generators(L::LieAlgebraModule) = dim(L)

gens(L::LieAlgebraModule) = basis(L)

gen(L::LieAlgebraModule, i::Int) = basis(L, i)

@doc raw"""
    dim(V::LieAlgebraModule{C}) -> Int

Return the dimension of the Lie algebra module `V`.
"""
dim(V::LieAlgebraModule) = vector_space_dim(V)
vector_space_dim(V::LieAlgebraModule) = V.dim

@doc raw"""
    basis(V::LieAlgebraModule{C}) -> Vector{LieAlgebraModuleElem{C}}

Return a basis of the Lie algebra module `V`.
"""
basis(V::LieAlgebraModule) = [basis(V, i)::elem_type(V) for i in 1:dim(V)]

@doc raw"""
    basis(V::LieAlgebraModule{C}, i::Int) -> LieAlgebraModuleElem{C}

Return the `i`-th basis element of the Lie algebra module `V`.
"""
function basis(V::LieAlgebraModule, i::Int)
  @req 1 <= i <= dim(V) "Index out of bounds."
  R = coefficient_ring(V)
  v = zero_matrix(R, 1, dim(V))
  v[1, i] = one(R)
  return V(v)
end

@doc raw"""
    zero(V::LieAlgebraModule{C}) -> LieAlgebraModuleElem{C}

Return the zero element of the Lie algebra module `V`.
"""
function zero(V::LieAlgebraModule)
  mat = zero_matrix(coefficient_ring(V), 1, dim(V))
  return elem_type(V)(V, mat)
end

@doc raw"""
    iszero(v::LieAlgebraModuleElem{C}) -> Bool

Check whether the Lie algebra module element `v` is zero.
"""
function iszero(v::LieAlgebraModuleElem)
  return iszero(coefficients(v))
end

@inline function _matrix(v::LieAlgebraModuleElem{C}) where {C<:FieldElem}
  return (v.mat)::dense_matrix_type(C)
end

@doc raw"""
    coefficients(v::LieAlgebraModuleElem{C}) -> Vector{C}

Return the coefficients of `v` with respect to [`basis(::LieAlgebraModule)`](@ref).
"""
function coefficients(v::LieAlgebraModuleElem)
  return collect(_matrix(v))[1, :]
end

@doc raw"""
    coeff(v::LieAlgebraModuleElem{C}, i::Int) -> C

Return the `i`-th coefficient of `v` with respect to [`basis(::LieAlgebraModule)`](@ref).
"""
function coeff(v::LieAlgebraModuleElem, i::Int)
  return _matrix(v)[1, i]
end

@doc raw"""
    getindex(v::LieAlgebraModuleElem{C}, i::Int) -> C

Return the `i`-th coefficient of `v` with respect to [`basis(::LieAlgebraModule)`](@ref).
"""
function getindex(v::LieAlgebraModuleElem, i::Int)
  return coeff(v, i)
end

function Base.deepcopy_internal(v::LieAlgebraModuleElem, dict::IdDict)
  return parent(v)(deepcopy_internal(_matrix(v), dict))
end

function check_parent(
  v1::LieAlgebraModuleElem{C}, v2::LieAlgebraModuleElem{C}
) where {C<:FieldElem}
  parent(v1) !== parent(v2) && error("Incompatible modules.")
end

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, mime::MIME"text/plain", V::LieAlgebraModule)
  @show_name(io, V)
  @show_special(io, mime, V)
  io = pretty(io)
  println(io, _module_type_to_string(V))
  println(io, Indent(), "of dimension $(dim(V))")
  if _is_dual(V)[1] ||
    _is_direct_sum(V)[1] ||
    _is_tensor_product(V)[1] ||
    _is_exterior_power(V)[1] ||
    _is_symmetric_power(V)[1] ||
    _is_tensor_power(V)[1]
    _show_inner(io, V)
  end
  print(io, Dedent())
  print(io, "over ")
  print(io, Lowercase(), base_lie_algebra(V))
end

function _show_inner(io::IO, V::LieAlgebraModule)
  if _is_standard_module(V)
    println(io, "standard module")
  elseif ((fl, W) = _is_dual(V); fl)
    println(io, "dual of ", Lowercase())
    print(io, Indent())
    _show_inner(io, W)
    print(io, Dedent())
  elseif ((fl, Ws) = _is_direct_sum(V); fl)
    println(io, "direct sum with direct summands")
    print(io, Indent())
    for W in Ws
      _show_inner(io, W)
    end
    print(io, Dedent())
  elseif ((fl, Ws) = _is_tensor_product(V); fl)
    println(io, "tensor product with tensor factors")
    print(io, Indent())
    for W in Ws
      _show_inner(io, W)
    end
    print(io, Dedent())
  elseif ((fl, W, k) = _is_exterior_power(V); fl)
    println(io, "$(ordinal_number_string(k)) exterior power of")
    print(io, Indent())
    _show_inner(io, W)
    print(io, Dedent())
  elseif ((fl, W, k) = _is_symmetric_power(V); fl)
    println(io, "$(ordinal_number_string(k)) symmetric power of")
    print(io, Indent())
    _show_inner(io, W)
    print(io, Dedent())
  elseif ((fl, W, k) = _is_tensor_power(V); fl)
    println(io, "$(ordinal_number_string(k)) tensor power of")
    print(io, Indent())
    _show_inner(io, W)
    print(io, Dedent())
  else
    println(io, "abstract module")
  end
end

function Base.show(io::IO, V::LieAlgebraModule)
  @show_name(io, V)
  @show_special(io, V)
  if is_terse(io)
    print(io, _module_type_to_string(V))
  else
    io = pretty(io)
    print(io, _module_type_to_string(V), " of dimension $(dim(V)) over ", Lowercase())
    print(terse(io), base_lie_algebra(V))
  end
end

function _module_type_to_string(V::LieAlgebraModule)
  if _is_standard_module(V)
    return "Standard module"
  elseif _is_dual(V)[1]
    return "Dual module"
  elseif _is_direct_sum(V)[1]
    return "Direct sum module"
  elseif _is_tensor_product(V)[1]
    return "Tensor product module"
  elseif _is_exterior_power(V)[1]
    return "Exterior power module"
  elseif _is_symmetric_power(V)[1]
    return "Symmetric power module"
  elseif _is_tensor_power(V)[1]
    return "Tensor power module"
  else
    return "Abstract Lie algebra module"
  end
end

@doc raw"""
    symbols(V::LieAlgebraModule{C}) -> Vector{Symbol}

Return the symbols used for printing basis elements of the Lie algebra module `V`.
"""
function symbols(V::LieAlgebraModule)
  return V.s
end

function expressify(v::LieAlgebraModuleElem, s=symbols(parent(v)); context=nothing)
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

@doc raw"""
    (V::LieAlgebraModule{C})() -> LieAlgebraModuleElem{C}

Return the zero element of the Lie algebra module `V`.
"""
function (V::LieAlgebraModule)()
  return zero(V)
end

@doc raw"""
    (V::LieAlgebraModule{C})(v::AbstractVector{Int}) -> LieAlgebraModuleElem{C}

Return the element of `V` with coefficient vector `v`.
Fails, if `Int` cannot be coerced into the base ring of `L`.
"""
function (V::LieAlgebraModule)(v::AbstractVector{Int})
  return V(coefficient_ring(V).(v))
end

@doc raw"""
    (V::LieAlgebraModule{C})(v::AbstractVector{C}) -> LieAlgebraModuleElem{C}

Return the element of `V` with coefficient vector `v`.
"""
function (V::LieAlgebraModule{C})(v::AbstractVector{C}) where {C<:FieldElem}
  @req length(v) == dim(V) "Length of vector does not match dimension."
  mat = matrix(coefficient_ring(V), 1, length(v), v)
  return elem_type(V)(V, mat)
end

@doc raw"""
    (V::LieAlgebraModule{C})(mat::MatElem{C}) -> LieAlgebraModuleElem{C}

Return the element of `V` with coefficient vector equivalent to
the $1 \times \dim(L)$ matrix `mat`.
"""
function (V::LieAlgebraModule{C})(v::MatElem{C}) where {C<:FieldElem}
  @req ncols(v) == dim(V) "Length of vector does not match dimension"
  @req nrows(v) == 1 "Not a vector in module constructor"
  return elem_type(V)(V, v)
end

@doc raw"""
    (V::LieAlgebraModule{C})(v::SRow{C}) -> LieAlgebraModuleElem{C}

Return the element of `V` with coefficient vector `v`.
"""
function (V::LieAlgebraModule{C})(v::SRow{C}) where {C<:FieldElem}
  mat = dense_row(v, dim(V))
  return elem_type(V)(V, mat)
end

@doc raw"""
    (V::LieAlgebraModule{C})(v::LieAlgebraModuleElem{C}) -> LieAlgebraModuleElem{C}

Return `v`. Fails, in general, if `v` is not an element of `V`.

If `V` is the dual module of the parent of `v`, return the dual of `v`.
"""
function (V::LieAlgebraModule{C})(v::LieAlgebraModuleElem{C}) where {C<:FieldElem}
  # handled by more general function below
  return V((v,))::elem_type(V)
end

@doc raw"""
    (V::LieAlgebraModule{C})(a::Vector{T}) where {T<:LieAlgebraModuleElem{C}}) -> LieAlgebraModuleElem{C}
    (V::LieAlgebraModule{C})(a::NTuple{T}) where {T<:LieAlgebraModuleElem{C}}) -> LieAlgebraModuleElem{C}
    (V::LieAlgebraModule{C})(a::T...) where {T<:LieAlgebraModuleElem{C}}) -> LieAlgebraModuleElem{C}

If `V` is a direct sum, return its element, where the $i$-th component is equal to `a[i]`.
If `V` is a tensor product, return the tensor product of the `a[i]`.
If `V` is a exterior (symmetric, tensor) power, return the wedge product
(product, tensor product) of the `a[i]`.

Requires that `a` has a suitable length, and that the `a[i]` are elements of the correct modules,
where _correct_ depends on the case above.
"""
function (V::LieAlgebraModule{C})(
  a::Vector{T}
) where {C<:FieldElem,T<:LieAlgebraModuleElem{C}}
  @req _is_allowed_input_length(V, length(a)) "Invalid input length." # Check here to not compile unnecessary dispatches on tuples
  return V(Tuple(a))::elem_type(V)
end

function (V::LieAlgebraModule{C})(
  a::Tuple{Vararg{LieAlgebraModuleElem{C}}}
) where {C<:FieldElem}
  if length(a) == 1 && V === parent(a[1])
    return a[1]
  elseif ((fl, W) = _is_dual(V); fl)
    @req length(a) == 1 "Invalid input length."
    @req W === parent(a) "Incompatible modules."
    return V(coefficients(a))
  elseif ((fl, Vs) = _is_direct_sum(V); fl)
    @req length(a) == length(Vs) "Invalid input length."
    @req all(i -> parent(a[i]) === Vs[i], 1:length(a)) "Incompatible modules."
    return sum(inji(ai) for (ai, inji) in zip(a, canonical_injections(V)); init=zero(V))
  elseif _is_tensor_product(V)[1]
    pure = get_attribute(V, :tensor_pure_function)
    return pure(a)::elem_type(V)
  elseif _is_exterior_power(V)[1]
    pure = get_attribute(V, :wedge_pure_function)
    return pure(a)::elem_type(V)
  elseif _is_symmetric_power(V)[1]
    pure = get_attribute(V, :mult_pure_function)
    return pure(a)::elem_type(V)
  elseif _is_tensor_power(V)[1]
    pure = get_attribute(V, :tensor_pure_function)
    return pure(a)::elem_type(V)
  else
    throw(ArgumentError("Invalid input."))
  end
end

function (V::LieAlgebraModule{C})(a::LieAlgebraModuleElem{C}...) where {C<:FieldElem}
  return V(a)
end

function _is_allowed_input_length(V::LieAlgebraModule, a::Int)
  if ((fl, Vs) = _is_direct_sum(V); fl)
    return a == length(Vs)
  elseif ((fl, Vs) = _is_tensor_product(V); fl)
    return a == length(Vs)
  elseif ((fl, W, k) = _is_exterior_power(V); fl)
    return a == k
  elseif ((fl, W, k) = _is_symmetric_power(V); fl)
    return a == k
  elseif ((fl, W, k) = _is_tensor_power(V); fl)
    return a == k
  else
    throw(ArgumentError("Invalid input."))
  end
end

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

function Base.:-(v::LieAlgebraModuleElem{C}) where {C<:FieldElem}
  return parent(v)(-_matrix(v))
end

function Base.:+(
  v1::LieAlgebraModuleElem{C}, v2::LieAlgebraModuleElem{C}
) where {C<:FieldElem}
  check_parent(v1, v2)
  return parent(v1)(_matrix(v1) + _matrix(v2))
end

function Base.:-(
  v1::LieAlgebraModuleElem{C}, v2::LieAlgebraModuleElem{C}
) where {C<:FieldElem}
  check_parent(v1, v2)
  return parent(v1)(_matrix(v1) - _matrix(v2))
end

function Base.:*(v::LieAlgebraModuleElem{C}, c::C) where {C<:FieldElem}
  coefficient_ring(v) != parent(c) && error("Incompatible rings.")
  return parent(v)(_matrix(v) * c)
end

function Base.:*(v::LieAlgebraModuleElem, c::U) where {U<:RationalUnion}
  return parent(v)(_matrix(v) * c)
end

function Base.:*(c::C, v::LieAlgebraModuleElem{C}) where {C<:FieldElem}
  coefficient_ring(v) != parent(c) && error("Incompatible rings.")
  return parent(v)(c * _matrix(v))
end

function Base.:*(c::U, v::LieAlgebraModuleElem) where {U<:RationalUnion}
  return parent(v)(c * _matrix(v))
end

###############################################################################
#
#   Comparison functions
#
###############################################################################

function Base.:(==)(V1::LieAlgebraModule{C}, V2::LieAlgebraModule{C}) where {C<:FieldElem}
  return dim(V1) == dim(V2) &&
         symbols(V1) == symbols(V2) &&
         base_lie_algebra(V1) == base_lie_algebra(V2) &&
         transformation_matrices(V1) == transformation_matrices(V2)
end

function Base.hash(V::LieAlgebraModule, h::UInt)
  b = 0x28b0c111e3ff8526 % UInt
  h = hash(dim(V), h)
  h = hash(symbols(V), h)
  h = hash(base_lie_algebra(V), h)
  h = hash(transformation_matrices(V), h)
  return xor(h, b)
end

function Base.:(==)(
  v1::LieAlgebraModuleElem{C}, v2::LieAlgebraModuleElem{C}
) where {C<:FieldElem}
  check_parent(v1, v2)
  return coefficients(v1) == coefficients(v2)
end

function Base.hash(v::LieAlgebraModuleElem, h::UInt)
  b = 0x723913014484513a % UInt
  h = hash(parent(v), h)
  h = hash(coefficients(v), h)
  return xor(h, b)
end

###############################################################################
#
#   Module action
#
###############################################################################

@doc raw"""
    action(x::LieAlgebraElem{C}, v::LieAlgebraModuleElem{C}) -> LieAlgebraModuleElem{C}
    *(x::LieAlgebraElem{C}, v::LieAlgebraModuleElem{C}) -> LieAlgebraModuleElem{C}

Apply the action of `x` on `v`.
"""
function Base.:*(x::LieAlgebraElem{C}, v::LieAlgebraModuleElem{C}) where {C<:FieldElem}
  return action(x, v)
end

function action(x::LieAlgebraElem{C}, v::LieAlgebraModuleElem{C}) where {C<:FieldElem}
  @req parent(x) === base_lie_algebra(parent(v)) "Incompatible Lie algebras."

  cx = coefficients(x)
  V = parent(v)
  return V(
    sum(
      cx[i] * _matrix(v) * transformation_matrix(V, i) for
      i in 1:dim(parent(x)) if !iszero(cx[i]);
      init=zero_matrix(coefficient_ring(V), 1, dim(V))::dense_matrix_type(C),
    ),
  )
end

function transformation_matrices(V::LieAlgebraModule{C}) where {C<:FieldElem}
  return (V.transformation_matrices)::Vector{dense_matrix_type(C)}
end

function transformation_matrix(V::LieAlgebraModule{C}, i::Int) where {C<:FieldElem}
  return transformation_matrices(V)[i]
end

###############################################################################
#
#   Attribute getters
#
###############################################################################

@doc raw"""
    _is_standard_module(V::LieAlgebraModule{C}) -> Bool

Check whether `V` has been constructed as a standard module.
"""
function _is_standard_module(V::LieAlgebraModule)
  if has_attribute(V, :is_standard_module)
    @assert get_attribute(V, :is_standard_module)::Bool === true
    return true
  end
  return false
end

@doc raw"""
    _is_dual(V::LieAlgebraModule{C}) -> Bool, LieAlgebraModule{C}

Check whether `V` has been constructed as a dual module.
If it has, return `true` and the base module.
If not, return `false` as the first return value, and an arbitrary value for the second.
"""
function _is_dual(V::LieAlgebraModule)
  if has_attribute(V, :is_dual)
    W = get_attribute(V, :is_dual)::typeof(V)
    return (true, W)
  end
  return (false, V)
end

@doc raw"""
    _is_direct_sum(V::LieAlgebraModule{C}) -> Bool, Vector{LieAlgebraModule{C}}

Check whether `V` has been constructed as a direct sum of modules.
If it has, return `true` and the summands.
If not, return `false` as the first return value, and an empty vector for the second.
"""
function _is_direct_sum(V::LieAlgebraModule)
  if has_attribute(V, :is_direct_sum)
    summands = get_attribute(V, :is_direct_sum)::Vector{typeof(V)}
    return (true, summands)
  end
  return (false, typeof(V)[])
end

@doc raw"""
    _is_tensor_product(V::LieAlgebraModule{C}) -> Bool, Vector{LieAlgebraModule{C}}

Check whether `V` has been constructed as a tensor product of modules.
If it has, return `true` and the tensor factors.
If not, return `false` as the first return value, and an empty vector for the second.
"""
function _is_tensor_product(V::LieAlgebraModule)
  if has_attribute(V, :is_tensor_product)
    factors = get_attribute(V, :is_tensor_product)::Vector{typeof(V)}
    return (true, factors)
  end
  return (false, typeof(V)[])
end

@doc raw"""
    _is_exterior_power(V::LieAlgebraModule{C}) -> Bool, LieAlgebraModule{C}, Int

Check whether `V` has been constructed as an exterior power of a module.
If it has, return `true`, the base module, and the power.
If not, return `false` as the first return value, and arbitrary values for the other two.
"""
function _is_exterior_power(V::LieAlgebraModule)
  if has_attribute(V, :is_exterior_power)
    W, k = get_attribute(V, :is_exterior_power)::Tuple{typeof(V),Int}
    return (true, W, k)
  end
  return (false, V, 0)
end

@doc raw"""
    _is_symmetric_power(V::LieAlgebraModule{C}) -> Bool, LieAlgebraModule{C}, Int

Check whether `V` has been constructed as an symmetric power of a module.
If it has, return `true`, the base module, and the power.
If not, return `false` as the first return value, and arbitrary values for the other two.
"""
function _is_symmetric_power(V::LieAlgebraModule)
  if has_attribute(V, :is_symmetric_power)
    W, k = get_attribute(V, :is_symmetric_power)::Tuple{typeof(V),Int}
    return (true, W, k)
  end
  return (false, V, 0)
end

@doc raw"""
    _is_tensor_power(V::LieAlgebraModule{C}) -> Bool, LieAlgebraModule{C}, Int

Check whether `V` has been constructed as a tensor power of a module.
If it has, return `true`, the base module, and the power.
If not, return `false` as the first return value, and arbitrary values for the other two.
"""
function _is_tensor_power(V::LieAlgebraModule)
  if has_attribute(V, :is_tensor_power)
    W, k = get_attribute(V, :is_tensor_power)::Tuple{typeof(V),Int}
    return (true, W, k)
  end
  return (false, V, 0)
end

###############################################################################
#
#   Constructor
#
###############################################################################

@doc raw"""
    abstract_module(L::LieAlgebra{C}, dimV::Int, transformation_matrices::Vector{<:MatElem{C}}, s::Vector{<:VarName}; check::Bool) -> LieAlgebraModule{C}

Construct the the Lie algebra module over `L` of dimension `dimV` given by
`transformation_matrices` and with basis element names `s`.

* `transformation_matrices`: The action of the $i$-th basis element of `L`
  on some element $v$ of the constructed module is given by right multiplication 
  of the matrix `transformation_matrices[i]` to the coefficient vector of $v$.
* `s`: A vector of basis element names. This is 
  `[Symbol("v_$i") for i in 1:dimV]` by default.
* `check`: If `true`, check that the structure constants are anti-symmetric and
  satisfy the Jacobi identity. This is `true` by default.
"""
function abstract_module(
  L::LieAlgebra{C},
  dimV::Int,
  transformation_matrices::Vector{<:MatElem{C}},
  s::Vector{<:VarName}=[Symbol("v_$i") for i in 1:dimV];
  check::Bool=true,
) where {C<:FieldElem}
  return LieAlgebraModule{C}(L, dimV, transformation_matrices, Symbol.(s); check)
end

@doc raw"""
    abstract_module(L::LieAlgebra{C}, dimV::Int, struct_consts::Matrix{sparse_row_type{C}}, s::Vector{<:VarName}; check::Bool) -> LieAlgebraModule{C}

Construct the the Lie algebra module over `L` of dimension `dimV` given by
structure constants `struct_consts` and with basis element names `s`.

The action on the newly constructed Lie algebra module `V` is determined by the structure
constants in `struct_consts` as follows: let $x_i$ denote the $i$-th standard basis vector
of `L`, and $v_i$ the $i$-th standard basis vector of `V`.
Then the entry `struct_consts[i,j][k]` is a scalar $a_{i,j,k}$
such that $x_i * v_j = \sum_k a_{i,j,k} v_k$.

* `s`: A vector of basis element names. This is 
  `[Symbol("v_$i") for i in 1:dimV]` by default.
* `check`: If `true`, check that the structure constants are anti-symmetric and
  satisfy the Jacobi identity. This is `true` by default.
"""
function abstract_module(
  L::LieAlgebra{C},
  dimV::Int,
  struct_consts::Matrix{<:SRow{C}},
  s::Vector{<:VarName}=[Symbol("v_$i") for i in 1:dimV];
  check::Bool=true,
) where {C<:FieldElem}
  @req dim(L) == size(struct_consts, 1) "Invalid structure constants dimensions."
  @req dimV == size(struct_consts, 2) "Invalid structure constants dimensions."
  @req dimV == length(s) "Invalid number of basis element names."

  transformation_matrices = [zero_matrix(coefficient_ring(L), dimV, dimV) for _ in 1:dim(L)]
  for i in 1:dim(L), j in 1:dimV
    transformation_matrices[i][j, :] = dense_row(
      struct_consts[i, j]::sparse_row_type(C), dimV
    )
  end

  return LieAlgebraModule{C}(L, dimV, transformation_matrices, Symbol.(s); check)
end

######################## Trivial module ########################

@doc raw"""
    trivial_module(L::LieAlgebra{C}, d=1) -> LieAlgebraModule{C}

Construct the `d`-dimensional module of the Lie algebra `L` with trivial action.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 3);

julia> trivial_module(L)
Abstract Lie algebra module
  of dimension 1
over special linear Lie algebra of degree 3 over QQ
```
"""
function trivial_module(L::LieAlgebra, d::IntegerUnion=1)
  @req d >= 0 "Dimension must be non-negative."
  dim_triv_V = Int(d)
  transformation_matrices = [
    zero_matrix(coefficient_ring(L), dim_triv_V, dim_triv_V) for _ in 1:dim(L)
  ]
  s = [Symbol("v_$(i)") for i in 1:dim_triv_V]
  triv_V = LieAlgebraModule{elem_type(coefficient_ring(L))}(
    L, dim_triv_V, transformation_matrices, s; check=false
  )
  return triv_V
end

######################## Standard module ########################

# cache storage
@attr LieAlgebraModule{elem_type(coefficient_ring(L))} function _standard_module(
  L::LinearLieAlgebra
)
  return standard_module(L; cached=false)
end

@doc raw"""
    standard_module(L::LinearLieAlgebra{C}; cached::Bool=true) -> LieAlgebraModule{C}

Construct the standard module of the linear Lie algebra `L`.
If `L` is a Lie subalgebra of $\mathfrak{gl}_n(R)$, then the standard module
is $R^n$ with the action of $L$ given by left multiplication.

!!! note
    This uses the left action of $L$, and converts that to internally use the equivalent
    right action.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 3);

julia> standard_module(L)
Standard module
  of dimension 3
over special linear Lie algebra of degree 3 over QQ
```
"""
function standard_module(L::LinearLieAlgebra; cached::Bool=true)
  cached && return _standard_module(L)

  dim_std_V = L.n
  transformation_matrices = transpose.(matrix_repr_basis(L))
  s = [Symbol("v_$(i)") for i in 1:dim_std_V]
  std_V = LieAlgebraModule{elem_type(coefficient_ring(L))}(
    L, dim_std_V, transformation_matrices, s; check=false
  )
  set_attribute!(std_V, :is_standard_module => true)

  return std_V
end

######################## Duals ########################

# cache storage
@attr typeof(V) function _dual(V::LieAlgebraModule)
  return dual(V; cached=false)
end

@doc raw"""
    dual(V::LieAlgebraModule{C}; cached::Bool=true) -> LieAlgebraModule{C}

Construct the dual module of `V`.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 3);

julia> V = exterior_power(standard_module(L), 2)[1]; # some module

julia> dual(V)
Dual module
  of dimension 3
  dual of
    2nd exterior power of
      standard module
over special linear Lie algebra of degree 3 over QQ
```
"""
function dual(V::LieAlgebraModule{C}; cached::Bool=true) where {C<:FieldElem}
  cached && return _dual(V)

  L = base_lie_algebra(V)
  dim_dual_V = dim(V)

  transformation_matrices = map(1:dim(L)) do i
    -transpose(transformation_matrix(V, i))
  end

  s = if _is_standard_module(V)
    [Symbol("$(s)*") for s in symbols(V)]
  else
    [Symbol("($(s))*") for s in symbols(V)]
  end

  pow_V = LieAlgebraModule{C}(L, dim_dual_V, transformation_matrices, s; check=false)
  set_attribute!(pow_V, :is_dual => V)

  return pow_V
end

Base.:^(V::LieAlgebraModule, ::typeof(Base.:*)) = dual(V)

######################## Direct sums ########################

# TODO: add caching
@doc raw"""
    direct_sum(V::LieAlgebraModule{C}...) -> LieAlgebraModule{C}
    ⊕(V::LieAlgebraModule{C}...) -> LieAlgebraModule{C}

Construct the direct sum of the modules `V...`.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 3);

julia> V1 = exterior_power(standard_module(L), 2)[1]; # some module

julia> V2 = symmetric_power(standard_module(L), 3)[1]; # some module

julia> direct_sum(V1, V2)
Direct sum module
  of dimension 13
  direct sum with direct summands
    2nd exterior power of
      standard module
    3rd symmetric power of
      standard module
over special linear Lie algebra of degree 3 over QQ
```
"""
function direct_sum(V::LieAlgebraModule{C}, Vs::LieAlgebraModule{C}...) where {C<:FieldElem}
  Vs = [V; Vs...]

  L = base_lie_algebra(Vs[1])
  @req all(x -> base_lie_algebra(x) === L, Vs) "All modules must have the same base Lie algebra."

  dim_direct_sum_V = sum(dim, Vs; init=0)
  transformation_matrices = map(1:dim(L)) do i
    block_diagonal_matrix([transformation_matrix(Vj, i) for Vj in Vs])
  end

  s = if length(Vs) == 1
    symbols(Vs[1])
  else
    [
      Symbol("$s^($j)") for (j, Vj) in enumerate(Vs) for
      s in (_is_standard_module(Vj) ? symbols(Vj) : (x -> "($x)").(symbols(Vj)))
    ]
  end

  direct_sum_V = LieAlgebraModule{C}(
    L, dim_direct_sum_V, transformation_matrices, s; check=false
  )
  set_attribute!(direct_sum_V, :is_direct_sum => Vs)
  return direct_sum_V
end

⊕(V::LieAlgebraModule{C}, Vs::LieAlgebraModule{C}...) where {C<:FieldElem} =
  direct_sum(V, Vs...)

######################## Tensor products ########################

# TODO: add caching
@doc raw"""
    tensor_product(Vs::LieAlgebraModule{C}...) -> LieAlgebraModule{C}
    ⊗(Vs::LieAlgebraModule{C}...) -> LieAlgebraModule{C}

Given modules $V_1,\dots,V_n$ over the same Lie algebra $L$,
construct their tensor product $V_1 \otimes \cdots \otimes V_n$.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 3);

julia> V1 = exterior_power(standard_module(L), 2)[1]; # some module

julia> V2 = symmetric_power(standard_module(L), 3)[1]; # some module

julia> tensor_product(V1, V2)
Tensor product module
  of dimension 30
  tensor product with tensor factors
    2nd exterior power of
      standard module
    3rd symmetric power of
      standard module
over special linear Lie algebra of degree 3 over QQ
```
"""
function tensor_product(
  V::LieAlgebraModule{C}, Vs::LieAlgebraModule{C}...
) where {C<:FieldElem}
  Vs = [V; Vs...]
  L = base_lie_algebra(Vs[1])
  @req all(x -> base_lie_algebra(x) === L, Vs) "All modules must have the same base Lie algebra."
  R = coefficient_ring(Vs[1])

  dim_tensor_product_V = prod(dim, Vs; init=1)
  ind_map = collect(reverse.(ProductIterator([1:dim(Vi) for Vi in reverse(Vs)])))

  transformation_matrices = map(1:dim(L)) do i
    ys = [transformation_matrix(Vj, i) for Vj in Vs]
    sum(
      kronecker_product(
        kronecker_product(identity_matrix(R, prod(dim, Vs[1:(j - 1)]; init=1)), ys[j]),
        identity_matrix(R, prod(dim, Vs[(j + 1):end]; init=1)),
      ) for j in 1:length(Vs)
    )
  end

  s = if length(Vs) == 1
    symbols(Vs[1])
  else
    [
      Symbol(join(s, (is_unicode_allowed() ? "⊗" : "(x)"))) for s in
      reverse.(
        ProductIterator([
          _is_standard_module(Vi) ? symbols(Vi) : (x -> "($x)").(symbols(Vi)) for
          Vi in reverse(Vs)
        ])
      )
    ]
  end

  tensor_product_V = LieAlgebraModule{C}(
    L, dim_tensor_product_V, transformation_matrices, s; check=false
  )

  function my_mult(as::Tuple{Vararg{LieAlgebraModuleElem{C}}})
    @req length(as) == length(Vs) "Length of vector does not match."
    @req all(i -> parent(as[i]) === Vs[i], 1:length(as)) "Incompatible modules."
    mat = zero_matrix(R, 1, dim_tensor_product_V)
    for (i, inds) in enumerate(ind_map)
      mat[1, i] += prod(as[j].mat[k]::elem_type(R) for (j, k) in enumerate(inds))
    end
    return tensor_product_V(mat)
  end
  function my_mult(as::LieAlgebraModuleElem{C}...)
    return my_mult(as)::LieAlgebraModuleElem{C}
  end

  function my_decomp(a::LieAlgebraModuleElem{C})
    @req parent(a) === tensor_product_V "Incompatible modules."
    if iszero(a)
      return Tuple(zero(Vi) for Vi in Vs)
    end
    nz = findall(!iszero, coefficients(a))
    @req length(nz) == 1 "Non-pure tensor product element."
    @req isone(coeff(a, nz[1])) "Non-pure tensor product element."
    inds = ind_map[nz[1]]
    return Tuple(basis(Vi, indi) for (Vi, indi) in zip(Vs, inds))
  end

  mult_map = MapFromFunc(
    Hecke.TupleParent(Tuple([zero(Vi) for Vi in Vs])), tensor_product_V, my_mult, my_decomp
  )
  inv_mult_map = MapFromFunc(tensor_product_V, domain(mult_map), my_decomp, my_mult)

  set_attribute!(
    tensor_product_V,
    :is_tensor_product => Vs,
    :tensor_pure_function => mult_map,
    :tensor_generator_decompose_function => inv_mult_map,
  )

  return tensor_product_V
end

⊗(V::LieAlgebraModule{C}, Vs::LieAlgebraModule{C}...) where {C<:FieldElem} =
  tensor_product(V, Vs...)

######################## Exterior powers ########################

# cache storage
@attr Dict{Int,Tuple{typeof(V),MapFromFunc}} function _exterior_powers(V::LieAlgebraModule)
  return Dict{Int,Tuple{typeof(V),MapFromFunc}}()
end

@doc raw"""
    exterior_power(V::LieAlgebraModule{C}, k::Int; cached::Bool=true) -> LieAlgebraModule{C}, Map

Construct the `k`-th exterior power $\bigwedge^k (V)$ of the module `V`,
together with a map that computes the wedge product of `k` elements of `V`.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 2);

julia> V = symmetric_power(standard_module(L), 2)[1]; # some module

julia> E, map = exterior_power(V, 2)
(Exterior power module of dimension 3 over L, Map: parent of tuples of type Tuple{LieAlgebraModuleElem{QQFieldElem, LinearLieAlgebraElem{QQFieldElem}}, LieAlgebraModuleElem{QQFieldElem, LinearLieAlgebraElem{QQFieldElem}}} -> E)

julia> E
Exterior power module
  of dimension 3
  2nd exterior power of
    2nd symmetric power of
      standard module
over special linear Lie algebra of degree 2 over QQ

julia> basis(E)
3-element Vector{LieAlgebraModuleElem{QQFieldElem, LinearLieAlgebraElem{QQFieldElem}}}:
 (v_1^2)^(v_1*v_2)
 (v_1^2)^(v_2^2)
 (v_1*v_2)^(v_2^2)
```
"""
function exterior_power(
  V::LieAlgebraModule{C}, k::Int; cached::Bool=true
) where {C<:FieldElem}
  @req k >= 0 "Non-negative exponent needed"
  @req characteristic(coefficient_ring(V)) == 0 "Characteristic must be zero."

  if cached
    powers = _exterior_powers(V)
    haskey(powers, k) && return powers[k]
  end

  L = base_lie_algebra(V)
  R = coefficient_ring(V)
  dim_E = binomial(dim(V), k)
  ind_map = collect(combinations(1:dim(V), k))

  T, _ = tensor_power(V, k; cached)
  E_to_T_mat = zero_matrix(R, dim_E, dim(T))
  T_to_E_mat = zero_matrix(R, dim(T), dim_E)
  for (i, _inds) in enumerate(ind_map), (inds, sgn) in permutations_with_sign(_inds)
    j = 1 + sum((ind - 1) * dim(V)^(k - 1) for (k, ind) in enumerate(reverse(inds)); init=0)
    E_to_T_mat[i, j] = R(sgn)//factorial(k)
    T_to_E_mat[j, i] = R(sgn)
  end
  transformation_matrices = map(1:dim(L)) do i
    E_to_T_mat * transformation_matrix(T, i) * T_to_E_mat
  end

  s = if k == 0
    [Symbol("1")]
  elseif k == 1
    symbols(V)
  elseif _is_standard_module(V)
    [Symbol(join(s, (is_unicode_allowed() ? "∧" : "^"))) for s in combinations(symbols(V), k)]
  else
    [
      Symbol(join((x -> "($x)").(s), (is_unicode_allowed() ? "∧" : "^"))) for
      s in combinations(symbols(V), k)
    ]
  end

  E = LieAlgebraModule{C}(L, dim_E, transformation_matrices, s; check=false)

  E_to_T = hom(E, T, E_to_T_mat; check=false)
  T_to_E = hom(T, E, T_to_E_mat; check=false)

  function my_mult(as::Tuple{Vararg{LieAlgebraModuleElem{C}}})
    @req length(as) == k "Length of vector does not match."
    @req all(a -> parent(a) === V, as) "Incompatible modules."
    mat = zero_matrix(R, 1, dim_E)
    for (i, _inds) in enumerate(ind_map), (inds, sgn) in permutations_with_sign(_inds)
      mat[1, i] +=
        sgn * prod(as[j].mat[k]::elem_type(R) for (j, k) in enumerate(inds); init=one(R))
    end
    return E(mat)
  end
  function my_mult(as::LieAlgebraModuleElem{C}...)
    return my_mult(as)::LieAlgebraModuleElem{C}
  end

  function my_decomp(a::LieAlgebraModuleElem{C})
    @req parent(a) === E "Incompatible modules."
    if iszero(a)
      return Tuple(zero(V) for _ in 1:k)
    end
    nz = findall(!iszero, coefficients(a))
    @req length(nz) == 1 "Non-pure exterior power element."
    @req isone(coeff(a, nz[1])) "Non-pure exterior power element."
    inds = ind_map[nz[1]]
    return Tuple(basis(V, i) for i in inds)
  end

  mult_map = MapFromFunc(
    Hecke.TupleParent(Tuple([zero(V) for _ in 1:k])), E, my_mult, my_decomp
  )
  inv_mult_map = MapFromFunc(E, domain(mult_map), my_decomp, my_mult)

  set_attribute!(
    E,
    :is_exterior_power => (V, k),
    :wedge_pure_function => mult_map,
    :wedge_generator_decompose_function => inv_mult_map,
    :embedding_tensor_power => T,
    :embedding_tensor_power_embedding => E_to_T,
    :embedding_tensor_power_projection => T_to_E,
  )

  cached && (_exterior_powers(V)[k] = (E, mult_map))

  return E, mult_map
end

######################## Symmetric powers ########################

# cache storage
@attr Dict{Int,Tuple{typeof(V),MapFromFunc}} function _symmetric_powers(V::LieAlgebraModule)
  return Dict{Int,Tuple{typeof(V),MapFromFunc}}()
end

@doc raw"""
    symmetric_power(V::LieAlgebraModule{C}, k::Int; cached::Bool=true) -> LieAlgebraModule{C}, Map

Construct the `k`-th symmetric power $S^k (V)$ of the module `V`,
together with a map that computes the product of `k` elements of `V`.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 4);

julia> V = exterior_power(standard_module(L), 3)[1]; # some module

julia> S, map = symmetric_power(V, 2)
(Symmetric power module of dimension 10 over L, Map: parent of tuples of type Tuple{LieAlgebraModuleElem{QQFieldElem, LinearLieAlgebraElem{QQFieldElem}}, LieAlgebraModuleElem{QQFieldElem, LinearLieAlgebraElem{QQFieldElem}}} -> S)

julia> S
Symmetric power module
  of dimension 10
  2nd symmetric power of
    3rd exterior power of
      standard module
over special linear Lie algebra of degree 4 over QQ

julia> basis(S)
10-element Vector{LieAlgebraModuleElem{QQFieldElem, LinearLieAlgebraElem{QQFieldElem}}}:
 (v_1^v_2^v_3)^2
 (v_1^v_2^v_3)*(v_1^v_2^v_4)
 (v_1^v_2^v_3)*(v_1^v_3^v_4)
 (v_1^v_2^v_3)*(v_2^v_3^v_4)
 (v_1^v_2^v_4)^2
 (v_1^v_2^v_4)*(v_1^v_3^v_4)
 (v_1^v_2^v_4)*(v_2^v_3^v_4)
 (v_1^v_3^v_4)^2
 (v_1^v_3^v_4)*(v_2^v_3^v_4)
 (v_2^v_3^v_4)^2
```
"""
function symmetric_power(
  V::LieAlgebraModule{C}, k::Int; cached::Bool=true
) where {C<:FieldElem}
  @req k >= 0 "Non-negative exponent needed"
  @req characteristic(coefficient_ring(V)) == 0 "Characteristic must be zero."

  if cached
    powers = _symmetric_powers(V)
    haskey(powers, k) && return powers[k]
  end

  L = base_lie_algebra(V)
  R = coefficient_ring(V)
  dim_S = binomial(dim(V) + k - 1, k)
  ind_map = collect(multicombinations(1:dim(V), k))

  T, _ = tensor_power(V, k; cached)
  S_to_T_mat = zero_matrix(R, dim_S, dim(T))
  T_to_S_mat = zero_matrix(R, dim(T), dim_S)
  for (i, _inds) in enumerate(ind_map), inds in permutations(_inds)
    j = 1 + sum((ind - 1) * dim(V)^(k - 1) for (k, ind) in enumerate(reverse(inds)); init=0)
    S_to_T_mat[i, j] += R(1)//factorial(k)
    T_to_S_mat[j, i] = R(1)
  end
  transformation_matrices = map(1:dim(L)) do i
    S_to_T_mat * transformation_matrix(T, i) * T_to_S_mat
  end

  s = if k == 0
    [Symbol("1")]
  elseif k == 1
    symbols(V)
  elseif _is_standard_module(V)
    [
      Symbol(
        join(
          (
            (e = count(==(i), inds)) == 1 ? s : "$(s)^$e" for
            (i, s) in enumerate(symbols(V)) if in(i, inds)
          ),
          "*",
        ),
      ) for inds in ind_map
    ]
  else
    [
      Symbol(
        join(
          (
            (e = count(==(i), inds)) == 1 ? "($(s))" : "($(s))^$e" for
            (i, s) in enumerate(symbols(V)) if in(i, inds)
          ),
          "*",
        ),
      ) for inds in ind_map
    ]
  end

  S = LieAlgebraModule{C}(L, dim_S, transformation_matrices, s; check=false)

  S_to_T = hom(S, T, S_to_T_mat; check=false)
  T_to_S = hom(T, S, T_to_S_mat; check=false)

  function my_mult(as::Tuple{Vararg{LieAlgebraModuleElem{C}}})
    @req length(as) == k "Length of vector does not match."
    @req all(a -> parent(a) === V, as) "Incompatible modules."
    mat = zero_matrix(R, 1, dim_S)
    for (i, _inds) in enumerate(ind_map), inds in unique(permutations(_inds))
      mat[1, i] += prod(
        as[j].mat[k]::elem_type(R) for (j, k) in enumerate(inds); init=one(R)
      )
    end
    return S(mat)
  end
  function my_mult(as::LieAlgebraModuleElem{C}...)
    return my_mult(as)::LieAlgebraModuleElem{C}
  end

  function my_decomp(a::LieAlgebraModuleElem{C})
    @req parent(a) === S "Incompatible modules."
    if iszero(a)
      return Tuple(zero(V) for _ in 1:k)
    end
    nz = findall(!iszero, coefficients(a))
    @req length(nz) == 1 "Non-pure exterior power element."
    @req isone(coeff(a, nz[1])) "Non-pure exterior power element."
    inds = ind_map[nz[1]]
    return Tuple(basis(V, i) for i in inds)
  end

  mult_map = MapFromFunc(
    Hecke.TupleParent(Tuple([zero(V) for _ in 1:k])), S, my_mult, my_decomp
  )
  inv_mult_map = MapFromFunc(S, domain(mult_map), my_decomp, my_mult)

  set_attribute!(
    S,
    :is_symmetric_power => (V, k),
    :mult_pure_function => mult_map,
    :mult_generator_decompose_function => inv_mult_map,
    :embedding_tensor_power => T,
    :embedding_tensor_power_embedding => S_to_T,
    :embedding_tensor_power_projection => T_to_S,
  )

  cached && (_symmetric_powers(V)[k] = (S, mult_map))

  return S, mult_map
end

######################## Tensor powers ########################

# cache storage
@attr Dict{Int,Tuple{typeof(V),MapFromFunc}} function _tensor_powers(V::LieAlgebraModule)
  return Dict{Int,Tuple{typeof(V),MapFromFunc}}()
end

@doc raw"""
    tensor_power(V::LieAlgebraModule{C}, k::Int; cached::Bool=true) -> LieAlgebraModule{C}, Map

Construct the `k`-th tensor power $T^k (V)$ of the module `V`,
together with a map that computes the tensor product of `k` elements of `V`.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 3);

julia> V = exterior_power(standard_module(L), 2)[1]; # some module

julia> T, map = tensor_power(V, 2)
(Tensor power module of dimension 9 over L, Map: parent of tuples of type Tuple{LieAlgebraModuleElem{QQFieldElem, LinearLieAlgebraElem{QQFieldElem}}, LieAlgebraModuleElem{QQFieldElem, LinearLieAlgebraElem{QQFieldElem}}} -> T)

julia> T
Tensor power module
  of dimension 9
  2nd tensor power of
    2nd exterior power of
      standard module
over special linear Lie algebra of degree 3 over QQ

julia> basis(T)
9-element Vector{LieAlgebraModuleElem{QQFieldElem, LinearLieAlgebraElem{QQFieldElem}}}:
 (v_1^v_2)(x)(v_1^v_2)
 (v_1^v_2)(x)(v_1^v_3)
 (v_1^v_2)(x)(v_2^v_3)
 (v_1^v_3)(x)(v_1^v_2)
 (v_1^v_3)(x)(v_1^v_3)
 (v_1^v_3)(x)(v_2^v_3)
 (v_2^v_3)(x)(v_1^v_2)
 (v_2^v_3)(x)(v_1^v_3)
 (v_2^v_3)(x)(v_2^v_3)
```
"""
function tensor_power(
  V::LieAlgebraModule{C}, k::Int; cached::Bool=true
) where {C<:FieldElem}
  @req k >= 0 "Non-negative exponent needed"

  if cached
    powers = _tensor_powers(V)
    haskey(powers, k) && return powers[k]
  end

  L = base_lie_algebra(V)
  R = coefficient_ring(V)
  dim_T = dim(V)^k
  ind_map = k > 0 ? reverse.(collect(ProductIterator(1:dim(V), k))) : [Int[]]

  transformation_matrices = map(1:dim(L)) do i
    y = transformation_matrix(V, i)
    sum(
      kronecker_product(
        kronecker_product(identity_matrix(R, dim(V)^(j - 1)), y),
        identity_matrix(R, dim(V)^(k - j)),
      ) for j in 1:k;
      init=zero_matrix(R, dim_T, dim_T),
    )
  end

  s = if k == 0
    [Symbol("1")]
  elseif k == 1
    symbols(V)
  elseif _is_standard_module(V)
    [
      Symbol(join(s, (is_unicode_allowed() ? "⊗" : "(x)"))) for
      s in reverse.(ProductIterator(symbols(V), k))
    ]
  else
    [
      Symbol(join((x -> "($x)").(s), (is_unicode_allowed() ? "⊗" : "(x)"))) for
      s in reverse.(ProductIterator(symbols(V), k))
    ]
  end

  T = LieAlgebraModule{C}(L, dim_T, transformation_matrices, s; check=false)

  function my_mult(as::Tuple{Vararg{LieAlgebraModuleElem{C}}})
    @req length(as) == k "Length of vector does not match."
    @req all(a -> parent(a) === V, as) "Incompatible modules."
    mat = zero_matrix(R, 1, dim_T)
    for (i, inds) in enumerate(ind_map)
      mat[1, i] += prod(
        as[j].mat[k]::elem_type(R) for (j, k) in enumerate(inds); init=one(R)
      )
    end
    return T(mat)
  end
  function my_mult(as::LieAlgebraModuleElem{C}...)
    return my_mult(as)::LieAlgebraModuleElem{C}
  end

  function my_decomp(a::LieAlgebraModuleElem{C})
    @req parent(a) === T "Incompatible modules."
    if iszero(a)
      return Tuple(zero(V) for _ in 1:k)
    end
    nz = findall(!iszero, coefficients(a))
    @req length(nz) == 1 "Non-pure tensor power element."
    @req isone(coeff(a, nz[1])) "Non-pure tensor power element."
    inds = ind_map[nz[1]]
    return Tuple(basis(V, i) for i in inds)
  end

  mult_map = MapFromFunc(
    Hecke.TupleParent(Tuple([zero(V) for _ in 1:k])), T, my_mult, my_decomp
  )
  inv_mult_map = MapFromFunc(T, domain(mult_map), my_decomp, my_mult)

  set_attribute!(
    T,
    :is_tensor_power => (V, k),
    :tensor_pure_function => mult_map,
    :tensor_generator_decompose_function => inv_mult_map,
  )

  cached && (_tensor_powers(V)[k] = (T, mult_map))

  return T, mult_map
end
