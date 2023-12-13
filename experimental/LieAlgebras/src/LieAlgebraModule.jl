@attributes mutable struct LieAlgebraModule{C<:FieldElem}
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
  ) where {C<:FieldElem}
    @req dimV == length(s) "Invalid number of basis element names."
    @req dim(L) == length(transformation_matrices) "Invalid number of transformation matrices."
    @req all(m -> size(m) == (dimV, dimV), transformation_matrices) "Invalid transformation matrix dimensions."

    V = new{C}(L, dimV, transformation_matrices, s)
    if check
      @req all(m -> all(e -> parent(e) === coefficient_ring(V), m), transformation_matrices) "Invalid transformation matrix entries."
      for xi in basis(L), xj in basis(L), v in basis(V)
        @req (xi * xj) * v == xi * (xj * v) - xj * (xi * v) "Transformation matrices do not define a module."
      end
    end
    return V
  end
end

struct LieAlgebraModuleElem{C<:FieldElem}
  parent::LieAlgebraModule{C}
  mat::MatElem{C}
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{LieAlgebraModuleElem{C}}) where {C<:FieldElem} = LieAlgebraModule{C}

elem_type(::Type{LieAlgebraModule{C}}) where {C<:FieldElem} = LieAlgebraModuleElem{C}

parent(v::LieAlgebraModuleElem) = v.parent

coefficient_ring(V::LieAlgebraModule) = coefficient_ring(base_lie_algebra(V))

coefficient_ring(v::LieAlgebraModuleElem) = coefficient_ring(parent(v))

@doc raw"""
    base_lie_algebra(V::LieAlgebraModule{C}) -> LieAlgebra{C}

Return the Lie algebra `V` is a module over.
"""
base_lie_algebra(V::LieAlgebraModule) = V.L

ngens(L::LieAlgebraModule) = dim(L)

gens(L::LieAlgebraModule) = basis(L)

gen(L::LieAlgebraModule, i::Int) = basis(L, i)

@doc raw"""
    dim(V::LieAlgebraModule{C}) -> Int

Return the dimension of the Lie algebra module `V`.
"""
dim(V::LieAlgebraModule) = V.dim

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
  return V([(j == i ? one(R) : zero(R)) for j in 1:dim(V)])
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

function Base.show(io::IO, ::MIME"text/plain", V::LieAlgebraModule)
  io = pretty(io)
  println(io, _module_type_to_string(get_attribute(V, :type, :unknown)))
  println(io, Indent(), "of dimension $(dim(V))")
  if is_dual(V) ||
    is_direct_sum(V) ||
    is_tensor_product(V) ||
    is_exterior_power(V) ||
    is_symmetric_power(V) ||
    is_tensor_power(V)
    _show_inner(io, V)
  end
  print(io, Dedent())
  print(io, "over ")
  print(io, Lowercase(), base_lie_algebra(V))
end

function _show_inner(io::IO, V::LieAlgebraModule)
  type = get_attribute(V, :type, :unknown)
  if type == :standard_module
    println(io, "standard module")
  elseif type == :dual
    println(io, "dual of ", Lowercase())
    print(io, Indent())
    _show_inner(io, base_module(V))
    print(io, Dedent())
  elseif type == :direct_sum
    println(io, "direct sum with direct summands")
    print(io, Indent())
    for W in base_modules(V)
      _show_inner(io, W)
    end
    print(io, Dedent())
  elseif type == :tensor_product
    println(io, "tensor product with tensor factors")
    print(io, Indent())
    for W in base_modules(V)
      _show_inner(io, W)
    end
    print(io, Dedent())
  elseif type == :exterior_power
    println(io, "$(ordinal_number_string(get_attribute(V, :power))) exterior power of")
    print(io, Indent())
    _show_inner(io, base_module(V))
    print(io, Dedent())
  elseif type == :symmetric_power
    println(io, "$(ordinal_number_string(get_attribute(V, :power))) symmetric power of")
    print(io, Indent())
    _show_inner(io, base_module(V))
    print(io, Dedent())
  elseif type == :tensor_power
    println(io, "$(ordinal_number_string(get_attribute(V, :power))) tensor power of")
    print(io, Indent())
    _show_inner(io, base_module(V))
    print(io, Dedent())
  else
    println(io, "abstract module")
  end
end

function Base.show(io::IO, V::LieAlgebraModule)
  if get(io, :supercompact, false)
    print(io, _module_type_to_string(get_attribute(V, :type, :unknown)))
  else
    io = pretty(io)
    print(
      io,
      _module_type_to_string(get_attribute(V, :type, :unknown)),
      " of dimension $(dim(V)) over ",
      Lowercase(),
    )
    print(IOContext(io, :supercompact => true), base_lie_algebra(V))
  end
end

function _module_type_to_string(type::Symbol)
  if type == :standard_module
    return "Standard module"
  elseif type == :dual
    return "Dual module"
  elseif type == :direct_sum
    return "Direct sum module"
  elseif type == :tensor_product
    return "Tensor product module"
  elseif type == :exterior_power
    return "Exterior power module"
  elseif type == :symmetric_power
    return "Symmetric power module"
  elseif type == :tensor_power
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
    (V::LieAlgebraModule{C})(v::Vector{Int}) -> LieAlgebraModuleElem{C}

Return the element of `V` with coefficient vector `v`.
Fails, if `Int` cannot be coerced into the base ring of `L`.
"""
function (V::LieAlgebraModule)(v::Vector{Int})
  return V(coefficient_ring(V).(v))
end

@doc raw"""
    (V::LieAlgebraModule{C})(v::Vector{C}) -> LieAlgebraModuleElem{C}

Return the element of `V` with coefficient vector `v`.
"""
function (V::LieAlgebraModule{C})(v::Vector{C}) where {C<:FieldElem}
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
  if is_dual(V) && base_module(V) === parent(v)
    return V(coefficients(v))
  end
  @req V === parent(v) "Incompatible modules."
  return v
end

@doc raw"""
    (V::LieAlgebraModule{C})(a::Vector{T}) where {T<:LieAlgebraModuleElem{C}}) -> LieAlgebraModuleElem{C}

If `V` is a direct sum, return its element, where the $i$-th component is equal to `a[i]`.
If `V` is a tensor product, return the tensor product of the `a[i]`.
If `V` is a exterior (symmetric, tensor) power, return the wedge product
(product, tensor product) of the `a[i]`.

Requires that `a` has a suitable length, and that the `a[i]` are elements of the correct modules,
where _correct_ depends on the case above.
"""
function (V::LieAlgebraModule{C})(
  a::Vector{T}
) where {T<:LieAlgebraModuleElem{C}} where {C<:FieldElem}
  if is_direct_sum(V)
    @req length(a) == length(base_modules(V)) "Length of vector does not match."
    @req all(i -> parent(a[i]) === base_modules(V)[i], 1:length(a)) "Incompatible modules."
    return sum(inji(ai) for (ai, inji) in zip(a, canonical_injections(V)); init=zero(V))
  elseif is_tensor_product(V)
    pure = get_attribute(V, :tensor_pure_function)
    return pure(a)::LieAlgebraModuleElem{C}
  elseif is_exterior_power(V)
    pure = get_attribute(V, :exterior_pure_function)
    return pure(a)::LieAlgebraModuleElem{C}
  elseif is_symmetric_power(V)
    pure = get_attribute(V, :symmetric_pure_function)
    return pure(a)::LieAlgebraModuleElem{C}
  elseif is_tensor_power(V)
    pure = get_attribute(V, :tensor_pure_function)
    return pure(a)::LieAlgebraModuleElem{C}
  else
    throw(ArgumentError("Invalid input."))
  end
end

function (V::LieAlgebraModule{C})(
  a::Tuple{T,Vararg{T}}
) where {T<:LieAlgebraModuleElem{C}} where {C<:FieldElem}
  return V(collect(a))
end

function (V::LieAlgebraModule{C})(_::Tuple{}) where {C<:FieldElem}
  return V(LieAlgebraModuleElem{C}[])
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
  return V1.dim == V2.dim &&
         V1.s == V2.s &&
         V1.L == V2.L &&
         V1.transformation_matrices == V2.transformation_matrices
end

function Base.hash(V::LieAlgebraModule, h::UInt)
  b = 0x28b0c111e3ff8526 % UInt
  h = hash(V.dim, h)
  h = hash(V.s, h)
  h = hash(V.L, h)
  h = hash(V.transformation_matrices, h)
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

function transformation_matrix(V::LieAlgebraModule{C}, i::Int) where {C<:FieldElem}
  return (V.transformation_matrices[i])::dense_matrix_type(C)
end

###############################################################################
#
#   Attribute getters
#
###############################################################################

@doc raw"""
    is_standard_module(V::LieAlgebraModule{C}) -> Bool

Check whether `V` has been constructed as a standard module.
"""
function is_standard_module(V::LieAlgebraModule)
  return get_attribute(V, :type, :fallback)::Symbol == :standard_module
end

@doc raw"""
    is_dual(V::LieAlgebraModule{C}) -> Bool

Check whether `V` has been constructed as a dual module.
"""
function is_dual(V::LieAlgebraModule)
  return get_attribute(V, :type, :fallback)::Symbol == :dual
end

@doc raw"""
    is_direct_sum(V::LieAlgebraModule{C}) -> Bool

Check whether `V` has been constructed as a direct sum of modules.
"""
function is_direct_sum(V::LieAlgebraModule)
  return get_attribute(V, :type, :fallback)::Symbol == :direct_sum
end

@doc raw"""
    is_tensor_product(V::LieAlgebraModule{C}) -> Bool

Check whether `V` has been constructed as a tensor product of modules.
"""
function is_tensor_product(V::LieAlgebraModule)
  return get_attribute(V, :type, :fallback)::Symbol == :tensor_product
end

@doc raw"""
    is_exterior_power(V::LieAlgebraModule{C}) -> Bool

Check whether `V` has been constructed as an exterior power of a module.
"""
function is_exterior_power(V::LieAlgebraModule)
  return get_attribute(V, :type, :fallback)::Symbol == :exterior_power
end

@doc raw"""
    is_symmetric_power(V::LieAlgebraModule{C}) -> Bool

Check whether `V` has been constructed as a symmetric power of a module.
"""
function is_symmetric_power(V::LieAlgebraModule)
  return get_attribute(V, :type, :fallback)::Symbol == :symmetric_power
end

@doc raw"""
    is_tensor_power(V::LieAlgebraModule{C}) -> Bool

Check whether `V` has been constructed as a tensor power of a module.
"""
function is_tensor_power(V::LieAlgebraModule)
  return get_attribute(V, :type, :fallback)::Symbol == :tensor_power
end

@doc raw"""
    base_module(V::LieAlgebraModule{C}) -> LieAlgebraModule{C}

Returns the base module of `V`, if `V` has been constructed as a power module.
"""
function base_module(V::LieAlgebraModule{C}) where {C<:FieldElem}
  @req is_dual(V) || is_exterior_power(V) || is_symmetric_power(V) || is_tensor_power(V) "Not a power module."
  return get_attribute(V, :base_module)::LieAlgebraModule{C}
end

@doc raw"""
    base_modules(V::LieAlgebraModule{C}) -> Vector{LieAlgebraModule{C}}

Returns the summands or tensor factors of `V`,
if `V` has been constructed as a direct sum or tensor product of modules.
"""
function base_modules(V::LieAlgebraModule{C}) where {C<:FieldElem}
  if is_tensor_product(V)
    return get_attribute(V, :tensor_product)::Vector{LieAlgebraModule{C}}
  elseif is_direct_sum(V)
    return get_attribute(V, :direct_sum)::Vector{LieAlgebraModule{C}}
  elseif is_tensor_power(V)
    return [base_module(V) for _ in 1:get_attribute(V, :power)]
  else
    error("Not a direct sum or tensor product module.")
  end
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
    abstract_module(L::LieAlgebra{C}, dimV::Int, struct_consts::Matrix{SRow{C}}, s::Vector{<:VarName}; check::Bool) -> LieAlgebraModule{C}

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
  struct_consts::Matrix{SRow{C}},
  s::Vector{<:VarName}=[Symbol("v_$i") for i in 1:dimV];
  check::Bool=true,
) where {C<:FieldElem}
  @req dim(L) == size(struct_consts, 1) "Invalid structure constants dimensions."
  @req dimV == size(struct_consts, 2) "Invalid structure constants dimensions."
  @req dimV == length(s) "Invalid number of basis element names."

  transformation_matrices = [zero_matrix(coefficient_ring(L), dimV, dimV) for _ in 1:dim(L)]
  for i in 1:dim(L), j in 1:dimV
    transformation_matrices[i][j, :] = dense_row(struct_consts[i, j], dimV)
  end

  return LieAlgebraModule{C}(L, dimV, transformation_matrices, Symbol.(s); check)
end

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
  # set_attribute!(triv_V, :type => :trivial)
  return triv_V
end

@doc raw"""
    standard_module(L::LinearLieAlgebra{C}) -> LieAlgebraModule{C}

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
function standard_module(L::LinearLieAlgebra)
  dim_std_V = L.n
  transformation_matrices = transpose.(matrix_repr_basis(L))
  s = [Symbol("v_$(i)") for i in 1:dim_std_V]
  std_V = LieAlgebraModule{elem_type(coefficient_ring(L))}(
    L, dim_std_V, transformation_matrices, s; check=false
  )
  set_attribute!(std_V, :type => :standard_module)
  return std_V
end

@doc raw"""
    dual(V::LieAlgebraModule{C}) -> LieAlgebraModule{C}

Construct the dual module of `V`.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 3);

julia> V = exterior_power(standard_module(L), 2); # some module

julia> dual(V)
Dual module
  of dimension 3
  dual of
    2nd exterior power of
      standard module
over special linear Lie algebra of degree 3 over QQ
```
"""
function dual(V::LieAlgebraModule{C}) where {C<:FieldElem}
  L = base_lie_algebra(V)
  dim_dual_V = dim(V)

  transformation_matrices = map(1:dim(L)) do i
    -transpose(transformation_matrix(V, i))
  end

  s = if is_standard_module(V)
    [Symbol("$(s)*") for s in symbols(V)]
  else
    [Symbol("($(s))*") for s in symbols(V)]
  end

  pow_V = LieAlgebraModule{C}(L, dim_dual_V, transformation_matrices, s; check=false)
  set_attribute!(pow_V, :type => :dual, :base_module => V)
  return pow_V
end

Base.:^(V::LieAlgebraModule, ::typeof(Base.:*)) = dual(V)

@doc raw"""
    direct_sum(V::LieAlgebraModule{C}...) -> LieAlgebraModule{C}
    ⊕(V::LieAlgebraModule{C}...) -> LieAlgebraModule{C}

Construct the direct sum of the modules `V...`.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 3);

julia> V1 = exterior_power(standard_module(L), 2); # some module

julia> V2 = symmetric_power(standard_module(L), 3); # some module

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
      s in (is_standard_module(Vj) ? symbols(Vj) : (x -> "($x)").(symbols(Vj)))
    ]
  end

  direct_sum_V = LieAlgebraModule{C}(
    L, dim_direct_sum_V, transformation_matrices, s; check=false
  )
  set_attribute!(direct_sum_V, :type => :direct_sum, :direct_sum => Vs)
  return direct_sum_V
end

⊕(V::LieAlgebraModule{C}, Vs::LieAlgebraModule{C}...) where {C<:FieldElem} =
  direct_sum(V, Vs...)

@doc raw"""
  tensor_product(Vs::LieAlgebraModule{C}...) -> LieAlgebraModule{C}
  ⊗(Vs::LieAlgebraModule{C}...) -> LieAlgebraModule{C}

Given modules $V_1,\dots,V_n$ over the same Lie algebra $L$,
construct their tensor product $V_1 \otimes \cdots \otimes \V_n$.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 3);

julia> V1 = exterior_power(standard_module(L), 2); # some module

julia> V2 = symmetric_power(standard_module(L), 3); # some module

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
      Symbol(join(s, " ⊗ ")) for s in
      reverse.(
        ProductIterator([
          is_standard_module(Vi) ? symbols(Vi) : (x -> "($x)").(symbols(Vi)) for
          Vi in reverse(Vs)
        ])
      )
    ]
  end

  tensor_product_V = LieAlgebraModule{C}(
    L, dim_tensor_product_V, transformation_matrices, s; check=false
  )

  function pure(as::LieAlgebraModuleElem{C}...)
    @req length(as) == length(Vs) "Length of vector does not match."
    @req all(i -> parent(as[i]) === Vs[i], 1:length(as)) "Incompatible modules."
    mat = zero_matrix(R, 1, dim_tensor_product_V)
    for (i, inds) in enumerate(ind_map)
      mat[1, i] += prod(as[j].mat[k]::elem_type(R) for (j, k) in enumerate(inds))
    end
    return tensor_product_V(mat)
  end
  function pure(as::Tuple)
    return pure(as...)::LieAlgebraModuleElem{C}
  end
  function pure(as::Vector{LieAlgebraModuleElem{C}})
    return pure(as...)::LieAlgebraModuleElem{C}
  end

  function inv_pure(a::LieAlgebraModuleElem{C})
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

  set_attribute!(
    tensor_product_V,
    :type => :tensor_product,
    :tensor_product => Vs,
    :tensor_pure_function => pure,
    :tensor_pure_preimage_function => inv_pure,
  )

  return tensor_product_V
end

⊗(V::LieAlgebraModule{C}, Vs::LieAlgebraModule{C}...) where {C<:FieldElem} =
  tensor_product(V, Vs...)

@doc raw"""
    exterior_power(V::LieAlgebraModule{C}, k::Int) -> LieAlgebraModule{C}

Construct the `k`-th exterior power $\bigwedge^k (V)$ of the module `V`.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 3);

julia> V = symmetric_power(standard_module(L), 3); # some module

julia> exterior_power(V, 2)
Exterior power module
  of dimension 45
  2nd exterior power of
    3rd symmetric power of
      standard module
over special linear Lie algebra of degree 3 over QQ
```
"""
function exterior_power(V::LieAlgebraModule{C}, k::Int) where {C<:FieldElem}
  @req k >= 0 "Non-negative exponent needed"
  L = base_lie_algebra(V)
  R = coefficient_ring(V)
  dim_E = binomial(dim(V), k)
  ind_map = collect(combinations(1:dim(V), k))

  T = tensor_power(V, k)
  E_to_T_mat = zero_matrix(coefficient_ring(V), dim_E, dim(T))
  T_to_E_mat = zero_matrix(coefficient_ring(V), dim(T), dim_E)
  for (i, _inds) in enumerate(ind_map), (inds, sgn) in permutations_with_sign(_inds)
    j = 1 + sum((ind - 1) * dim(V)^(k - 1) for (k, ind) in enumerate(reverse(inds)); init=0)
    E_to_T_mat[i, j] = sgn//factorial(k)
    T_to_E_mat[j, i] = sgn
  end
  transformation_matrices = map(1:dim(L)) do i
    E_to_T_mat * transformation_matrix(T, i) * T_to_E_mat
  end

  s = if k == 0
    [Symbol("1")]
  elseif k == 1
    symbols(V)
  elseif is_standard_module(V)
    [Symbol(join(s, " ∧ ")) for s in combinations(symbols(V), k)]
  else
    [Symbol(join((x -> "($x)").(s), " ∧ ")) for s in combinations(symbols(V), k)]
  end

  E = LieAlgebraModule{C}(L, dim_E, transformation_matrices, s; check=false)

  E_to_T = hom(E, T, E_to_T_mat; check=false)
  T_to_E = hom(T, E, T_to_E_mat; check=false)

  function pure(as::LieAlgebraModuleElem{C}...)
    @req length(as) == k "Length of vector does not match."
    @req all(a -> parent(a) === V, as) "Incompatible modules."
    mat = zero_matrix(R, 1, dim_E)
    for (i, _inds) in enumerate(ind_map), (inds, sgn) in permutations_with_sign(_inds)
      mat[1, i] +=
        sgn * prod(as[j].mat[k]::elem_type(R) for (j, k) in enumerate(inds); init=one(R))
    end
    return E(mat)
  end
  function pure(as::Tuple)
    return pure(as...)::LieAlgebraModuleElem{C}
  end
  function pure(as::Vector{LieAlgebraModuleElem{C}})
    return pure(as...)::LieAlgebraModuleElem{C}
  end

  function inv_pure(a::LieAlgebraModuleElem{C})
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

  set_attribute!(
    E,
    :type => :exterior_power,
    :power => k,
    :base_module => V,
    :exterior_pure_function => pure,
    :exterior_pure_preimage_function => inv_pure,
    :embedding_tensor_power => T,
    :embedding_tensor_power_embedding => E_to_T,
    :embedding_tensor_power_projection => T_to_E,
  )
  return E
end

@doc raw"""
    symmetric_power(V::LieAlgebraModule{C}, k::Int) -> LieAlgebraModule{C}

Construct the `k`-th symmetric power $S^k (V)$ of the module `V`.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 3);

julia> V = exterior_power(standard_module(L), 2); # some module

julia> symmetric_power(V, 3)
Symmetric power module
  of dimension 10
  3rd symmetric power of
    2nd exterior power of
      standard module
over special linear Lie algebra of degree 3 over QQ
```
"""
function symmetric_power(V::LieAlgebraModule{C}, k::Int) where {C<:FieldElem}
  @req k >= 0 "Non-negative exponent needed"
  L = base_lie_algebra(V)
  R = coefficient_ring(V)
  dim_S = binomial(dim(V) + k - 1, k)
  ind_map = collect(multicombinations(1:dim(V), k))

  T = tensor_power(V, k)
  S_to_T_mat = zero_matrix(coefficient_ring(V), dim_S, dim(T))
  T_to_S_mat = zero_matrix(coefficient_ring(V), dim(T), dim_S)
  for (i, _inds) in enumerate(ind_map), inds in permutations(_inds)
    j = 1 + sum((ind - 1) * dim(V)^(k - 1) for (k, ind) in enumerate(reverse(inds)); init=0)
    S_to_T_mat[i, j] += 1//factorial(k)
    T_to_S_mat[j, i] = 1
  end
  transformation_matrices = map(1:dim(L)) do i
    S_to_T_mat * transformation_matrix(T, i) * T_to_S_mat
  end

  s = if k == 0
    [Symbol("1")]
  elseif k == 1
    symbols(V)
  elseif is_standard_module(V)
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

  function pure(as::LieAlgebraModuleElem{C}...)
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
  function pure(as::Tuple)
    return pure(as...)::LieAlgebraModuleElem{C}
  end
  function pure(as::Vector{LieAlgebraModuleElem{C}})
    return pure(as...)::LieAlgebraModuleElem{C}
  end

  function inv_pure(a::LieAlgebraModuleElem{C})
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

  set_attribute!(
    S,
    :type => :symmetric_power,
    :power => k,
    :base_module => V,
    :symmetric_pure_function => pure,
    :symmetric_pure_preimage_function => inv_pure,
    :embedding_tensor_power => T,
    :embedding_tensor_power_embedding => S_to_T,
    :embedding_tensor_power_projection => T_to_S,
  )
  return S
end

@doc raw"""
    tensor_power(V::LieAlgebraModule{C}, k::Int) -> LieAlgebraModule{C}

Construct the `k`-th tensor power $T^k (V)$ of the module `V`.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 3);

julia> V = exterior_power(standard_module(L), 2); # some module

julia> tensor_power(V, 3)
Tensor power module
  of dimension 27
  3rd tensor power of
    2nd exterior power of
      standard module
over special linear Lie algebra of degree 3 over QQ
```
"""
function tensor_power(V::LieAlgebraModule{C}, k::Int) where {C<:FieldElem}
  @req k >= 0 "Non-negative exponent needed"
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
  elseif is_standard_module(V)
    [Symbol(join(s, " ⊗ ")) for s in reverse.(ProductIterator(symbols(V), k))]
  else
    [Symbol(join((x -> "($x)").(s), " ⊗ ")) for s in reverse.(ProductIterator(symbols(V), k))]
  end

  T = LieAlgebraModule{C}(L, dim_T, transformation_matrices, s; check=false)

  function pure(as::LieAlgebraModuleElem{C}...)
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
  function pure(as::Tuple)
    return pure(as...)::LieAlgebraModuleElem{C}
  end
  function pure(as::Vector{LieAlgebraModuleElem{C}})
    return pure(as...)::LieAlgebraModuleElem{C}
  end

  function inv_pure(a::LieAlgebraModuleElem{C})
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

  set_attribute!(
    T,
    :type => :tensor_power,
    :power => k,
    :base_module => V,
    :tensor_pure_function => pure,
    :tensor_pure_preimage_function => inv_pure,
  )
  return T
end

###############################################################################
#
#   Simple modules (via highest weight) of semisimple Lie algebras
#
###############################################################################

# TODO: check semisimplicity check once that is available

function is_dominant_weight(hw::Vector{Int})
  return all(>=(0), hw)
end

@doc raw"""
    simple_module(L::LieAlgebra{C}, hw::Vector{Int}) -> LieAlgebraModule{C}

Construct the simple module of the Lie algebra `L` with highest weight `hw`.
"""
function simple_module(L::LieAlgebra, hw::Vector{Int})
  @req is_dominant_weight(hw) "Not a dominant weight."
  struct_consts = lie_algebra_simple_module_struct_consts_gap(L, hw)
  dimV = size(struct_consts, 2)
  V = abstract_module(L, dimV, struct_consts; check=false)
  # TODO: set appropriate attributes
  return V
end

@doc raw"""
    dim_of_simple_module([T = Int], L::LieAlgebra{C}, hw::Vector{Int}) -> T

Computes the dimension of the simple module of the Lie algebra `L` with highest weight `hw`.
The return value is of type `T`.

# Example
```jldoctest
julia> L = lie_algebra(QQ, :A, 3);

julia> dim_of_simple_module(L, [1, 1, 1])
64
```
"""
function dim_of_simple_module(T::Type, L::LieAlgebra, hw::Vector{Int})
  @req is_dominant_weight(hw) "Not a dominant weight."
  return T(
    GAPWrap.DimensionOfHighestWeightModule(codomain(Oscar.iso_oscar_gap(L)), GAP.Obj(hw))
  )
end

function dim_of_simple_module(L::LieAlgebra, hw::Vector{Int})
  return dim_of_simple_module(Int, L, hw)
end

@doc raw"""
    dominant_character(L::LieAlgebra{C}, hw::Vector{Int}) -> Dict{Vector{Int}, Int}

Computes the dominant weights occurring in the simple module of the Lie algebra `L` with highest weight `hw`,
together with their multiplicities.

# Example
```jldoctest
julia> L = lie_algebra(QQ, :A, 3);

julia> dominant_character(L, [2, 1, 0])
Dict{Vector{Int64}, Int64} with 4 entries:
  [2, 1, 0] => 1
  [1, 0, 1] => 2
  [0, 0, 0] => 3
  [0, 2, 0] => 1
```
"""
function dominant_character(L::LieAlgebra, hw::Vector{Int})
  @req is_dominant_weight(hw) "Not a dominant weight."
  return Dict{Vector{Int},Int}(
    Vector{Int}(w) => d for (w, d) in
    zip(GAPWrap.DominantCharacter(codomain(Oscar.iso_oscar_gap(L)), GAP.Obj(hw))...)
  )
end

@doc raw"""
    character(L::LieAlgebra{C}, hw::Vector{Int}) -> Dict{Vector{Int}, Int}

Computes all weights occurring in the simple module of the Lie algebra `L` with highest weight `hw`,
together with their multiplicities.

# Example
```jldoctest
julia> L = lie_algebra(QQ, :A, 3);

julia> character(L, [2, 0, 0])
Dict{Vector{Int64}, Int64} with 10 entries:
  [0, 1, 0]   => 1
  [0, -2, 2]  => 1
  [0, 0, -2]  => 1
  [-1, 1, -1] => 1
  [-2, 2, 0]  => 1
  [1, -1, 1]  => 1
  [-1, 0, 1]  => 1
  [1, 0, -1]  => 1
  [0, -1, 0]  => 1
  [2, 0, 0]   => 1
```
"""
function character(L::LieAlgebra, hw::Vector{Int})
  @req is_dominant_weight(hw) "Not a dominant weight."
  dc = dominant_character(L, hw)
  c = Dict{Vector{Int},Int}()
  W = GAPWrap.WeylGroup(GAPWrap.RootSystem(codomain(Oscar.iso_oscar_gap(L))))
  for (w, d) in dc
    it = GAPWrap.WeylOrbitIterator(W, GAP.Obj(w))
    while !GAPWrap.IsDoneIterator(it)
      push!(c, Vector{Int}(GAPWrap.NextIterator(it)) => d)
    end
  end
  return c
end

@doc raw"""
    tensor_product_decomposition(L::LieAlgebra, hw1::Vector{Int}, hw2::Vector{Int}) -> MSet{Vector{Int}}

Computes the decomposition of the tensor product of the simple modules of the Lie algebra `L` with highest weights `hw1` and `hw2`
into simple modules with their multiplicities.

# Example
```jldoctest
julia> L = lie_algebra(QQ, :A, 2);

julia> tensor_product_decomposition(L, [1, 0], [0, 1])
MSet{Vector{Int64}} with 2 elements:
  [0, 0]
  [1, 1]

julia> tensor_product_decomposition(L, [1, 1], [1, 1])
MSet{Vector{Int64}} with 6 elements:
  [0, 0]
  [1, 1] : 2
  [2, 2]
  [3, 0]
  [0, 3]
```
"""
function tensor_product_decomposition(L::LieAlgebra, hw1::Vector{Int}, hw2::Vector{Int})
  @req is_dominant_weight(hw1) && is_dominant_weight(hw2) "Both weights must be dominant."
  return multiset(
    Tuple{Vector{Vector{Int}},Vector{Int}}(
      GAPWrap.DecomposeTensorProduct(
        codomain(Oscar.iso_oscar_gap(L)), GAP.Obj(hw1), GAP.Obj(hw2)
      ),
    )...,
  )
end
