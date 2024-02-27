#################################################
#
# Abstract parent type
#
#################################################

abstract type LieAlgebra{C<:FieldElem} end

abstract type LieAlgebraElem{C<:FieldElem} end

# To be implemented by subtypes:
# parent_type(::Type{MyLieAlgebraElem{C}})
# elem_type(::Type{MyLieAlgebra{C}})
# parent(x::MyLieAlgebraElem{C})
# coefficient_ring(L::MyLieAlgebra{C})
# dim(L::MyLieAlgebra{C})
# Base.show(io::IO, x::MyLieAlgebra{C})
# symbols(L::MyLieAlgebra{C})
# bracket(x::MyLieAlgebraElem{C}, y::MyLieAlgebraElem{C})

###############################################################################
#
#   Basic manipulation
#
###############################################################################

coefficient_ring(x::LieAlgebraElem) = coefficient_ring(parent(x))

number_of_generators(L::LieAlgebra) = dim(L)

gens(L::LieAlgebra) = basis(L)

gen(L::LieAlgebra, i::Int) = basis(L, i)

@doc raw"""
    dim(L::LieAlgebra) -> Int

Return the dimension of the Lie algebra `L`.
"""
dim(_::LieAlgebra) = error("Should be implemented by subtypes.")

@doc raw"""
    basis(L::LieAlgebra{C}) -> Vector{LieAlgebraElem{C}}

Return a basis of the Lie algebra `L`.
"""
basis(L::LieAlgebra) = [basis(L, i)::elem_type(L) for i in 1:dim(L)]

@doc raw"""
    basis(L::LieAlgebra{C}, i::Int) -> LieAlgebraElem{C}

Return the `i`-th basis element of the Lie algebra `L`.
"""
function basis(L::LieAlgebra, i::Int)
  @req 1 <= i <= dim(L) "Index out of bounds."
  R = coefficient_ring(L)
  return L([(j == i ? one(R) : zero(R)) for j in 1:dim(L)])
end

@doc raw"""
    characteristic(L::LieAlgebra) -> Int

Return the characteristic of the coefficient ring of the Lie algebra `L`.
"""
characteristic(L::LieAlgebra) = characteristic(coefficient_ring(L))

@doc raw"""
    zero(L::LieAlgebra{C}) -> LieAlgebraElem{C}

Return the zero element of the Lie algebra `L`.
"""
function zero(L::LieAlgebra)
  mat = zero_matrix(coefficient_ring(L), 1, dim(L))
  return elem_type(L)(L, mat)
end

@doc raw"""
    iszero(x::LieAlgebraElem{C}) -> Bool

Check whether the Lie algebra element `x` is zero.
"""
function iszero(x::LieAlgebraElem)
  return iszero(coefficients(x))
end

@inline function _matrix(x::LieAlgebraElem{C}) where {C<:FieldElem}
  return (x.mat)::dense_matrix_type(C)
end

@doc raw"""
    coefficients(x::LieAlgebraElem{C}) -> Vector{C}

Return the coefficients of `x` with respect to [`basis(::LieAlgebra)`](@ref).
"""
function coefficients(x::LieAlgebraElem)
  return collect(_matrix(x))[1, :]
end

@doc raw"""
    coeff(x::LieAlgebraElem{C}, i::Int) -> C

Return the `i`-th coefficient of `x` with respect to [`basis(::LieAlgebra)`](@ref).
"""
function coeff(x::LieAlgebraElem, i::Int)
  return _matrix(x)[1, i]
end

@doc raw"""
    getindex(x::LieAlgebraElem{C}, i::Int) -> C

Return the `i`-th coefficient of `x` with respect to [`basis(::LieAlgebra)`](@ref).
"""
function getindex(x::LieAlgebraElem, i::Int)
  return coeff(x, i)
end

function Base.deepcopy_internal(x::LieAlgebraElem, dict::IdDict)
  return parent(x)(deepcopy_internal(_matrix(x), dict))
end

function check_parent(x1::LieAlgebraElem{C}, x2::LieAlgebraElem{C}) where {C<:FieldElem}
  parent(x1) !== parent(x2) && error("Incompatible Lie algebras.")
end

###############################################################################
#
#   String I/O
#
###############################################################################

@doc raw"""
    symbols(L::LieAlgebra{C}) -> Vector{Symbol}

Return the symbols used for printing basis elements of the Lie algebra `L`.
"""
symbols(_::LieAlgebra) = error("Should be implemented by subtypes.")

function expressify(v::LieAlgebraElem, s=symbols(parent(v)); context=nothing)
  sum = Expr(:call, :+)
  for (i, c) in enumerate(coefficients(v))
    push!(sum.args, Expr(:call, :*, expressify(c; context=context), s[i]))
  end
  return sum
end

@enable_all_show_via_expressify LieAlgebraElem

###############################################################################
#
#   Parent object call overload
#
###############################################################################

@doc raw"""
    (L::LieAlgebra{C})() -> LieAlgebraElem{C}

Return the zero element of the Lie algebra `L`.
"""
function (L::LieAlgebra)()
  return zero(L)
end

@doc raw"""
    (L::LieAlgebra{C})(v::Vector{Int}) -> LieAlgebraElem{C}

Return the element of `L` with coefficient vector `v`.
Fail, if `Int` cannot be coerced into the base ring of `L`.
"""
function (L::LieAlgebra)(v::Vector{Int})
  return L(coefficient_ring(L).(v))
end

@doc raw"""
    (L::LieAlgebra{C})(v::Vector{C}) -> LieAlgebraElem{C}

Return the element of `L` with coefficient vector `v`.
"""
function (L::LieAlgebra{C})(v::Vector{C}) where {C<:FieldElem}
  @req length(v) == dim(L) "Length of vector does not match dimension."
  mat = matrix(coefficient_ring(L), 1, length(v), v)
  return elem_type(L)(L, mat)
end

@doc raw"""
    (L::LieAlgebra{C})(mat::MatElem{C}) -> LieAlgebraElem{C}

Return the element of `L` with coefficient vector equivalent to
the $1 \times \dim(L)$ matrix `mat`.
"""
function (L::LieAlgebra{C})(mat::MatElem{C}) where {C<:FieldElem}
  @req size(mat) == (1, dim(L)) "Invalid matrix dimensions."
  return elem_type(L)(L, mat)
end

@doc raw"""
    (L::LieAlgebra{C})(v::SRow{C}) -> LieAlgebraElem{C}

Return the element of `L` with coefficient vector `v`.
"""
function (L::LieAlgebra{C})(v::SRow{C}) where {C<:FieldElem}
  mat = dense_row(v, dim(L))
  return elem_type(L)(L, mat)
end

@doc raw"""
    (L::LieAlgebra{C})(x::LieAlgebraElem{C}) -> LieAlgebraElem{C}

Return `x`. Fails if `x` is not an element of `L`.
"""
function (L::LieAlgebra{C})(x::LieAlgebraElem{C}) where {C<:FieldElem}
  @req L === parent(x) "Incompatible modules."
  return x
end

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

function Base.:-(x::LieAlgebraElem{C}) where {C<:FieldElem}
  return parent(x)(-_matrix(x))
end

function Base.:+(x1::LieAlgebraElem{C}, x2::LieAlgebraElem{C}) where {C<:FieldElem}
  check_parent(x1, x2)
  return parent(x1)(_matrix(x1) + _matrix(x2))
end

function Base.:-(x1::LieAlgebraElem{C}, x2::LieAlgebraElem{C}) where {C<:FieldElem}
  check_parent(x1, x2)
  return parent(x1)(_matrix(x1) - _matrix(x2))
end

function Base.:*(x::LieAlgebraElem{C}, c::C) where {C<:FieldElem}
  coefficient_ring(x) != parent(c) && error("Incompatible rings.")
  return parent(x)(_matrix(x) * c)
end

function Base.:*(x::LieAlgebraElem, c::U) where {U<:RationalUnion}
  return parent(x)(_matrix(x) * c)
end

function Base.:*(c::C, x::LieAlgebraElem{C}) where {C<:FieldElem}
  coefficient_ring(x) != parent(c) && error("Incompatible rings.")
  return parent(x)(c * _matrix(x))
end

function Base.:*(c::U, x::LieAlgebraElem) where {U<:RationalUnion}
  return parent(x)(c * _matrix(x))
end

function Base.:*(x::LieAlgebraElem{C}, y::LieAlgebraElem{C}) where {C<:FieldElem}
  return bracket(x, y)
end

###############################################################################
#
#   Comparison functions
#
###############################################################################

function Base.:(==)(x1::LieAlgebraElem{C}, x2::LieAlgebraElem{C}) where {C<:FieldElem}
  check_parent(x1, x2)
  return coefficients(x1) == coefficients(x2)
end

function Base.hash(x::LieAlgebraElem, h::UInt)
  b = 0x6724cbedbd860982 % UInt
  h = hash(parent(x), h)
  h = hash(coefficients(x), h)
  return xor(h, b)
end

###############################################################################
#
#   Important ideals and subalgebras
#
###############################################################################

@doc raw"""
    derived_algebra(L::LieAlgebra) -> LieAlgebraIdeal

Return the derived algebra of `L`, i.e. $[L, L]$.
"""
function derived_algebra(L::LieAlgebra)
  L_ideal = ideal(L)
  return bracket(L_ideal, L_ideal)
end

@doc raw"""
    center(L::LieAlgebra) -> LieAlgebraIdeal

Return the center of `L`, i.e. $\{x \in L \mid [x, L] = 0\}$
"""
function center(L::LieAlgebra)
  dim(L) == 0 && return ideal(L, [])

  mat = zero_matrix(coefficient_ring(L), dim(L), dim(L)^2)
  for (i, bi) in enumerate(basis(L))
    for (j, bj) in enumerate(basis(L))
      mat[i, ((j - 1) * dim(L) + 1):(j * dim(L))] = _matrix(bi * bj)
    end
  end

  c_basis = kernel(mat; side = :left)
  c_dim = nrows(c_basis)
  return ideal(L, [L(c_basis[i, :]) for i in 1:c_dim]; is_basis=true)
end

@doc raw"""
    centralizer(L::LieAlgebra, xs::AbstractVector{<:LieAlgebraElem}) -> LieSubalgebra

Return the centralizer of `xs` in `L`, i.e. $\{y \in L \mid [x, y] = 0 \forall x \in xs\}$.
"""
function centralizer(L::LieAlgebra, xs::AbstractVector{<:LieAlgebraElem})
  @req all(x -> parent(x) === L, xs) "Incompatible Lie algebras."

  mat = zero_matrix(coefficient_ring(L), dim(L), dim(L) * length(xs))
  for (i, bi) in enumerate(basis(L))
    for (j, xj) in enumerate(xs)
      mat[i, ((j - 1) * dim(L) + 1):(j * dim(L))] = _matrix(bracket(bi, xj))
    end
  end

  c_basis = kernel(mat; side = :left)
  c_dim = nrows(c_basis)
  return sub(L, [L(c_basis[i, :]) for i in 1:c_dim]; is_basis=true)
end

@doc raw"""
    centralizer(L::LieAlgebra, x::LieAlgebraElem) -> LieSubalgebra

Return the centralizer of `x` in `L`, i.e. $\{y \in L \mid [x, y] = 0\}$.
"""
function centralizer(L::LieAlgebra, x::LieAlgebraElem)
  return centralizer(L, [x])
end

###############################################################################
#
#   Derived and central series
#
###############################################################################

@doc raw"""
    derived_series(L::LieAlgebra) -> Vector{LieAlgebraIdeal}

Return the derived series of `L`, i.e. the sequence of ideals 
$L^{(0)} = L$, $L^{(i + 1)} = [L^{(i)}, L^{(i)}]$.
"""
function derived_series(L::LieAlgebra)
  curr = ideal(L)
  series = [curr]
  while true
    next = bracket(curr, curr)
    dim(next) == dim(curr) && break
    push!(series, next)
    curr = next
  end
  return series
end

@doc raw"""
    lower_central_series(L::LieAlgebra) -> Vector{LieAlgebraIdeal}

Return the lower central series of `L`, i.e. the sequence of ideals
$L^{(0)} = L$, $L^{(i + 1)} = [L, L^{(i)}]$.
"""
function lower_central_series(L::LieAlgebra)
  curr = ideal(L)
  series = [curr]
  while true
    next = bracket(L, curr)
    dim(next) == dim(curr) && break
    push!(series, next)
    curr = next
  end
  return series
end

###############################################################################
#
#   Properties
#
###############################################################################

@doc raw"""
    is_abelian(L::LieAlgebra) -> Bool

Return `true` if `L` is abelian, i.e. $[L, L] = 0$.
"""
@attr Bool function is_abelian(L::LieAlgebra)
  b = basis(L)
  n = length(b)
  return all(iszero, b[i] * b[j] for i in 1:n for j in (i + 1):n)
end

@doc raw"""
    is_nilpotent(L::LieAlgebra) -> Bool

Return `true` if `L` is nilpotent, i.e. the lower central series of `L` terminates in $0$.
"""
@attr Bool function is_nilpotent(L::LieAlgebra)
  return dim(lower_central_series(L)[end]) == 0
end

@doc raw"""
    is_perfect(L::LieAlgebra) -> Bool

Return `true` if `L` is perfect, i.e. $[L, L] = L$.
"""
@attr Bool function is_perfect(L::LieAlgebra)
  return dim(derived_algebra(L)) == dim(L)
end

@doc raw"""
    is_simple(L::LieAlgebra) -> Bool

Return `true` if `L` is simple, i.e. `L` is not abelian and has no non-trivial ideals.

!!! warning
    This function is not implemented yet.
"""
@attr Bool function is_simple(L::LieAlgebra)
  is_abelian(L) && return false
  !is_perfect(L) && return false
  error("Not implemented.") # TODO
end

@doc raw"""
    is_solvable(L::LieAlgebra) -> Bool

Return `true` if `L` is solvable, i.e. the derived series of `L` terminates in $0$.
"""
@attr Bool function is_solvable(L::LieAlgebra)
  return dim(derived_series(L)[end]) == 0
end

###############################################################################
#
#   Universal enveloping algebra
#
###############################################################################

@doc raw"""
    universal_enveloping_algebra(L::LieAlgebra; ordering::Symbol=:lex) -> PBWAlgRing, Map

Return the universal enveloping algebra `U(L)` of `L` with the given monomial ordering,
together with a map from `L` into the filtered component of degree 1 of `U(L)`.
"""
function universal_enveloping_algebra(L::LieAlgebra; ordering::Symbol=:lex)
  R, gensR = polynomial_ring(coefficient_ring(L), symbols(L))
  n = dim(L)
  b = basis(L)

  to_R(x::LieAlgebraElem) =
    sum(c * g for (c, g) in zip(coefficients(x), gensR); init=zero(R))

  rel = strictly_upper_triangular_matrix([
    to_R(b[i]) * to_R(b[j]) - to_R(b[i] * b[j]) for i in 1:(n - 1) for j in (i + 1):n
  ])
  U, gensU = pbw_algebra(R, rel, monomial_ordering(R, ordering); check=true)

  L_to_U = MapFromFunc(
    L, U, function (x::LieAlgebraElem)
      sum(c * g for (c, g) in zip(coefficients(x), gensU); init=zero(U))
    end
  )
  return U, L_to_U
end

###############################################################################
#
#   Constructor
#
###############################################################################

@doc raw"""
    lie_algebra(gapL::GAP.GapObj, s::Vector{<:VarName}; cached::Bool) -> LieAlgebra{elem_type(R)}

Construct a Lie algebra isomorphic to the GAP Lie algebra `gapL`. Its basis element are named by `s`,
or by `x_i` by default.
We require `gapL` to be a finite-dimensional GAP Lie algebra. The return type is dependent on
properties of `gapL`, in particular, whether GAP knows about a matrix representation.

If `cached` is `true`, the constructed Lie algebra is cached.
"""
function lie_algebra(
  gapL::GAP.GapObj,
  s::Vector{<:VarName}=[Symbol("x_$i") for i in 1:GAPWrap.Dimension(gapL)];
  cached::Bool=true,
)
  @req GAPWrap.IsLieAlgebra(gapL) "gapL must be a Lie algebra."
  if GAPWrap.IsFiniteDimensional(gapL)
    if GAPWrap.IsLieObjectCollection(gapL)
      return codomain(_iso_gap_oscar_linear_lie_algebra(gapL, s; cached))
    else
      return codomain(_iso_gap_oscar_abstract_lie_algebra(gapL, s; cached))
    end
  end
  error("Not implemented.")
end
