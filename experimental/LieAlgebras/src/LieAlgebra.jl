#################################################
#
# Abstract parent type
#
#################################################

abstract type LieAlgebra{C<:RingElement} end

abstract type LieAlgebraElem{C<:RingElement} end

# To be implemented by subtypes:
# parent_type(::Type{MyLieAlgebraElem{C}})
# elem_type(::Type{MyLieAlgebra{C}})
# parent(x::MyLieAlgebraElem{C})
# base_ring(L::MyLieAlgebra{C})
# dim(L::MyLieAlgebra{C})
# Base.show(io::IO, x::MyLieAlgebra{C})
# symbols(L::MyLieAlgebra{C})
# bracket(x::MyLieAlgebraElem{C}, y::MyLieAlgebraElem{C})

###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_ring(x::LieAlgebraElem{C}) where {C<:RingElement} = base_ring(parent(x))

ngens(L::LieAlgebra{C}) where {C<:RingElement} = dim(L)

gens(L::LieAlgebra{C}) where {C<:RingElement} = basis(L)

gen(L::LieAlgebra{C}, i::Int) where {C<:RingElement} = basis(L, i)

@doc raw"""
    dim(L::LieAlgebra{C}) -> Int

Return the dimension of the Lie algebra `L`.
"""
dim(_::LieAlgebra{C}) where {C<:RingElement} = error("Should be implemented by subtypes.")

@doc raw"""
    basis(L::LieAlgebra{C}) -> Vector{LieAlgebraElem{C}}

Return a basis of the Lie algebra `L`.
"""
basis(L::LieAlgebra{C}) where {C<:RingElement} =
  [basis(L, i)::elem_type(L) for i in 1:dim(L)]

@doc raw"""
    basis(L::LieAlgebra{C}, i::Int) -> LieAlgebraElem{C}

Return the `i`-th basis element of the Lie algebra `L`.
"""
function basis(L::LieAlgebra{C}, i::Int) where {C<:RingElement}
  R = base_ring(L)
  return L([(j == i ? one(R) : zero(R)) for j in 1:dim(L)])
end

@doc raw"""
    zero(L::LieAlgebra{C}) -> LieAlgebraElem{C}

Return the zero element of the Lie algebra `L`.
"""
function zero(L::LieAlgebra{C}) where {C<:RingElement}
  mat = zero_matrix(base_ring(L), 1, dim(L))
  return elem_type(L)(L, mat)
end

@doc raw"""
    iszero(x::LieAlgebraElem{C}) -> Bool

Check whether the Lie algebra element `x` is zero.
"""
function iszero(x::LieAlgebraElem{C}) where {C<:RingElement}
  return iszero(coefficients(x))
end

@inline function _matrix(x::LieAlgebraElem{C}) where {C<:RingElement}
  return (x.mat)::dense_matrix_type(C)
end

@doc raw"""
    coefficients(x::LieAlgebraElem{C}) -> Vector{C}

Return the coefficients of `x` with respect to [`basis(::LieAlgebra)`](@ref).
"""
function coefficients(x::LieAlgebraElem{C}) where {C<:RingElement}
  return collect(_matrix(x))[1, :]
end

@doc raw"""
    coeff(x::LieAlgebraElem{C}, i::Int) -> C

Return the `i`-th coefficient of `x` with respect to [`basis(::LieAlgebra)`](@ref).
"""
function coeff(x::LieAlgebraElem{C}, i::Int) where {C<:RingElement}
  return _matrix(x)[1, i]
end

@doc raw"""
    getindex(x::LieAlgebraElem{C}, i::Int) -> C

Return the `i`-th coefficient of `x` with respect to [`basis(::LieAlgebra)`](@ref).
"""
function getindex(x::LieAlgebraElem{C}, i::Int) where {C<:RingElement}
  return coeff(x, i)
end

function Base.deepcopy_internal(x::LieAlgebraElem{C}, dict::IdDict) where {C<:RingElement}
  return parent(x)(deepcopy_internal(_matrix(x), dict))
end

function check_parent(x1::LieAlgebraElem{C}, x2::LieAlgebraElem{C}) where {C<:RingElement}
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
symbols(_::LieAlgebra{C}) where {C<:RingElement} =
  error("Should be implemented by subtypes.")

function expressify(
  v::LieAlgebraElem{C}, s=symbols(parent(v)); context=nothing
) where {C<:RingElement}
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
function (L::LieAlgebra{C})() where {C<:RingElement}
  return zero(L)
end

@doc raw"""
    (L::LieAlgebra{C})(v::Vector{Int}) -> LieAlgebraElem{C}

Return the element of `L` with coefficent vector `v`.
Fail, if `Int` cannot be coerced into the base ring of `L`.
"""
function (L::LieAlgebra{C})(v::Vector{Int}) where {C<:RingElement}
  return L(base_ring(L).(v))
end

@doc raw"""
    (L::LieAlgebra{C})(v::Vector{C}) -> LieAlgebraElem{C}

Return the element of `L` with coefficent vector `v`.
"""
function (L::LieAlgebra{C})(v::Vector{C}) where {C<:RingElement}
  @req length(v) == dim(L) "Length of vector does not match dimension."
  mat = matrix(base_ring(L), 1, length(v), v)
  return elem_type(L)(L, mat)
end

@doc raw"""
    (L::LieAlgebra{C})(mat::MatElem{C}) -> LieAlgebraElem{C}

Return the element of `L` with coefficient vector equivalent to
the $1 \times \dim(L)$ matrix `mat`.
"""
function (L::LieAlgebra{C})(mat::MatElem{C}) where {C<:RingElement}
  @req size(mat) == (1, dim(L)) "Invalid matrix dimensions."
  return elem_type(L)(L, mat)
end

@doc raw"""
    (L::LieAlgebra{C})(v::SRow{C}) -> LieAlgebraElem{C}

Return the element of `L` with coefficent vector `v`.
"""
function (L::LieAlgebra{C})(v::SRow{C}) where {C<:RingElement}
  mat = dense_row(v, dim(L))
  return elem_type(L)(L, mat)
end

@doc raw"""
    (L::LieAlgebra{C})(x::LieAlgebraElem{C}) -> LieAlgebraElem{C}

Return `x`. Fails if `x` is not an element of `L`.
"""
function (L::LieAlgebra{C})(x::LieAlgebraElem{C}) where {C<:RingElement}
  @req L == parent(x) "Incompatible modules."
  return x
end

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

function Base.:-(x::LieAlgebraElem{C}) where {C<:RingElement}
  return parent(x)(-_matrix(x))
end

function Base.:+(x1::LieAlgebraElem{C}, x2::LieAlgebraElem{C}) where {C<:RingElement}
  check_parent(x1, x2)
  return parent(x1)(_matrix(x1) + _matrix(x2))
end

function Base.:-(x1::LieAlgebraElem{C}, x2::LieAlgebraElem{C}) where {C<:RingElement}
  check_parent(x1, x2)
  return parent(x1)(_matrix(x1) - _matrix(x2))
end

function Base.:*(x::LieAlgebraElem{C}, c::C) where {C<:RingElem}
  base_ring(x) != parent(c) && error("Incompatible rings.")
  return parent(x)(_matrix(x) * c)
end

function Base.:*(
  x::LieAlgebraElem{C}, c::U
) where {C<:RingElement,U<:Union{Rational,Integer}}
  return parent(x)(_matrix(x) * c)
end

function Base.:*(c::C, x::LieAlgebraElem{C}) where {C<:RingElem}
  base_ring(x) != parent(c) && error("Incompatible rings.")
  return parent(x)(c * _matrix(x))
end

function Base.:*(
  c::U, x::LieAlgebraElem{C}
) where {C<:RingElement,U<:Union{Rational,Integer}}
  return parent(x)(c * _matrix(x))
end

function Base.:*(x::LieAlgebraElem{C}, y::LieAlgebraElem{C}) where {C<:RingElement}
  return bracket(x, y)
end

###############################################################################
#
#   Comparison functions
#
###############################################################################

function Base.:(==)(x1::LieAlgebraElem{C}, x2::LieAlgebraElem{C}) where {C<:RingElement}
  check_parent(x1, x2)
  return coefficients(x1) == coefficients(x2)
end

function Base.hash(x::LieAlgebraElem{C}, h::UInt) where {C<:RingElement}
  b = 0x6724cbedbd860982 % UInt
  h = hash(parent(x), h)
  h = hash(coefficients(x), h)
  return xor(h, b)
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
