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

basis(L::LieAlgebra{C}) where {C<:RingElement} =
  [basis(L, i)::elem_type(L) for i in 1:dim(L)]

function basis(L::LieAlgebra{C}, i::Int) where {C<:RingElement}
  R = base_ring(L)
  return L([(j == i ? one(R) : zero(R)) for j in 1:dim(L)])
end

function zero(L::LieAlgebra{C}) where {C<:RingElement}
  mat = zero_matrix(base_ring(L), 1, dim(L))
  return elem_type(L)(L, mat)
end

function iszero(x::LieAlgebraElem{C}) where {C<:RingElement}
  return iszero(coefficients(x))
end

@inline function Generic._matrix(x::LieAlgebraElem{C}) where {C<:RingElement}
  return (x.mat)::dense_matrix_type(C)
end

function coefficients(x::LieAlgebraElem{C}) where {C<:RingElement}
  return collect(Generic._matrix(x))[1, :]
end

function coeff(x::LieAlgebraElem{C}, i::Int) where {C<:RingElement}
  return Generic._matrix(x)[1, i]
end

@doc raw"""
    getindex(x::LieAlgebraElem{C}, i::Int) where C <: RingElement

Return the $i$-th coefficient of the module element $x$.
"""
function getindex(x::LieAlgebraElem{C}, i::Int) where {C<:RingElement}
  return coeff(x, i)
end

function Base.deepcopy_internal(x::LieAlgebraElem{C}, dict::IdDict) where {C<:RingElement}
  return parent(x)(deepcopy_internal(Generic._matrix(x), dict))
end

function check_parent(x1::LieAlgebraElem{C}, x2::LieAlgebraElem{C}) where {C<:RingElement}
  parent(x1) !== parent(x2) && error("Incompatible Lie algebras.")
end

###############################################################################
#
#   String I/O
#
###############################################################################

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

function (L::LieAlgebra{C})() where {C<:RingElement}
  return zero(L)
end

function (L::LieAlgebra{C})(v::Vector{Int}) where {C<:RingElement}
  return L(base_ring(L).(v))
end

function (L::LieAlgebra{C})(v::Vector{C}) where {C<:RingElement}
  @req length(v) == dim(L) "Length of vector does not match dimension."
  mat = matrix(base_ring(L), 1, length(v), v)
  return elem_type(L)(L, mat)
end

function (L::LieAlgebra{C})(mat::MatElem{C}) where {C<:RingElement}
  @req size(mat) == (1, dim(L)) "Invalid matrix dimensions."
  return elem_type(L)(L, mat)
end

function (L::LieAlgebra{C})(v::SRow{C}) where {C<:RingElement}
  mat = dense_row(v, dim(L))
  return elem_type(L)(L, mat)
end

function (L::LieAlgebra{C})(v::LieAlgebraElem{C}) where {C<:RingElement}
  @req L == parent(v) "Incompatible modules."
  return v
end

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

function Base.:-(x::LieAlgebraElem{C}) where {C<:RingElement}
  return parent(x)(-Generic._matrix(x))
end

function Base.:+(x1::LieAlgebraElem{C}, x2::LieAlgebraElem{C}) where {C<:RingElement}
  check_parent(x1, x2)
  return parent(x1)(Generic._matrix(x1) + Generic._matrix(x2))
end

function Base.:-(x1::LieAlgebraElem{C}, x2::LieAlgebraElem{C}) where {C<:RingElement}
  check_parent(x1, x2)
  return parent(x1)(Generic._matrix(x1) - Generic._matrix(x2))
end

function Base.:*(x::LieAlgebraElem{C}, c::C) where {C<:RingElem}
  base_ring(x) != parent(c) && error("Incompatible rings.")
  return parent(x)(Generic._matrix(x) * c)
end

function Base.:*(
  x::LieAlgebraElem{C}, c::U
) where {C<:RingElement,U<:Union{Rational,Integer}}
  return parent(x)(Generic._matrix(x) * c)
end

function Base.:*(c::C, x::LieAlgebraElem{C}) where {C<:RingElem}
  base_ring(x) != parent(c) && error("Incompatible rings.")
  return parent(x)(c * Generic._matrix(x))
end

function Base.:*(
  c::U, x::LieAlgebraElem{C}
) where {C<:RingElement,U<:Union{Rational,Integer}}
  return parent(x)(c * Generic._matrix(x))
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
  return xor(hash(coefficients(x), hash(parent(x), h)), b)
end

###############################################################################
#
#   Isomorphism to GAP accessors
#
###############################################################################

# currently, the user should only use _iso_oscar_gap(LO) to get the isomorphism

function _get_iso_oscar_gap!(f, LO::LieAlgebra{C}) where {C<:RingElement}
  get_attribute!(f, LO, :iso_oscar_gap)::Map{typeof(LO),GAP.GapObj}
end

function _set_iso_oscar_gap!(
  LO::LieAlgebra{C}, iso::Map{<:LieAlgebra{C},GAP.GapObj}
) where {C<:RingElement}
  set_attribute!(LO, :iso_oscar_gap => iso)
end

###############################################################################
#
#   Constructor
#
###############################################################################

function lie_algebra(
  gapL::GAP.GapObj,
  s::Vector{<:VarName}=[Symbol("x_$i") for i in 1:GAPWrap.Dimension(gapL)];
  cached::Bool=true,
)
  @req GAP.Globals.IsLieAlgebra(gapL) "gapL must be a Lie algebra."
  if GAP.Globals.IsLieObjectCollection(gapL)
    return codomain(_iso_gap_oscar_linear_lie_algebra(gapL, s; cached))
  else
    return codomain(_iso_gap_oscar_abstract_lie_algebra(gapL, s; cached))
  end
end
