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
  return iszero(Generic._matrix(x))
end

@inline function Generic._matrix(x::LieAlgebraElem{C}) where {C<:RingElement}
  return (x.mat)::dense_matrix_type(C)
end

@doc raw"""
    getindex(x::LieAlgebraElem{C}, i::Int) where C <: RingElement

Return the $i$-th coefficient of the module element $x$.
"""
function getindex(x::LieAlgebraElem{C}, i::Int) where {C<:RingElement}
  return Generic._matrix(x)[1, i]
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
  for (i, c) in enumerate(Generic._matrix(v))
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
  return Generic._matrix(x1) == Generic._matrix(x2)
end

function Base.hash(x::LieAlgebraElem{C}, h::UInt) where {C<:RingElement}
  b = 0x6724cbedbd860982 % UInt
  return xor(hash(Generic._matrix(x), hash(parent(x), h)), b)
end

###############################################################################
#
#   Attribute accessors
#
###############################################################################

function _gap_object(L::LieAlgebra{C}) where {C<:RingElement}
  # later change to storing an isomorphism instead
  get_attribute!(L, :gap_object) do
    gap_lie_algebra_by_struct_consts(L)
  end
end

function _set_gap_object!(L::LieAlgebra{C}, gapL::GAP.Obj) where {C<:RingElement}
  set_attribute!(L, :gap_object => gapL)
end

###############################################################################
#
#   Constructor
#
###############################################################################

function general_linear_lie_algebra(R::Ring, n::Int)
  basis = [(b = zero_matrix(R, n, n); b[i, j] = 1; b) for i in 1:n for j in 1:n]
  s = ["x_$(i)_$(j)" for i in 1:n for j in 1:n]
  L = lie_algebra(R, n, basis, s)
  set_attribute!(L, :type => :general_linear)
  return L
end

function special_linear_lie_algebra(R::Ring, n::Int)
  basis_e = [(b = zero_matrix(R, n, n); b[i, j] = 1; b) for i in 1:n for j in (i + 1):n]
  basis_f = [(b = zero_matrix(R, n, n); b[j, i] = 1; b) for i in 1:n for j in (i + 1):n]
  basis_h = [
    (b = zero_matrix(R, n, n); b[i, i] = 1; b[i + 1, i + 1] = -1; b) for i in 1:(n - 1)
  ]
  s_e = ["e_$(i)_$(j)" for i in 1:n for j in (i + 1):n]
  s_f = ["f_$(i)_$(j)" for i in 1:n for j in (i + 1):n]
  s_h = ["h_$(i)" for i in 1:(n - 1)]
  L = lie_algebra(R, n, [basis_e; basis_f; basis_h], [s_e; s_f; s_h])
  set_attribute!(L, :type => :special_linear)
  return L
end

function special_orthogonal_lie_algebra(R::Ring, n::Int)
  basis = [
    (b = zero_matrix(R, n, n); b[i, j] = 1; b[j, i] = -1; b) for i in 1:n for j in (i + 1):n
  ]
  s = ["x_$(i)_$(j)" for i in 1:n for j in (i + 1):n]
  L = lie_algebra(R, n, basis, s)
  set_attribute!(L, :type => :special_orthogonal)
  return L
end
