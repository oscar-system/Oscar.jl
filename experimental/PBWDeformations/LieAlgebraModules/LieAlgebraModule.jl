#################################################
#
# Abstract parent type
#
#################################################

abstract type LieAlgebraModule{C<:RingElement} end

abstract type LieAlgebraModuleElem{C<:RingElement} end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_ring(v::LieAlgebraModuleElem{C}) where {C<:RingElement} = base_ring(parent(v))

ngens(L::LieAlgebraModule{C}) where {C<:RingElement} = dim(L)

gens(L::LieAlgebraModule{C}) where {C<:RingElement} = basis(L)

gen(L::LieAlgebraModule{C}, i::Int) where {C<:RingElement} = basis(L, i)

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

@doc Markdown.doc"""
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
  return V1 === V2
end

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

function action(
  x::LieAlgebraElem{C}, v::ElemT
) where {ElemT<:LieAlgebraModuleElem{C}} where {C<:RingElement}
  @req parent(x) == base_liealgebra(parent(v)) "Incompatible Lie algebras."

  cx = Generic._matrix(x)

  return parent(v)(
    sum(
      cx[i] *
      Generic._matrix(v) *
      transpose(transformation_matrix_by_basisindex(parent(v), i)) for
      i in 1:dim(parent(x)) if !iszero(cx[i]);
      init=zero_matrix(base_ring(parent(v)), 1, dim(parent(v)))::dense_matrix_type(C),
    ), # equivalent to (x * v^T)^T, since we work with row vectors
  )::ElemT
end

function Base.:*(x::LieAlgebraElem{C}, v::LieAlgebraModuleElem{C}) where {C<:RingElement}
  return action(x, v)
end

function transformation_matrix_by_basisindex(
  V::LieAlgebraModule{C}, _::Int
) where {C<:RingElement}
  error("Not implemented for $(typeof(V))")
end
