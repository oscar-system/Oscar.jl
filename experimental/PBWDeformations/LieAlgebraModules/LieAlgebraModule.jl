#################################################
#
# Abstract parent type
#
#################################################

abstract type LieAlgebraModule{C<:RingElement} <: FPModule{C} end

abstract type LieAlgebraModuleElem{C<:RingElement} <: FPModuleElem{C} end

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

@inline function Generic._matrix(v::LieAlgebraModuleElem{C}) where {C<:RingElement}
  return (v.mat)::dense_matrix_type(C)
end

function Generic.rels(_::LieAlgebraModule{C}) where {C<:RingElement}
  # there are no relations in a vector space
  return Vector{dense_matrix_type(C)}(undef, 0)
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
  mat = zero_matrix(base_ring(V), 1, dim(V))
  return elem_type(V)(V, mat)
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

# Vector space operations get inherited from FPModule

###############################################################################
#
#   Comparison functions
#
###############################################################################

# Overwrite the equality of FPModule to be used for CacheDicts
function Base.:(==)(V1::LieAlgebraModule{C}, V2::LieAlgebraModule{C}) where {C<:RingElement}
  return V1 === V2
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
