@attributes mutable struct LieAlgebraIdeal{C<:RingElement,LieT<:LieAlgebraElem{C}}
  base_lie_algebra::LieAlgebra{C}
  gens::Vector{LieT}
  basis::Vector{LieT}
  basis_matrix::MatElem{C}

  function LieAlgebraIdeal{C,LieT}(
    L::LieAlgebra{C}, gens::Vector{LieT}; is_basis::Bool=false
  ) where {C<:RingElement,LieT<:LieAlgebraElem{C}}
    @req all(g -> parent(g) === L, gens) "Parent mismatch."
    L::parent_type(LieT)
    if is_basis
      basis = gens
      basis_matrix = matrix([coefficients(g) for g in gens])
      return new{C,LieT}(L, gens, basis, basis_matrix)
    else
      return new{C,LieT}(L, gens)
    end
  end

  function LieAlgebraIdeal{C,LieT}(
    L::LieAlgebra{C}, gens::Vector; kwargs...
  ) where {C<:RingElement,LieT<:LieAlgebraElem{C}}
    return LieAlgebraIdeal{C,LieT}(L, Vector{LieT}(map(L, gens)); kwargs...)
  end
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_lie_algebra(
  I::LieAlgebraIdeal{C,LieT}
) where {C<:RingElement,LieT<:LieAlgebraElem{C}} = I.base_lie_algebra::parent_type(LieT)

function gens(I::LieAlgebraIdeal)
  return I.gens
end

function ngens(I::LieAlgebraIdeal)
  return length(gens(I))
end

function basis_matrix(I::LieAlgebraIdeal)
  if !isdefined(I, :basis_matrix)
    rank, mat = rref(
      matrix(
        [
          [coefficients(g) for g in gens(I) if !iszero(g)]
          [
            coefficients(x * g) for x in basis(base_lie_algebra(I)) for
            g in gens(I) if !iszero(x * g)
          ]
          [coefficients(zero(base_lie_algebra(I)))]
        ],
      ),
    )
    I.basis_matrix = mat[1:rank, :]
  end
  return I.basis_matrix
end

function has_basis(I::LieAlgebraIdeal)
  return isdefined(I, :basis)
end

function basis(I::LieAlgebraIdeal)
  if !has_basis(I)
    mat = basis_matrix(I)
    I.basis = [base_lie_algebra(I)(mat[i, :]) for i in 1:nrows(mat)]
  end
  return I.basis
end

dim(I::LieAlgebraIdeal) = length(basis(I))

###############################################################################
#
#   Comparisons
#
###############################################################################

function Base.issubset(
  I1::LieAlgebraIdeal{C,LieT}, I2::LieAlgebraIdeal{C,LieT}
) where {C<:RingElement,LieT<:LieAlgebraElem{C}}
  @req base_lie_algebra(I1) === base_lie_algebra(I2) "Incompatible Lie algebras."
  return all(in(I2), gens(I1))
end

function Base.:(==)(
  I1::LieAlgebraIdeal{C,LieT}, I2::LieAlgebraIdeal{C,LieT}
) where {C<:RingElement,LieT<:LieAlgebraElem{C}}
  base_lie_algebra(I1) === base_lie_algebra(I2) || return false
  gens(I1) == gens(I2) && return true
  return basis_matrix(I1) == basis_matrix(I2)
end

function Base.hash(I::LieAlgebraIdeal, h::UInt)
  b = 0xd745d725d0b66431 % UInt
  h = hash(base_lie_algebra(I), h)
  h = hash(gens(I), h)
  h = hash(basis_matrix(I), h)
  return xor(h, b)
end

###############################################################################
#
#   Ideal membership
#
###############################################################################

function Base.in(x::LieAlgebraElem, I::LieAlgebraIdeal)
  return can_solve(basis_matrix(I), _matrix(x); side=:left)
end

###############################################################################
#
#   Constructions
#
###############################################################################

function Base.:+(
  I1::LieAlgebraIdeal{C,LieT}, I2::LieAlgebraIdeal{C,LieT}
) where {C<:RingElement,LieT<:LieAlgebraElem{C}}
  @req base_lie_algebra(I1) === base_lie_algebra(I2) "Incompatible Lie algebras."
  return ideal(base_lie_algebra(I1), [gens(I1); gens(I2)])
end

function bracket(
  I1::LieAlgebraIdeal{C,LieT}, I2::LieAlgebraIdeal{C,LieT}
) where {C<:RingElement,LieT<:LieAlgebraElem{C}}
  @req base_lie_algebra(I1) === base_lie_algebra(I2) "Incompatible Lie algebras."
  return ideal(base_lie_algebra(I1), [x * y for x in gens(I1) for y in gens(I2)])
end

###############################################################################
#
#   Constructor
#
###############################################################################

function ideal(L::LieAlgebra, gens::Vector; is_basis::Bool=false)
  return LieAlgebraIdeal{elem_type(coefficient_ring(L)),elem_type(L)}(L, gens; is_basis)
end

function ideal(L::LieAlgebra{C}, gen::LieAlgebraElem{C}) where {C<:RingElement}
  return ideal(L, [gen])
end

function ideal(L::LieAlgebra)
  return ideal(L, basis(L); is_basis=true)
end
