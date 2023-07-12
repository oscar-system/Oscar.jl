@attributes mutable struct LieAlgebraIdeal{C<:RingElement,LieT<:LieAlgebraElem{C}}
  base_lie_algebra::LieAlgebra{C}
  gens::Vector{LieT}
  basis_elems::Vector{LieT}
  basis_matrix::MatElem{C}

  function LieAlgebraIdeal{C,LieT}(
    L::LieAlgebra{C}, gens::Vector{LieT}; is_basis::Bool=false
  ) where {C<:RingElement,LieT<:LieAlgebraElem{C}}
    @req all(g -> parent(g) === L, gens) "Parent mismatch."
    L::parent_type(LieT)
    if is_basis
      basis_elems = gens
      basis_matrix = if length(gens) == 0
        matrix(coefficient_ring(L), 0, dim(L), C[])
      else
        matrix(coefficient_ring(L), [coefficients(g) for g in gens])
      end
      return new{C,LieT}(L, gens, basis_elems, basis_matrix)
    else
      basis_matrix = matrix(coefficient_ring(L), 0, dim(L), C[])
      todo = copy(gens)
      while !isempty(todo)
        g = pop!(todo)
        can_solve(basis_matrix, _matrix(g); side=:left) && continue
        basis_matrix = vcat(basis_matrix, _matrix(g))
        for b in basis(L)
          push!(todo, b * g)
        end
        rank = rref!(basis_matrix)
        basis_matrix = basis_matrix[1:rank, :]
      end
      basis_elems = [L(basis_matrix[i, :]) for i in 1:nrows(basis_matrix)]
      return new{C,LieT}(L, gens, basis_elems, basis_matrix)
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

function basis_matrix(
  I::LieAlgebraIdeal{C,LieT}
) where {C<:RingElement,LieT<:LieAlgebraElem{C}}
  return I.basis_matrix::dense_matrix_type(C)
end

function basis(I::LieAlgebraIdeal)
  return I.basis_elems
end

dim(I::LieAlgebraIdeal) = length(basis(I))

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, ::MIME"text/plain", I::LieAlgebraIdeal)
  io = pretty(io)
  println(io, LowercaseOff(), "Lie algebra ideal")
  println(io, Indent(), "of dimension $(dim(I))", Dedent())
  if dim(I) != ngens(I)
    println(io, Indent(), "with $(ngens(I)) generators", Dedent())
  end
  print(io, "over ")
  print(io, Lowercase(), base_lie_algebra(I))
end

function Base.show(io::IO, I::LieAlgebraIdeal)
  io = pretty(io)
  if get(io, :supercompact, false)
    print(io, LowercaseOff(), "Lie algebra ideal")
  else
    print(io, LowercaseOff(), "Lie algebra ideal over ", Lowercase())
    print(IOContext(io, :supercompact => true), base_lie_algebra(I))
  end
end

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
