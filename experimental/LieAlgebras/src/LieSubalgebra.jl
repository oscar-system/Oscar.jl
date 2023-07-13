@attributes mutable struct LieSubalgebra{C<:RingElement,LieT<:LieAlgebraElem{C}}
  base_lie_algebra::LieAlgebra{C}
  gens::Vector{LieT}
  basis_elems::Vector{LieT}
  basis_matrix::MatElem{C}

  function LieSubalgebra{C,LieT}(
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
        for i in 1:nrows(basis_matrix)
          push!(todo, g * L(basis_matrix[i, :]))
        end
        basis_matrix = vcat(basis_matrix, _matrix(g))
        rank = rref!(basis_matrix)
        basis_matrix = basis_matrix[1:rank, :]
      end
      basis_elems = [L(basis_matrix[i, :]) for i in 1:nrows(basis_matrix)]
      return new{C,LieT}(L, gens, basis_elems, basis_matrix)
    end
  end

  function LieSubalgebra{C,LieT}(
    L::LieAlgebra{C}, gens::Vector; kwargs...
  ) where {C<:RingElement,LieT<:LieAlgebraElem{C}}
    return LieSubalgebra{C,LieT}(L, Vector{LieT}(map(L, gens)); kwargs...)
  end
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_lie_algebra(S::LieSubalgebra{C,LieT}) where {C<:RingElement,LieT<:LieAlgebraElem{C}} =
  S.base_lie_algebra::parent_type(LieT)

function gens(S::LieSubalgebra)
  return S.gens
end

function gen(S::LieSubalgebra, i::Int)
  return S.gens[i]
end

function ngens(S::LieSubalgebra)
  return length(gens(S))
end

function basis_matrix(
  S::LieSubalgebra{C,LieT}
) where {C<:RingElement,LieT<:LieAlgebraElem{C}}
  return S.basis_matrix::dense_matrix_type(C)
end

@doc raw"""
    basis(S::LieSubalgebra{C}) -> Vector{LieAlgebraElem{C}}

Return a basis of the Lie subalgebra `S`.
"""
function basis(S::LieSubalgebra)
  return S.basis_elems
end

@doc raw"""
    basis(S::LieSubalgebra{C}, i::Int) -> LieAlgebraElem{C}

Return the `i`-th basis element of the Lie subalgebra `S`.
"""
function basis(S::LieSubalgebra, i::Int)
  return S.basis_elems[i]
end

@doc raw"""
  dim(S::LieSubalgebra) -> Int

Return the dimension of the Lie subalgebra `S`.
"""
dim(S::LieSubalgebra) = length(basis(S))

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, ::MIME"text/plain", S::LieSubalgebra)
  io = pretty(io)
  println(io, LowercaseOff(), "Lie subalgebra")
  println(io, Indent(), "of dimension $(dim(S))", Dedent())
  if dim(S) != ngens(S)
    println(io, Indent(), "with $(ItemQuantity(ngens(S), "generator"))", Dedent())
  end
  print(io, "of ")
  print(io, Lowercase(), base_lie_algebra(S))
end

function Base.show(io::IO, S::LieSubalgebra)
  io = pretty(io)
  if get(io, :supercompact, false)
    print(io, LowercaseOff(), "Lie subalgebra")
  else
    print(io, LowercaseOff(), "Lie subalgebra of ", Lowercase())
    print(IOContext(io, :supercompact => true), base_lie_algebra(S))
  end
end

###############################################################################
#
#   Comparisons
#
###############################################################################

function Base.issubset(
  S1::LieSubalgebra{C,LieT}, S2::LieSubalgebra{C,LieT}
) where {C<:RingElement,LieT<:LieAlgebraElem{C}}
  @req base_lie_algebra(S1) === base_lie_algebra(S2) "Incompatible Lie algebras."
  return all(in(S2), gens(S1))
end

function Base.:(==)(
  S1::LieSubalgebra{C,LieT}, S2::LieSubalgebra{C,LieT}
) where {C<:RingElement,LieT<:LieAlgebraElem{C}}
  base_lie_algebra(S1) === base_lie_algebra(S2) || return false
  gens(S1) == gens(S2) && return true
  return basis_matrix(S1) == basis_matrix(S2)
end

function Base.hash(S::LieSubalgebra, h::UInt)
  b = 0x712b12e671184d1a % UInt
  h = hash(base_lie_algebra(S), h)
  h = hash(gens(S), h)
  h = hash(basis_matrix(S), h)
  return xor(h, b)
end

###############################################################################
#
#   Subalgebra membership
#
###############################################################################

@doc raw"""
    in(x::LieAlgebraElem, S::LieSubalgebra) -> Bool

Return `true` if `x` is in the Lie subalgebra `S`, `false` otherwise.
"""
function Base.in(x::LieAlgebraElem, S::LieSubalgebra)
  return can_solve(basis_matrix(S), _matrix(x); side=:left)
end

###############################################################################
#
#   Constructions
#
###############################################################################

@doc raw"""
    bracket(S1::LieSubalgebra, S2::LieSubalgebra) -> LieAlgebraIdeal

Return $[S_1, S_2]$.
"""
function bracket(
  S1::LieSubalgebra{C,LieT}, S2::LieSubalgebra{C,LieT}
) where {C<:RingElement,LieT<:LieAlgebraElem{C}}
  @req base_lie_algebra(S1) === base_lie_algebra(S2) "Incompatible Lie algebras."
  return ideal(base_lie_algebra(S1), [x * y for x in gens(S1) for y in gens(S2)])
end

###############################################################################
#
#   Important ideals and subalgebras
#
###############################################################################

@doc raw"""
    normalizer(L::LieAlgebra, S::LieSubalgebra) -> LieSubalgebra

Return the normalizer of `S` in `L`, i.e. $\{x \in L \mid [x, S] \subseteq S\}$.
"""
function normalizer(L::LieAlgebra, S::LieSubalgebra)
  @req base_lie_algebra(S) == L "Incompatible Lie algebras."

  mat = zero_matrix(coefficient_ring(L), dim(L) + dim(S)^2, dim(L) * dim(S))
  for (i, bi) in enumerate(basis(L))
    for (j, sj) in enumerate(basis(S))
      mat[i, ((j - 1) * dim(L) + 1):(j * dim(L))] = _matrix(bracket(bi, sj))
    end
  end
  for i in 1:dim(S)
    mat[(dim(L) + (i - 1) * dim(S) + 1):(dim(L) + i * dim(S)), ((i - 1) * dim(L) + 1):(i * dim(L))] = basis_matrix(
      S
    )
  end
  show(stdout, MIME"text/plain"(), mat)
  println("")
  sol_dim, sol = left_kernel(mat)
  sol = sol[:, 1:dim(L)]
  c_dim, c_basis = rref(sol)
  println(c_dim)
  show(stdout, MIME"text/plain"(), c_basis)
  println("")
  return sub(L, [L(c_basis[i, :]) for i in 1:c_dim]; is_basis=true)
end

@doc raw"""
    centralizer(L::LieAlgebra, S::LieSubalgebra) -> LieSubalgebra
  
Return the centralizer of `S` in `L`, i.e. $\{x \in L \mid [x, S] = 0\}$.
"""
function centralizer(L::LieAlgebra, S::LieSubalgebra)
  return centralizer(L, basis(S))
end

###############################################################################
#
#   Properties
#
###############################################################################

@doc raw"""
    is_self_normalizing(S::LieSubalgebra) -> Bool

Return `true` if `S` is self-normalizing, i.e. if its normalizer is `S`.
"""
function is_self_normalizing(S::LieSubalgebra)
  return normalizer(base_lie_algebra(S), S) == S
end

###############################################################################
#
#   Conversion
#
###############################################################################

@doc raw"""
    lie_algebra(S::LieSubalgebra) -> LieAlgebra

Return `S` as a Lie algebra.
"""
function lie_algebra(S::LieSubalgebra)
  return lie_algebra(basis(S)) #, embedding_hom   # TODO
end

###############################################################################
#
#   Constructor
#
###############################################################################

@doc raw"""
    sub(L::LieAlgebra, gens::Vector{LieAlgebraElem}; is_basis::Bool=false) -> LieSubalgebra

Return the smallest Lie subalgebra of `L` containing `gens`.
If `is_basis` is `true`, then `gens` is assumed to be a basis of the subalgebra.
"""
function sub(L::LieAlgebra, gens::Vector; is_basis::Bool=false)
  return LieSubalgebra{elem_type(coefficient_ring(L)),elem_type(L)}(L, gens; is_basis)
end

@doc raw"""
    sub(L::LieAlgebra, gen::LieAlgebraElem) -> LieSubalgebra

Return the smallest Lie subalgebra of `L` containing `gen`.
"""
function sub(L::LieAlgebra{C}, gen::LieAlgebraElem{C}) where {C<:RingElement}
  return sub(L, [gen])
end

@doc raw"""
    sub(L::LieAlgebra) -> LieSubalgebra

Return `L` as a Lie subalgebra of itself.
"""
function sub(L::LieAlgebra)
  return sub(L, basis(L); is_basis=true)
end
