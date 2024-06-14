@attributes mutable struct LieAlgebraIdeal{C<:FieldElem,LieT<:LieAlgebraElem{C}}
  base_lie_algebra::LieAlgebra{C}
  gens::Vector{LieT}
  basis_elems::Vector{LieT}
  basis_matrix::MatElem{C}

  function LieAlgebraIdeal{C,LieT}(
    L::LieAlgebra{C}, gens::Vector{LieT}; is_basis::Bool=false
  ) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
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
      gens = unique(g for g in gens if !iszero(g))
      left = copy(gens)
      while !isempty(left)
        g = pop!(left)
        can_solve(basis_matrix, _matrix(g); side=:left) && continue
        for b in basis(L)
          push!(left, b * g)
        end
        basis_matrix = vcat(basis_matrix, _matrix(g))
        rank = rref!(basis_matrix)
        basis_matrix = basis_matrix[1:rank, :]
      end
      basis_elems = [L(basis_matrix[i, :]) for i in 1:nrows(basis_matrix)]
      return new{C,LieT}(L, gens, basis_elems, basis_matrix)
    end
  end

  function LieAlgebraIdeal{C,LieT}(
    L::LieAlgebra{C}, gens::Vector; kwargs...
  ) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
    return LieAlgebraIdeal{C,LieT}(L, Vector{LieT}(map(L, gens)); kwargs...)
  end
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_lie_algebra(I::LieAlgebraIdeal{C,LieT}) where {C<:FieldElem,LieT<:LieAlgebraElem{C}} =
  I.base_lie_algebra::parent_type(LieT)

function gens(I::LieAlgebraIdeal)
  return I.gens
end

function gen(I::LieAlgebraIdeal, i::Int)
  return I.gens[i]
end

function number_of_generators(I::LieAlgebraIdeal)
  return length(gens(I))
end

function basis_matrix(
  I::LieAlgebraIdeal{C,LieT}
) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
  return I.basis_matrix::dense_matrix_type(C)
end

@doc raw"""
    basis(I::LieAlgebraIdeal) -> Vector{LieAlgebraElem}

Return the basis of the ideal `I`.
"""
function basis(I::LieAlgebraIdeal)
  return I.basis_elems
end

@doc raw"""
    basis(I::LieAlgebraIdeal, i::Int) -> LieAlgebraElem

Return the `i`-th basis element of the ideal `I`.
"""
function basis(I::LieAlgebraIdeal, i::Int)
  return I.basis_elems[i]
end

@doc raw"""
    dim(I::LieAlgebraIdeal) -> Int

Return the dimension of the ideal `I`.
"""
dim(I::LieAlgebraIdeal) = length(basis(I))

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, mime::MIME"text/plain", I::LieAlgebraIdeal)
  @show_name(io, I)
  @show_special(io, mime, I)
  io = pretty(io)
  println(io, LowercaseOff(), "Lie algebra ideal")
  println(io, Indent(), "of dimension $(dim(I))", Dedent())
  if dim(I) != ngens(I)
    println(io, Indent(), "with $(ItemQuantity(ngens(I), "generator"))", Dedent())
  end
  print(io, "over ")
  print(io, Lowercase(), base_lie_algebra(I))
end

function Base.show(io::IO, I::LieAlgebraIdeal)
  @show_name(io, I)
  @show_special(io, I)
  io = pretty(io)
  if is_terse(io)
    print(io, LowercaseOff(), "Lie algebra ideal")
  else
    print(io, LowercaseOff(), "Lie algebra ideal of dimension $(dim(I)) over ", Lowercase())
    print(terse(io), base_lie_algebra(I))
  end
end

###############################################################################
#
#   Comparisons
#
###############################################################################

function Base.issubset(
  I1::LieAlgebraIdeal{C,LieT}, I2::LieAlgebraIdeal{C,LieT}
) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
  @req base_lie_algebra(I1) === base_lie_algebra(I2) "Incompatible Lie algebras."
  return all(in(I2), gens(I1))
end

function Base.:(==)(
  I1::LieAlgebraIdeal{C,LieT}, I2::LieAlgebraIdeal{C,LieT}
) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
  base_lie_algebra(I1) === base_lie_algebra(I2) || return false
  gens(I1) == gens(I2) && return true
  return issubset(I1, I2) && issubset(I2, I1)
end

function Base.hash(I::LieAlgebraIdeal, h::UInt)
  b = 0xd745d725d0b66431 % UInt
  h = hash(base_lie_algebra(I), h)
  h = hash(rref(basis_matrix(I)), h)
  return xor(h, b)
end

###############################################################################
#
#   Ideal membership
#
###############################################################################

@doc raw"""
    in(x::LieAlgebraElem, I::LieAlgebraIdeal) -> Bool

Return `true` if `x` is in the ideal `I`, `false` otherwise.
"""
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
) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
  @req base_lie_algebra(I1) === base_lie_algebra(I2) "Incompatible Lie algebras."
  return ideal(base_lie_algebra(I1), [gens(I1); gens(I2)])
end

@doc raw"""
    bracket(I1::LieAlgebraIdeal, I2::LieAlgebraIdeal) -> LieAlgebraIdeal

Return $[I_1,I_2]$.
"""
function bracket(
  I1::LieAlgebraIdeal{C,LieT}, I2::LieAlgebraIdeal{C,LieT}
) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
  @req base_lie_algebra(I1) === base_lie_algebra(I2) "Incompatible Lie algebras."
  return ideal(base_lie_algebra(I1), [x * y for x in gens(I1) for y in gens(I2)])
end

function bracket(
  L::LieAlgebra{C}, I::LieAlgebraIdeal{C,LieT}
) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
  @req L === base_lie_algebra(I) "Incompatible Lie algebras."
  return bracket(ideal(L), I)
end

###############################################################################
#
#   Important ideals and subalgebras
#
###############################################################################

@doc raw"""
    normalizer(L::LieAlgebra, I::LieAlgebraIdeal) -> LieSubalgebra

Return the normalizer of `I` in `L`, i.e. $\{x \in L \mid [x, I] \subseteq I\} = L$.
As `I` is an ideal in `L`, this is just `L`.
"""
function normalizer(L::LieAlgebra, I::LieAlgebraIdeal)
  @req base_lie_algebra(I) === L "Incompatible Lie algebras."
  return sub(L)
end

function normalizer(I::LieAlgebraIdeal)
  return normalizer(base_lie_algebra(I), I)
end

@doc raw"""
    centralizer(L::LieAlgebra, I::LieAlgebraIdeal) -> LieSubalgebra
  
Return the centralizer of `I` in `L`, i.e. $\{x \in L \mid [x, I] = 0\}$.
"""
function centralizer(L::LieAlgebra, I::LieAlgebraIdeal)
  return centralizer(L, basis(I))
end

function centralizer(I::LieAlgebraIdeal)
  return centralizer(base_lie_algebra(I), basis(I))
end

###############################################################################
#
#   Conversion
#
###############################################################################

@doc raw"""
    lie_algebra(I::LieAlgebraIdeal) -> LieAlgebra

Return `I` as a Lie algebra `LI`, together with an embedding `LI -> L`,
where `L` is the Lie algebra where `I` lives in.
"""
function lie_algebra(I::LieAlgebraIdeal)
  LI = lie_algebra(basis(I))
  L = base_lie_algebra(I)
  emb = hom(LI, L, basis(I))
  return LI, emb
end

@doc raw"""
    sub(L::LieAlgebra, I::LieAlgebraIdeal) -> LieSubalgebra

Return `I` as a subalgebra of `L`.
"""
function sub(L::LieAlgebra{C}, I::LieAlgebraIdeal{C}) where {C<:FieldElem}
  @req base_lie_algebra(I) === L "Incompatible Lie algebras."
  return sub(L, basis(I); is_basis=true)
end

###############################################################################
#
#   Constructor
#
###############################################################################

@doc raw"""
    ideal(L::LieAlgebra, gens::Vector{LieAlgebraElem}; is_basis::Bool=false) -> LieAlgebraIdeal

Return the smallest ideal of `L` containing `gens`.
If `is_basis` is `true`, then `gens` is assumed to be a basis of the ideal.
"""
function ideal(L::LieAlgebra, gens::Vector; is_basis::Bool=false)
  return LieAlgebraIdeal{elem_type(coefficient_ring(L)),elem_type(L)}(L, gens; is_basis)
end

@doc raw"""
    ideal(L::LieAlgebra, gen::LieAlgebraElem) -> LieAlgebraIdeal

Return the smallest ideal of `L` containing `gen`.
"""
function ideal(L::LieAlgebra{C}, gen::LieAlgebraElem{C}) where {C<:FieldElem}
  return ideal(L, [gen])
end

@doc raw"""
    ideal(L::LieAlgebra) -> LieAlgebraIdeal

Return `L` as an ideal of itself.
"""
function ideal(L::LieAlgebra)
  return ideal(L, basis(L); is_basis=true)
end
