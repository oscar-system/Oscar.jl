###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_lie_algebra(S::LieSubalgebra{C,LieT}) where {C<:FieldElem,LieT<:LieAlgebraElem{C}} =
  S.base_lie_algebra::parent_type(LieT)

function gens(S::LieSubalgebra)
  return S.gens
end

function gen(S::LieSubalgebra, i::Int)
  return S.gens[i]
end

function number_of_generators(S::LieSubalgebra)
  return length(gens(S))
end

function basis_matrix(S::LieSubalgebra{C,LieT}) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
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

coefficient_ring(S::LieSubalgebra) = coefficient_ring(base_lie_algebra(S))

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, mime::MIME"text/plain", S::LieSubalgebra)
  @show_name(io, S)
  @show_special(io, mime, S)
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
  @show_name(io, S)
  @show_special(io, S)
  io = pretty(io)
  if is_terse(io)
    print(io, LowercaseOff(), "Lie subalgebra")
  else
    print(io, LowercaseOff(), "Lie subalgebra of dimension $(dim(S)) of ", Lowercase())
    print(terse(io), base_lie_algebra(S))
  end
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

@doc raw"""
    (S::LieSubalgebra{C})() -> LieAlgebraElem{C}

Return the zero element of the Lie subalgebra `S`.
"""
function (S::LieSubalgebra)()
  return zero(base_lie_algebra(S))
end

@doc raw"""
    (S::LieSubalgebra{C})(v::AbstractVector{Int}) -> LieAlgebraElem{C}

Return the element of `S` with coefficient vector `v`.
Fail, if `Int` cannot be coerced into the base ring of `S`.
"""
function (S::LieSubalgebra)(v::AbstractVector{Int})
  return S(coefficient_ring(S).(v))
end

@doc raw"""
    (S::LieSubalgebra{C})(v::AbstractVector{C}) -> LieAlgebraElem{C}

Return the element of `S` with coefficient vector `v`.
"""
function (S::LieSubalgebra{C})(v::AbstractVector{C}) where {C<:FieldElem}
  @req length(v) == dim(S) "Length of vector does not match dimension."
  mat = matrix(coefficient_ring(S), 1, length(v), v)
  L = base_lie_algebra(S)
  return elem_type(L)(L, mat * basis_matrix(S))
end

@doc raw"""
    (S::LieSubalgebra{C})(mat::MatElem{C}) -> LieAlgebraElem{C}

Return the element of `S` with coefficient vector equivalent to
the $1 \times \dim(S)$ matrix `mat`.
"""
function (S::LieSubalgebra{C})(mat::MatElem{C}) where {C<:FieldElem}
  @req size(mat) == (1, dim(S)) "Invalid matrix dimensions."
  L = base_lie_algebra(S)
  return elem_type(L)(L, mat * basis_matrix(S))
end

###############################################################################
#
#   Comparisons
#
###############################################################################

function Base.issubset(
  S1::LieSubalgebra{C,LieT}, S2::LieSubalgebra{C,LieT}
) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
  @req base_lie_algebra(S1) === base_lie_algebra(S2) "Incompatible Lie algebras."
  return all(in(S2), gens(S1))
end

function Base.:(==)(
  S1::LieSubalgebra{C,LieT}, S2::LieSubalgebra{C,LieT}
) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
  base_lie_algebra(S1) === base_lie_algebra(S2) || return false
  gens(S1) == gens(S2) && return true
  return issubset(S1, S2) && issubset(S2, S1)
end

function Base.hash(S::LieSubalgebra, h::UInt)
  b = 0x712b12e671184d1a % UInt
  h = hash(base_lie_algebra(S), h)
  h = hash(rref(basis_matrix(S)), h)
  return xor(h, b)
end

###############################################################################
#
#   Subalgebra membership
#
###############################################################################

@doc raw"""
    in(x::LieAlgebraElem{C}, S::LieSubalgebra{C}) -> Bool

Return `true` if `x` is in the Lie subalgebra `S`, `false` otherwise.
"""
function Base.in(x::LieAlgebraElem{C}, S::LieSubalgebra{C}) where {C<:FieldElem}
  @req parent(x) === base_lie_algebra(S) "Incompatible Lie algebras"
  return can_solve(basis_matrix(S), _matrix(x); side=:left)
end

@doc raw"""
    Oscar.LieAlgebras.coefficient_vector(x::LieAlgebraElem{C}, S::LieSubalgebra{C}) -> Vector{C}

Return the coefficient vector of `x` in the basis of `S`.
This function will throw an error if `x` is not in `S`.
"""
function coefficient_vector(x::LieAlgebraElem{C}, S::LieSubalgebra{C}) where {C<:FieldElem}
  @req parent(x) === base_lie_algebra(S) "Incompatible Lie algebras"
  return solve(basis_matrix(S), _matrix(x); side=:left)
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
) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
  @req base_lie_algebra(S1) === base_lie_algebra(S2) "Incompatible Lie algebras."
  return ideal(base_lie_algebra(S1), [x * y for x in gens(S1) for y in gens(S2)])
end

function bracket(
  L::LieAlgebra{C}, S::LieSubalgebra{C,LieT}
) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
  @req L === base_lie_algebra(S) "Incompatible Lie algebras."
  return bracket(ideal(L), S)
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
  @req base_lie_algebra(S) === L "Incompatible Lie algebras."

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
  sol = kernel(mat; side=:left)
  sol_dim = nrows(sol)
  sol = sol[:, 1:dim(L)]
  c_dim, c_basis = rref(sol)
  return sub(L, [L(c_basis[i, :]) for i in 1:c_dim]; is_basis=true)
end

function normalizer(S::LieSubalgebra)
  return normalizer(base_lie_algebra(S), S)
end

@doc raw"""
    centralizer(L::LieAlgebra, S::LieSubalgebra) -> LieSubalgebra
  
Return the centralizer of `S` in `L`, i.e. $\{x \in L \mid [x, S] = 0\}$.
"""
function centralizer(L::LieAlgebra, S::LieSubalgebra)
  return centralizer(L, basis(S))
end

function centralizer(S::LieSubalgebra)
  return centralizer(base_lie_algebra(S), basis(S))
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
#   More misc stuff
#
###############################################################################

function any_non_ad_nilpotent_element(S::LieSubalgebra)
  LS, emb = lie_algebra(S)
  x = any_non_ad_nilpotent_element(LS)
  return emb(x)
end

###############################################################################
#
#   Conversion
#
###############################################################################

@doc raw"""
    lie_algebra(S::LieSubalgebra) -> LieAlgebra, LieAlgebraHom

Return `S` as a Lie algebra `LS`, together with an embedding `LS -> L`,
where `L` is the Lie algebra where `S` lives in.
"""
function lie_algebra(S::LieSubalgebra)
  L = base_lie_algebra(S)
  LS = lie_algebra(L, basis(S))
  emb = hom(LS, L, basis(S); check=false)
  return LS, emb
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
function sub(L::LieAlgebra{C}, gen::LieAlgebraElem{C}) where {C<:FieldElem}
  return sub(L, [gen])
end

@doc raw"""
    sub(L::LieAlgebra) -> LieSubalgebra

Return `L` as a Lie subalgebra of itself.
"""
function sub(L::LieAlgebra)
  return sub(L, basis(L); is_basis=true)
end
