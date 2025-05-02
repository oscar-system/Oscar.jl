# To be implemented by subtypes:
# Mandatory:
#   parent_type(::Type{MyLieAlgebraElem{C}}) = MyLieAlgebra{C}
#   elem_type(::Type{MyLieAlgebra{C}}) = MyLieAlgebraElem{C}
#   parent(x::MyLieAlgebraElem{C}) -> MyLieAlgebra{C}
#   coefficient_ring(L::MyLieAlgebra{C}) -> parent_type(C)
#   dim(L::MyLieAlgebra) -> Int
#   symbols(L::MyLieAlgebra) -> Vector{Symbol}
#   bracket(x::MyLieAlgebraElem{C}, y::MyLieAlgebraElem{C}) -> MyLieAlgebraElem{C}
#   Base.show(io::IO, x::MyLieAlgebra)
#   has_root_system(::MyLieAlgebra) -> Bool
#   root_system(::MyLieAlgebra) -> RootSystem
#   chevalley_basis(L::MyLieAlgebra{C}) -> NTuple{3,Vector{MyLieAlgebraElem{C}}}
#   set_root_system_and_chevalley_basis!(L::MyLieAlgebra{C}, R::RootSystem, chev::NTuple{3,Vector{MyLieAlgebraElem{C}}}}})

###############################################################################
#
#   Basic manipulation
#
###############################################################################

coefficient_ring(x::LieAlgebraElem) = coefficient_ring(parent(x))

number_of_generators(L::LieAlgebra) = dim(L)

gens(L::LieAlgebra) = basis(L)

gen(L::LieAlgebra, i::Int) = basis(L, i)

@doc raw"""
    dim(L::LieAlgebra) -> Int

Return the dimension of the Lie algebra `L`.
"""
dim(_::LieAlgebra) = error("Should be implemented by subtypes.")

@doc raw"""
    basis(L::LieAlgebra{C}) -> Vector{LieAlgebraElem{C}}

Return a basis of the Lie algebra `L`.
"""
basis(L::LieAlgebra) = [basis(L, i)::elem_type(L) for i in 1:dim(L)]

@doc raw"""
    basis(L::LieAlgebra{C}, i::Int) -> LieAlgebraElem{C}

Return the `i`-th basis element of the Lie algebra `L`.
"""
function basis(L::LieAlgebra, i::Int)
  @req 1 <= i <= dim(L) "Index out of bounds."
  R = coefficient_ring(L)
  v = zero_matrix(R, 1, dim(L))
  v[1, i] = one(R)
  return L(v)
end

@doc raw"""
    characteristic(L::LieAlgebra) -> Int

Return the characteristic of the coefficient ring of the Lie algebra `L`.
"""
characteristic(L::LieAlgebra) = characteristic(coefficient_ring(L))

@doc raw"""
    zero(L::LieAlgebra{C}) -> LieAlgebraElem{C}

Return the zero element of the Lie algebra `L`.
"""
function zero(L::LieAlgebra)
  mat = zero_matrix(coefficient_ring(L), 1, dim(L))
  return elem_type(L)(L, mat)
end

@doc raw"""
    iszero(x::LieAlgebraElem{C}) -> Bool

Check whether the Lie algebra element `x` is zero.
"""
function iszero(x::LieAlgebraElem)
  return iszero(_matrix(x))
end

@inline function _matrix(x::LieAlgebraElem{C}) where {C<:FieldElem}
  return (x.mat)::dense_matrix_type(C)
end

@doc raw"""
    coefficients(x::LieAlgebraElem{C}) -> Vector{C}

Return the coefficients of `x` with respect to [`basis(::LieAlgebra)`](@ref).
"""
function coefficients(x::LieAlgebraElem)
  return collect(_matrix(x))[1, :]
end

@doc raw"""
    coeff(x::LieAlgebraElem{C}, i::Int) -> C

Return the `i`-th coefficient of `x` with respect to [`basis(::LieAlgebra)`](@ref).
"""
function coeff(x::LieAlgebraElem, i::Int)
  return _matrix(x)[1, i]
end

@doc raw"""
    getindex(x::LieAlgebraElem{C}, i::Int) -> C

Return the `i`-th coefficient of `x` with respect to [`basis(::LieAlgebra)`](@ref).
"""
function getindex(x::LieAlgebraElem, i::Int)
  return coeff(x, i)
end

function Base.deepcopy_internal(x::LieAlgebraElem, dict::IdDict)
  return parent(x)(deepcopy_internal(_matrix(x), dict))
end

function check_parent(x1::LieAlgebraElem{C}, x2::LieAlgebraElem{C}) where {C<:FieldElem}
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
symbols(_::LieAlgebra) = error("Should be implemented by subtypes.")

function expressify(v::LieAlgebraElem, s=symbols(parent(v)); context=nothing)
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
function (L::LieAlgebra)()
  return zero(L)
end

@doc raw"""
    (L::LieAlgebra{C})(v::AbstractVector{Int}) -> LieAlgebraElem{C}

Return the element of `L` with coefficient vector `v`.
Fail, if `Int` cannot be coerced into the base ring of `L`.
"""
function (L::LieAlgebra)(v::AbstractVector{Int})
  return L(coefficient_ring(L).(v))
end

@doc raw"""
    (L::LieAlgebra{C})(v::AbstractVector{C}) -> LieAlgebraElem{C}

Return the element of `L` with coefficient vector `v`.
"""
function (L::LieAlgebra{C})(v::AbstractVector{C}) where {C<:FieldElem}
  @req length(v) == dim(L) "Length of vector does not match dimension."
  mat = matrix(coefficient_ring(L), 1, length(v), v)
  return elem_type(L)(L, mat)
end

@doc raw"""
    (L::LieAlgebra{C})(mat::MatElem{C}) -> LieAlgebraElem{C}

Return the element of `L` with coefficient vector equivalent to
the $1 \times \dim(L)$ matrix `mat`.
"""
function (L::LieAlgebra{C})(mat::MatElem{C}) where {C<:FieldElem}
  @req size(mat) == (1, dim(L)) "Invalid matrix dimensions."
  return elem_type(L)(L, mat)
end

@doc raw"""
    (L::LieAlgebra{C})(v::SRow{C}) -> LieAlgebraElem{C}

Return the element of `L` with coefficient vector `v`.
"""
function (L::LieAlgebra{C})(v::SRow{C}) where {C<:FieldElem}
  mat = dense_row(v, dim(L))
  return elem_type(L)(L, mat)
end

@doc raw"""
    (L::LieAlgebra{C})(x::LieAlgebraElem{C}) -> LieAlgebraElem{C}

Return `x`. Fails if `x` is not an element of `L`.
"""
function (L::LieAlgebra{C})(x::LieAlgebraElem{C}) where {C<:FieldElem}
  @req L === parent(x) "Incompatible modules."
  return x
end

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

function Base.:-(x::LieAlgebraElem{C}) where {C<:FieldElem}
  return parent(x)(-_matrix(x))
end

function Base.:+(x1::LieAlgebraElem{C}, x2::LieAlgebraElem{C}) where {C<:FieldElem}
  check_parent(x1, x2)
  return parent(x1)(_matrix(x1) + _matrix(x2))
end

function Base.:-(x1::LieAlgebraElem{C}, x2::LieAlgebraElem{C}) where {C<:FieldElem}
  check_parent(x1, x2)
  return parent(x1)(_matrix(x1) - _matrix(x2))
end

function Base.:*(x::LieAlgebraElem{C}, c::C) where {C<:FieldElem}
  coefficient_ring(x) != parent(c) && error("Incompatible rings.")
  return parent(x)(_matrix(x) * c)
end

function Base.:*(x::LieAlgebraElem, c::U) where {U<:RationalUnion}
  return parent(x)(_matrix(x) * c)
end

function Base.:*(c::C, x::LieAlgebraElem{C}) where {C<:FieldElem}
  coefficient_ring(x) != parent(c) && error("Incompatible rings.")
  return parent(x)(c * _matrix(x))
end

function Base.:*(c::U, x::LieAlgebraElem) where {U<:RationalUnion}
  return parent(x)(c * _matrix(x))
end

function Base.:*(x::LieAlgebraElem{C}, y::LieAlgebraElem{C}) where {C<:FieldElem}
  return bracket(x, y)
end

###############################################################################
#
#   Comparison functions
#
###############################################################################

function Base.:(==)(x1::LieAlgebraElem{C}, x2::LieAlgebraElem{C}) where {C<:FieldElem}
  check_parent(x1, x2)
  return _matrix(x1) == _matrix(x2)
end

function Base.hash(x::LieAlgebraElem, h::UInt)
  b = 0x6724cbedbd860982 % UInt
  h = hash(parent(x), h)
  h = hash(_matrix(x), h)
  return xor(h, b)
end

###############################################################################
#
#   Important ideals and subalgebras
#
###############################################################################

@doc raw"""
    derived_algebra(L::LieAlgebra) -> LieAlgebraIdeal

Return the derived algebra of `L`, i.e. $[L, L]$.
"""
function derived_algebra(L::LieAlgebra)
  L_ideal = ideal(L)
  return bracket(L_ideal, L_ideal)
end

@doc raw"""
    center(L::LieAlgebra) -> LieAlgebraIdeal

Return the center of `L`, i.e. $\{x \in L \mid [x, L] = 0\}$
"""
function center(L::LieAlgebra)
  dim(L) == 0 && return ideal(L, elem_type(L)[])

  mat = zero_matrix(coefficient_ring(L), dim(L), dim(L)^2)
  for (i, bi) in enumerate(basis(L))
    for (j, bj) in enumerate(basis(L))
      mat[i, ((j - 1) * dim(L) + 1):(j * dim(L))] = _matrix(bi * bj)
    end
  end

  ker = kernel(mat; side=:left)
  return ideal(L, L.(eachrow(ker)); is_basis=true)
end

@doc raw"""
    centralizer(L::LieAlgebra, xs::AbstractVector{<:LieAlgebraElem}) -> LieSubalgebra

Return the centralizer of `xs` in `L`, i.e. $\{y \in L \mid [x, y] = 0 \forall x \in xs\}$.
"""
function centralizer(L::LieAlgebra, xs::AbstractVector{<:LieAlgebraElem})
  @req all(x -> parent(x) === L, xs) "Incompatible Lie algebras."

  mat = zero_matrix(coefficient_ring(L), dim(L), dim(L) * length(xs))
  for (i, bi) in enumerate(basis(L))
    for (j, xj) in enumerate(xs)
      mat[i, ((j - 1) * dim(L) + 1):(j * dim(L))] = _matrix(bracket(bi, xj))
    end
  end

  c_basis = kernel(mat; side=:left)
  c_dim = nrows(c_basis)
  return sub(L, [L(c_basis[i, :]) for i in 1:c_dim]; is_basis=true)
end

@doc raw"""
    centralizer(L::LieAlgebra, x::LieAlgebraElem) -> LieSubalgebra

Return the centralizer of `x` in `L`, i.e. $\{y \in L \mid [x, y] = 0\}$.
"""
function centralizer(L::LieAlgebra, x::LieAlgebraElem)
  return centralizer(L, [x])
end

###############################################################################
#
#   Derived and central series
#
###############################################################################

@doc raw"""
    derived_series(L::LieAlgebra) -> Vector{LieAlgebraIdeal}

Return the derived series of `L`, i.e. the sequence of ideals 
$L^{(0)} = L$, $L^{(i + 1)} = [L^{(i)}, L^{(i)}]$.
"""
function derived_series(L::LieAlgebra)
  curr = ideal(L)
  series = [curr]
  while true
    next = bracket(curr, curr)
    dim(next) == dim(curr) && break
    push!(series, next)
    curr = next
  end
  return series
end

@doc raw"""
    lower_central_series(L::LieAlgebra) -> Vector{LieAlgebraIdeal}

Return the lower central series of `L`, i.e. the sequence of ideals
$L^{(0)} = L$, $L^{(i + 1)} = [L, L^{(i)}]$.
"""
function lower_central_series(L::LieAlgebra)
  curr = ideal(L)
  series = [curr]
  while true
    next = bracket(L, curr)
    dim(next) == dim(curr) && break
    push!(series, next)
    curr = next
  end
  return series
end

###############################################################################
#
#   Properties
#
###############################################################################

@doc raw"""
    is_abelian(L::LieAlgebra) -> Bool

Return `true` if `L` is abelian, i.e. $[L, L] = 0$.
"""
@attr Bool function is_abelian(L::LieAlgebra)
  b = basis(L)
  n = length(b)
  return all(iszero, b[i] * b[j] for i in 1:n for j in (i + 1):n)
end

@doc raw"""
    is_nilpotent(L::LieAlgebra) -> Bool

Return `true` if `L` is nilpotent, i.e. the lower central series of `L` terminates in $0$.
"""
@attr Bool function is_nilpotent(L::LieAlgebra)
  return dim(lower_central_series(L)[end]) == 0
end

@doc raw"""
    is_perfect(L::LieAlgebra) -> Bool

Return `true` if `L` is perfect, i.e. $[L, L] = L$.
"""
@attr Bool function is_perfect(L::LieAlgebra)
  return dim(derived_algebra(L)) == dim(L)
end

@doc raw"""
    is_semisimple(L::LieAlgebra) -> Bool

Return `true` if `L` is semisimple, i.e. the solvable radical of `L` is zero.

!!! warning
    This function is not implemented yet for all cases in positive characteristic.
"""
@attr Bool function is_semisimple(L::LieAlgebra)
  has_attribute(L, :is_simple) && is_simple(L) && return true
  is_invertible(killing_matrix(L)) && return true
  characteristic(L) == 0 && return false
  error("Not implemented.") # TODO
end

@doc raw"""
    is_simple(L::LieAlgebra) -> Bool

Return `true` if `L` is simple, i.e. `L` is not abelian and has no non-trivial ideals.

!!! warning
    This function is not implemented yet.
"""
@attr Bool function is_simple(L::LieAlgebra)
  is_abelian(L) && return false
  !is_perfect(L) && return false
  error("Not implemented.") # TODO
end

@doc raw"""
    is_solvable(L::LieAlgebra) -> Bool

Return `true` if `L` is solvable, i.e. the derived series of `L` terminates in $0$.
"""
@attr Bool function is_solvable(L::LieAlgebra)
  return dim(derived_series(L)[end]) == 0
end

###############################################################################
#
#   Adjoint elements
#
###############################################################################

@attr Vector{dense_matrix_type(C)} function adjoint_matrices(
  L::LieAlgebra{C}
) where {C<:FieldElem}
  return map(1:dim(L)) do i
    x = basis(L, i)
    A = zero_matrix(coefficient_ring(L), dim(L), dim(L))
    for (j, bj) in enumerate(basis(L))
      A[j, :] = _matrix(bracket(x, bj))
    end
    return A
  end
end

function adjoint_matrix(x::LieAlgebraElem{C}) where {C<:FieldElem}
  L = parent(x)
  mat = zero_matrix(coefficient_ring(L), dim(L), dim(L))
  tmp = zero(mat)
  for (c, g) in zip(_matrix(x), adjoint_matrices(L))
    mat = addmul!(mat, g, c, tmp)
  end
  return mat
end

function _adjoint_matrix(S::LieSubalgebra{C}, x::LieAlgebraElem{C}) where {C<:FieldElem}
  L = parent(x)
  @req base_lie_algebra(S) === L "Incompatible Lie algebras"
  A = zero_matrix(coefficient_ring(L), dim(S), dim(S))
  for (i, bi) in enumerate(basis(S))
    A[i, :] = coefficient_vector(bracket(x, bi), S)
  end
  return A
end

@attr dense_matrix_type(C) function killing_matrix(L::LieAlgebra{C}) where {C<:FieldElem}
  R = coefficient_ring(L)
  A = zero_matrix(R, dim(L), dim(L))
  for (i, adxi) in enumerate(adjoint_matrices(L))
    for (j, adxj) in enumerate(adjoint_matrices(L))
      i > j && continue # killing form is symmetric
      val = tr(adxi * adxj)
      A[j, i] = A[i, j] = val
    end
  end
  return A
end

@doc raw"""
    is_ad_nilpotent(x::LieAlgebraElem{C}) -> Bool

Return whether `x` is ad-nilpotent, i.e. whether the linear operator $\mathrm{ad}(x)$ is nilpotent.
"""
function is_ad_nilpotent(x::LieAlgebraElem{C}) where {C<:FieldElem}
  return is_nilpotent(adjoint_matrix(x))
end

###############################################################################
#
#   Root system detection
#
###############################################################################

@doc raw"""
    any_non_ad_nilpotent_element(L::LieAlgebra{C}) -> LieAlgebraElem{C}

Return an element of `L` that is not ad-nilpotent, or the zero element if all elements are ad-nilpotent.

The used algorithm is described in [GIR96; Ch. 3](@cite).
"""
function any_non_ad_nilpotent_element(L::LieAlgebra{C}) where {C<:FieldElem}
  if dim(L) <= 1
    # L is abelian and hence nilpotent
  elseif characteristic(L) == 0
    for x in basis(L)
      !is_ad_nilpotent(x) && return x
    end
    for (i, x) in enumerate(basis(L))
      for (j, y) in enumerate(basis(L))
        i > j && continue
        xy = x * y
        !is_ad_nilpotent(xy) && return xy
      end
    end
  else # characteristic > 0
    x = basis(L, 1)
    !is_ad_nilpotent(x) && return x
    K = sub(L, [x]; is_basis=true)
    while dim(K) < dim(L)
      # find an element b in L \ K with [b,K]⊆K
      N = normalizer(L, K)
      b = basis(N, findfirst(!in(K), basis(N)))
      !is_ad_nilpotent(b) && return b
      K = sub(L, [basis(K); b]; is_basis=true)
    end
  end
  set_attribute!(L, :is_nilpotent, true)
  return zero(L)
end

@doc raw"""
    engel_subalgebra(x::LieAlgebraElem{C}) -> LieSubalgebra{C,elem_type(parent(x))}

Return the Engel subalgebra of `x`, i.e. the generalized eigenspace of the linear operator $\mathrm{ad}(x)$.
"""
function engel_subalgebra(x::LieAlgebraElem{C}) where {C<:FieldElem}
  L = parent(x)
  n = dim(L)
  A = adjoint_matrix(x)^n
  ker = kernel(A; side=:left)
  basis = L.(eachrow(ker))
  L0adx = sub(L, basis; is_basis=true)
  return L0adx
end

function _cartan_subalgebra(L::LieAlgebra{C}) where {C<:FieldElem}
  F = coefficient_ring(L)
  n = dim(L)
  @req is_infinite(F) || length(F) > dim(L) "The implemented algorithm requires a large field"
  x = any_non_ad_nilpotent_element(L)
  if is_zero(x) # L is nilpotent
    return sub(L)
  end

  L0adx = engel_subalgebra(x)
  while true # decreasing variant is dim(L0adx)
    y = any_non_ad_nilpotent_element(L0adx)
    if is_zero(y) # L0adx is nilpotent
      return L0adx
    end

    c_itr =
      characteristic(F) == 0 ? (F(i) for i in 1:(n + 1)) : Iterators.filter(!iszero, F)
    z = x
    L0adz = L0adx
    for c in c_itr # at most n+1 iterations
      z = x + c * (y - x)
      L0adz = engel_subalgebra(z)
      if dim(L0adz) < dim(L0adx) && is_subset(L0adz, L0adx)
        break
      end
    end
    x = z
    L0adx = L0adz
  end
end

function _root_system_and_chevalley_basis(
  L::LieAlgebra{C}, H::LieSubalgebra{C}=_cartan_subalgebra(L)
) where {C<:FieldElem}
  @req base_lie_algebra(H) === L "Incompatible Lie algebras."
  # we just assume that H is indeed a Cartan subalgebra
  @req is_invertible(killing_matrix(L)) "The Killing form is degenerate"

  R = coefficient_ring(L)

  # compute the common eigenspaces of the adjoint action of H.
  # B is a list of subspaces of L that gets refined in each iteration.
  # to exploit existing functionality, we use the LieSubalgebra type even
  # though the subspaces are in general not Lie subalgebras. With setting
  # is_basis=true, we ensure that the subspace generators do not get
  # extended to a subalgebra.
  B = [sub(L, basis(L); is_basis=true)]
  for h in basis(H)
    B_new = empty(B)
    for B_j in B
      A = _adjoint_matrix(B_j, h)
      facs = factor(minimal_polynomial(A))
      for (f, k) in facs
        @assert k == 1 # TODO: is this always the case?
        ker = kernel((f^k)(A); side=:left)
        basis = B_j.(eachrow(ker))
        push!(B_new, sub(L, basis; is_basis=true))
      end
    end
    B = B_new
  end
  filter!(!=(H), B)
  @req all(B_j -> dim(B_j) == 1, B) "The Cartan subalgebra is not split"

  # compute the roots, i.e. the list of eigenvalues of basis(H) on each B_j
  root_spaces = Dict(
    begin
      b = only(basis(B_j))
      root = [only(solve(_matrix(b), _matrix(bracket(h, b)); side=:left)) for h in basis(H)]
      root => B_j
    end for B_j in B
  )
  roots = collect(keys(root_spaces))

  # compute an R-basis of the root space, s.t. the corresponding co-roots are a basis of H
  roots_basis = empty(roots)
  basis_mat_H = zero_matrix(R, 0, dim(L))
  for root in roots
    nrows(basis_mat_H) == dim(H) && break
    x_j = only(basis(root_spaces[root]))
    y_j = only(basis(root_spaces[-root]))
    h_j = bracket(x_j, y_j)
    if !can_solve(basis_mat_H, _matrix(h_j); side=:left)
      basis_mat_H = vcat(basis_mat_H, _matrix(h_j))
      push!(roots_basis, root)
    end
  end

  function CartanInt(roots::Vector{RootType}, a::RootType, b::RootType) where {RootType}
    # `a` and `b` are two roots in `roots`.
    a == b && return 2
    a == -b && return -2
    # If a != ±b, the Cartan integer of `a` and `b` is `s-t`, where
    # `s` and `t` are the largest integers such that `b-s*a` and `b+t*a` are still roots.
    rt = b - a
    s = 0
    while rt in roots
      s += 1
      rt -= a
    end
    rt = b + a
    t = 0
    while rt in roots
      t += 1
      rt += a
    end
    return s - t
  end

  # we define a root to be positive if the first of its non-zero Cartan integers with the R-basis is positive
  roots_positive = empty(roots)
  for root in roots
    2 * length(roots_positive) == length(roots) && break
    root in roots_positive && continue
    -root in roots_positive && continue
    c = first(
      Iterators.dropwhile(
        iszero, Iterators.map(b_j -> CartanInt(roots, root, b_j), roots_basis)
      ),
    )
    if c > 0
      push!(roots_positive, root)
    else
      push!(roots_positive, -root)
    end
  end

  # a positive root is simple if it is not the sum of two positive roots
  roots_simple = empty(roots)
  roots_positive_sums = Set(
    alpha_i + alpha_j for (i, alpha_i) in enumerate(roots_positive) for
    (j, alpha_j) in enumerate(roots_positive) if i < j
  )
  for root in roots_positive
    if !(root in roots_positive_sums)
      push!(roots_simple, root)
    end
  end
  @assert length(roots_simple) == dim(H)

  # compute the Cartan matrix and abstract root system
  cm = matrix(
    ZZ,
    [
      CartanInt(roots, alpha_i, alpha_j) for alpha_i in roots_simple,
      alpha_j in roots_simple
    ],
  )
  type, ordering = cartan_type_with_ordering(cm; check=false)
  permute!(roots_simple, ordering)
  rs = root_system(type)

  # compute conrete root vectors of L
  root_vectors = [
    begin
      concrete_root = sum(
        Int(c) .* root_simple for
        (c, root_simple) in zip(coefficients(abstract_root), roots_simple)
      )
      @assert concrete_root in roots
      root_vector = only(basis(root_spaces[concrete_root]))
    end for abstract_root in Oscar.roots(rs)
  ]

  # scale the root vectors to get canonical generators, i.e.
  # 1. [h_i,h_j] = 0            2. [x_i,y_j] = delta_ij * h_i
  # 3. [h_i,x_j] = a_ij * x_j   4. [h_i,y_j] = -a_ij * y_j
  canonical_gens_xs = root_vectors[1:n_simple_roots(rs)]
  canonical_gens_ys = Vector{elem_type(L)}([
    begin
      x = canonical_gens_xs[i]
      y = root_vectors[n_positive_roots(rs) + i]
      y *= 2//only(solve(_matrix(x), _matrix(bracket(bracket(x, y), x)); side=:left))
      y
    end for i in 1:n_simple_roots(rs)
  ])
  canonical_gens_hs = elem_type(L)[
    canonical_gens_xs[i] * canonical_gens_ys[i] for i in 1:n_simple_roots(rs)
  ]

  # Construct an automorphism `f` of `L` with `f(L_alpha) = L_-alpha`, and `f(h) = -h` for `h` in `H`.
  # This is already defined by mapping `canonical_gens_xs` to `canonical_gens_ys` and vice versa.
  # To save computation time, we only compute f on L_alpha for alpha a positive root
  domain_basis = copy(canonical_gens_xs)
  codomain_basis = copy(canonical_gens_ys)
  for (i, alpha) in enumerate(positive_roots(rs))
    i <= n_simple_roots(rs) && continue # already set above
    j, k = 0, 0
    found = false
    for outer j in 1:rank(rs)
      (fl, k) = is_positive_root_with_index(alpha - simple_root(rs, j))
      if fl
        found = true
        break
      end
    end
    @assert found
    push!(domain_basis, canonical_gens_xs[j] * domain_basis[k])
    push!(codomain_basis, canonical_gens_ys[j] * codomain_basis[k])
  end

  # For each positive root vector `x` we set `y = -f(x)`.
  # We compute a scalar `coeff` such that `h = 2/coeff * [x,y]`.
  # We need to multiply `x` and `y` by `sqrt(2/coeff)` to get elements of a Chevally basis.
  pos_root_vectors = elem_type(L)[]
  neg_root_vectors = elem_type(L)[]
  scaled_cartan_elems = elem_type(L)[]
  coeffs = C[]
  for i in 1:n_positive_roots(rs)
    x = root_vectors[i]
    y = -only(solve(_matrix(domain_basis[i]), _matrix(x); side=:left)) * codomain_basis[i] # y = -f(x)
    h = bracket(x, y)
    coeff = only(solve(_matrix(x), _matrix(bracket(h, x)); side=:left))
    if i <= n_simple_roots(rs)
      push!(scaled_cartan_elems, (2//coeff) * h)
    end
    push!(coeffs, 2//coeff)
    push!(pos_root_vectors, x)
    push!(neg_root_vectors, y)
  end

  # In general, `coeffs` may not be squares in the coefficient field.
  # Construct a field extension `F` where they are.
  if all(is_square, coeffs)
    F = R
    coeffs_F_sqrt = sqrt.(coeffs)
  else
    F, coeffs_F_sqrt = _field_ext_with_sqrts(R, coeffs)
  end

  # Construct a Lie algebra `L_F` over `F` with the same structure constants as `L`,
  # and construct the Chevalley basis of `L_F`. 
  L_F = change_base_ring(F, L)
  L_F_chev_basis_mat = reduce(vcat,
    [
      (coeff .* _matrix(rv) for (coeff, rv) in zip(coeffs_F_sqrt, pos_root_vectors))...,
      (coeff .* _matrix(rv) for (coeff, rv) in zip(coeffs_F_sqrt, neg_root_vectors))...,
      (F.(_matrix(h)) for h in scaled_cartan_elems)...,
    ]; init=zero_matrix(F, 0, dim(L_F)))
  # The structure constants of `L_F` w.r.t. the Chevalley basis lie in `R`,
  # so we can construct the Lie algebra `K` over `R` with the these structure constants.
  # The basis elements of `K` form a Chevalley basis.
  K = lie_algebra(
    R, map(e -> change_base_ring(R, e), _structure_constant_table(L_F, L_F_chev_basis_mat))
  )

  # Construct an isomorphism f from K to L.
  # To save computation time, we only compute f on L_alpha for alpha a root.
  domain_basis_xs = [basis(K, i) for i in 1:n_simple_roots(rs)]
  domain_basis_ys = [basis(K, n_positive_roots(rs) + i) for i in 1:n_simple_roots(rs)]
  codomain_basis_xs = copy(canonical_gens_xs)
  codomain_basis_ys = copy(canonical_gens_ys)
  for (i, alpha) in enumerate(positive_roots(rs))
    i <= n_simple_roots(rs) && continue # already set above
    j, k = 0, 0
    found = false
    for outer j in 1:rank(rs)
      (fl, k) = is_positive_root_with_index(alpha - simple_root(rs, j))
      if fl
        found = true
        break
      end
    end
    @assert found
    push!(domain_basis_xs, domain_basis_xs[j] * domain_basis_xs[k])
    push!(domain_basis_ys, domain_basis_ys[j] * domain_basis_ys[k])
    push!(codomain_basis_xs, codomain_basis_xs[j] * codomain_basis_xs[k])
    push!(codomain_basis_ys, codomain_basis_ys[j] * codomain_basis_ys[k])
  end

  # The image of the Chevalley basis of K under f is a Chevalley basis of L.
  xs = [
    only(solve(_matrix(domain_basis_xs[i]), _matrix(basis(K, i)); side=:left)) *
    codomain_basis_xs[i] for i in 1:n_positive_roots(rs)
  ]
  ys = [
    only(
      solve(
        _matrix(domain_basis_ys[i]), _matrix(basis(K, n_positive_roots(rs) + i)); side=:left
      ),
    ) * codomain_basis_ys[i] for i in 1:n_positive_roots(rs)
  ]

  return rs, (xs, ys, scaled_cartan_elems)
end

# TODO: clean the following up, once we have a unified interface for field extensions
function _field_ext_with_sqrts(F::QQField, elems::Vector{QQFieldElem})
  i = findfirst(!is_square, elems)
  if isnothing(i)
    return F, sqrt.(elems)
  end
  Fx, x = polynomial_ring(F; cached=false)
  F, _ = number_field(x^2 - elems[i]; cached=false)
  return _field_ext_with_sqrts(F, [F(c) for c in elems])
end

function _field_ext_with_sqrts(F::NumField, elems::Vector{<:NumFieldElem})
  @req all(e -> parent(e) === F, elems) "Incompatible parent fields"
  elems_F = elems
  i = findfirst(!is_square, elems_F)
  while !isnothing(i)
    Fx, x = polynomial_ring(F; cached=false)
    F, _ = number_field(x^2 - elems_F[i]; cached=false)
    elems_F = [F(c) for c in elems_F]
    i = findfirst(!is_square, elems_F)
  end
  return F, sqrt.(elems_F)
end

function _field_ext_with_sqrts(F::FqField, elems::Vector{<:FqFieldElem})
  @req all(e -> parent(e) === F, elems) "Incompatible parent fields"
  elems_F = elems
  i = findfirst(!is_square, elems_F)
  while !isnothing(i)
    Fx, x = polynomial_ring(F; cached=false)
    F, _ = finite_field(x^2 - elems_F[i]; cached=false)
    elems_F = [F(c) for c in elems_F]
    i = findfirst(!is_square, elems_F)
  end
  return F, sqrt.(elems_F)
end

function assure_root_system(L::LieAlgebra{C}) where {C<:FieldElem}
  if !has_root_system(L)
    R, chev = _root_system_and_chevalley_basis(L)
    set_root_system_and_chevalley_basis!(L, R, chev)
  end
  @assert has_root_system(L)
end

###############################################################################
#
#   Root system getters
#
###############################################################################

@doc raw"""
    has_root_system(L::LieAlgebra) -> Bool

Return whether a root system for `L` is known.
"""
has_root_system(L::LieAlgebra) = false # to be implemented by subtypes

@doc raw"""
    root_system(L::LieAlgebra) -> RootSystem

Return the root system of `L`.

This function will error if no root system is known and none can be computed.
"""
function root_system(L::LieAlgebra) # to be implemented by subtypes
  throw(Hecke.NotImplemented())
end

@doc raw"""
    chevalley_basis(L::LieAlgebra) -> NTuple{3,Vector{elem_type(L)}}

Return the Chevalley basis of the Lie algebra `L` in three vectors, stating first the positive root vectors, 
then the negative root vectors, and finally the basis of the Cartan subalgebra. The order of root vectors corresponds
to the order of the roots in [`root_system(::LieAlgebra)`](@ref).

This function will error if no root system is known and none can be computed.
"""
function chevalley_basis(L::LieAlgebra) # to be implemented by subtypes
  throw(Hecke.NotImplemented())
end

@doc raw"""
    cartan_matrix(L::LieAlgebra) -> ZZMatrix

Return the Cartan matrix of the root system of `L`.
"""
function cartan_matrix(L::LieAlgebra{C}) where {C<:FieldElem}
  return cartan_matrix(root_system(L))
end

function cartan_matrix_inv(L::LieAlgebra{C}) where {C<:FieldElem}
  return cartan_matrix_inv(root_system(L))
end

@doc raw"""
    cartan_subalgebra(L::LieAlgebra{C}) where {C<:FieldElem} -> LieSubalgebra{C,elem_type(L)}

Return a Cartan subalgebra of `L`.

If `L` knows its root system, this function uses the Chevalley basis to construct a Cartan subalgebra.
Otherwise, it uses the algorithm described in [Gra00; Ch. 3.2](@cite).
The return value of this function may change when the root system of `L` is first computed.
"""
function cartan_subalgebra(L::LieAlgebra{C}) where {C<:FieldElem}
  if has_root_system(L)
    return sub(L, chevalley_basis(L)[3]; is_basis=true)
  else
    return _cartan_subalgebra(L)
  end
end

###############################################################################
#
#   Structure constant table
#
###############################################################################

function structure_constant_table(L::LieAlgebra; copy::Bool=true)
  # copy gets ignored as the generic implementation creates a new table anyway

  R = coefficient_ring(L)
  struct_consts = Matrix{sparse_row_type(R)}(undef, dim(L), dim(L))

  basis = Oscar.basis(L)

  for i in 1:dim(L)
    struct_consts[i, i] = sparse_row(R)
  end
  for i in 1:dim(L), j in (i + 1):dim(L)
    struct_consts[i, j] = sparse_row(_matrix(bracket(basis[i], basis[j])))
    struct_consts[j, i] = -struct_consts[i, j]
  end

  return struct_consts
end

function structure_constant_table(
  L::LieAlgebra{C}, basis::Vector{<:LieAlgebraElem{C}}
) where {C}
  @req all(b -> parent(b) === L, basis) "Incompatible Lie algebras"
  basis_mat = reduce(
    vcat, _matrix.(basis); init=zero_matrix(coefficient_ring(L), 0, dim(L))
  )
  return _structure_constant_table(L, basis_mat)
end

function _structure_constant_table(L::LieAlgebra{C}, basis_mat::MatElem{C}) where {C}
  @req is_invertible(basis_mat) "input is not a basis"
  basis_mat_ctx = solve_init(basis_mat)

  R = coefficient_ring(L)
  struct_consts = Matrix{sparse_row_type(R)}(undef, dim(L), dim(L))
  sc_standard_basis = structure_constant_table(L; copy=false)

  for i in 1:dim(L)
    struct_consts[i, i] = sparse_row(R)
  end
  for i in 1:dim(L), j in (i + 1):dim(L)
    # lookup entry in `basis_mat * sc_standard_basis * transpose(basis_mat)`
    entry = zero_matrix(R, 1, dim(L))
    for k in 1:dim(L), l in 1:dim(L)
      coeff = basis_mat[i, k] * basis_mat[j, l]
      if !is_zero(coeff)
        entry = addmul!(entry, coeff, dense_row(sc_standard_basis[k, l], dim(L)))
      end
    end
    struct_consts[i, j] = sparse_row(solve(basis_mat_ctx, entry; side=:left))
    struct_consts[j, i] = -struct_consts[i, j]
  end

  return struct_consts
end

###############################################################################
#
#   Universal enveloping algebra
#
###############################################################################

@doc raw"""
    universal_enveloping_algebra(L::LieAlgebra; ordering::Symbol=:lex) -> PBWAlgRing, Map

Return the universal enveloping algebra `U(L)` of `L` with the given monomial ordering,
together with a map from `L` into the filtered component of degree 1 of `U(L)`.
"""
function universal_enveloping_algebra(L::LieAlgebra; ordering::Symbol=:lex)
  R, gensR = polynomial_ring(coefficient_ring(L), symbols(L))
  n = dim(L)
  b = basis(L)

  to_R(x::LieAlgebraElem) =
    sum(c * g for (c, g) in zip(_matrix(x), gensR); init=zero(R))

  rel = strictly_upper_triangular_matrix([
    to_R(b[i]) * to_R(b[j]) - to_R(b[i] * b[j]) for i in 1:(n - 1) for j in (i + 1):n
  ])
  U, gensU = pbw_algebra(R, rel, monomial_ordering(R, ordering); check=true)

  L_to_U = MapFromFunc(
    L, U, function (x::LieAlgebraElem)
      sum(c * g for (c, g) in zip(_matrix(x), gensU); init=zero(U))
    end
  )
  return U, L_to_U
end

###############################################################################
#
#   Constructor
#
###############################################################################

@doc raw"""
    abelian_lie_algebra(R::Field, n::Int) -> LinearLieAlgebra{elem_type(R)}
    abelian_lie_algebra(::Type{LinearLieAlgebra}, R::Field, n::Int) -> LinearLieAlgebra{elem_type(R)}
    abelian_lie_algebra(::Type{AbstractLieAlgebra}, R::Field, n::Int) -> AbstractLieAlgebra{elem_type(R)}

Return the abelian Lie algebra of dimension `n` over the field `R`.
The first argument can be optionally provided to specify the type of the returned
Lie algebra.

# Examples
```jldoctest
julia> abelian_lie_algebra(LinearLieAlgebra, QQ, 3)
Linear Lie algebra with 3x3 matrices
  of dimension 3
over rational field

julia> abelian_lie_algebra(AbstractLieAlgebra, QQ, 3)
Abstract Lie algebra
  of dimension 3
over rational field
```
"""
function abelian_lie_algebra(R::Field, n::Int)
  @req n >= 0 "Dimension must be non-negative."
  return abelian_lie_algebra(LinearLieAlgebra, R, n)
end

@doc raw"""
    lie_algebra(gapL::GapObj, s::Vector{<:VarName}) -> LieAlgebra{elem_type(R)}

Construct a Lie algebra isomorphic to the GAP Lie algebra `gapL`. Its basis element are named by `s`,
or by `x_i` by default.
We require `gapL` to be a finite-dimensional GAP Lie algebra. The return type is dependent on
properties of `gapL`, in particular, whether GAP knows about a matrix representation.

"""
function lie_algebra(
  gapL::GapObj,
  s::Vector{<:VarName}=[Symbol("x_$i") for i in 1:GAPWrap.Dimension(gapL)];
)
  @req GAPWrap.IsLieAlgebra(gapL) "input must be a Lie algebra."
  if GAPWrap.IsFiniteDimensional(gapL)
    if GAPWrap.IsLieObjectCollection(gapL)
      return codomain(_iso_gap_oscar_linear_lie_algebra(gapL, s))
    else
      return codomain(_iso_gap_oscar_abstract_lie_algebra(gapL, s))
    end
  end
  error("Not implemented.")
end
