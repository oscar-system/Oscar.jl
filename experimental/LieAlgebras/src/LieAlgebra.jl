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
  return iszero(coefficients(x))
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
  return coefficients(x1) == coefficients(x2)
end

function Base.hash(x::LieAlgebraElem, h::UInt)
  b = 0x6724cbedbd860982 % UInt
  h = hash(parent(x), h)
  h = hash(coefficients(x), h)
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
  for (c, g) in zip(coefficients(x), adjoint_matrices(L))
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

  F = coefficient_ring(L)

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
  basis_mat_H = zero_matrix(F, 0, dim(L))
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
  R = root_system(type)

  # compute a Chevalley basis of L
  root_vectors = [
    begin
      concrete_root = sum(
        c .* root_simple for
        (c, root_simple) in zip(coefficients(abstract_root), roots_simple)
      )
      @assert concrete_root in roots
      root_vector = only(basis(root_spaces[concrete_root]))
    end for abstract_root in Oscar.roots(R)
  ]
  xs = root_vectors[1:n_positive_roots(R)]
  ys = Vector{elem_type(L)}([
    begin
      x = xs[i]
      y = root_vectors[n_positive_roots(R) + i]
      y *= 2//only(solve(_matrix(x), _matrix(bracket(bracket(x, y), x)); side=:left))
      y
    end for i in 1:n_positive_roots(R)
  ])
  hs = elem_type(L)[xs[i] * ys[i] for i in 1:n_simple_roots(R)]

  return R, (xs, ys, hs)
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
    sum(c * g for (c, g) in zip(coefficients(x), gensR); init=zero(R))

  rel = strictly_upper_triangular_matrix([
    to_R(b[i]) * to_R(b[j]) - to_R(b[i] * b[j]) for i in 1:(n - 1) for j in (i + 1):n
  ])
  U, gensU = pbw_algebra(R, rel, monomial_ordering(R, ordering); check=true)

  L_to_U = MapFromFunc(
    L, U, function (x::LieAlgebraElem)
      sum(c * g for (c, g) in zip(coefficients(x), gensU); init=zero(U))
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
