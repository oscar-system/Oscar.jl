
###############################################################################
#
#  Accessors
#
###############################################################################

@doc raw"""
    space(Vf::QuadSpaceWithIsom) -> QuadSpace

Given a quadratic space with isometry $(V, f)$, return the underlying space $V$.

# Examples
```jldoctest
julia> V = quadratic_space(QQ, 2);

julia> Vf = quadratic_space_with_isometry(V; neg=true);

julia> space(Vf) === V
true
```
"""
space(Vf::QuadSpaceWithIsom) = Vf.V

@doc raw"""
    isometry(Vf::QuadSpaceWithIsom) -> QQMatrix

Given a quadratic space with isometry $(V, f)$, return the underlying isometry
$f$.

# Examples
```jldoctest
julia> V = quadratic_space(QQ, 2);

julia> Vf = quadratic_space_with_isometry(V; neg=true);

julia> isometry(Vf)
[-1    0]
[ 0   -1]
```
"""
isometry(Vf::QuadSpaceWithIsom) = Vf.f

@doc raw"""
    order_of_isometry(Vf::QuadSpaceWithIsom) -> IntExt

Given a quadratic space with isometry $(V, f)$, return the order of the
underlying isometry $f$.

# Examples
```jldoctest
julia> V = quadratic_space(QQ, 2);

julia> Vf = quadratic_space_with_isometry(V; neg=true);

julia> order_of_isometry(Vf) == 2
true
```
"""
order_of_isometry(Vf::QuadSpaceWithIsom) = Vf.n

###############################################################################
#
#  Attributes
#
###############################################################################

@doc raw"""
    rank(Vf::QuadSpaceWithIsom) -> Integer

Given a quadratic space with isometry $(V, f)$, return the rank of the
underlying space $V$.

See [`rank(::AbstractSpace)`](@ref).

# Examples
```jldoctest
julia> V = quadratic_space(QQ, 2);

julia> Vf = quadratic_space_with_isometry(V);

julia> rank(Vf) == 2
true
```
"""
rank(Vf::QuadSpaceWithIsom) = rank(space(Vf))

@doc raw"""
    dim(Vf::QuadSpaceWithIsom) -> Integer

Given a quadratic space with isometry $(V, f)$, return the dimension of the
underlying space of $V$.

See [`dim(::AbstractSpace)`](@ref).

# Examples
```jldoctest
julia> V = quadratic_space(QQ, 2);

julia> Vf = quadratic_space_with_isometry(V);

julia> dim(Vf) == 2
true
```
"""
dim(Vf::QuadSpaceWithIsom) = dim(space(Vf))

@doc raw"""
    characteristic_polynomial(Vf::QuadSpaceWithIsom) -> QQPolyRingElem

Given a quadratic space with isometry $(V, f)$, return the characteristic
polynomial of the underlying isometry $f$.

# Examples
```jldoctest
julia> V = quadratic_space(QQ, 2);

julia> Vf = quadratic_space_with_isometry(V; neg=true);

julia> characteristic_polynomial(Vf)
x^2 + 2*x + 1
```
"""
characteristic_polynomial(Vf::QuadSpaceWithIsom) = characteristic_polynomial(isometry(Vf))

@doc raw"""
    minimal_polynomial(Vf::QuadSpaceWithIsom) -> QQPolyRingElem

Given a quadratic space with isometry $(V, f)$, return the minimal
polynomial of the underlying isometry $f$.

# Examples
```jldoctest
julia> V = quadratic_space(QQ, 2);

julia> Vf = quadratic_space_with_isometry(V; neg=true);

julia> minimal_polynomial(Vf)
x + 1
```
"""
minimal_polynomial(Vf::QuadSpaceWithIsom) = minimal_polynomial(isometry(Vf))

@doc raw"""
    gram_matrix(Vf::QuadSpaceWithIsom) -> QQMatrix

Given a quadratic space with isometry $(V, f)$, return the Gram matrix
of the underlying space $V$ with respect to its standard basis.

# Examples
```jldoctest
julia> V = quadratic_space(QQ, 2);

julia> Vf = quadratic_space_with_isometry(V; neg=true);

julia> is_one(gram_matrix(Vf))
true
```
"""
gram_matrix(Vf::QuadSpaceWithIsom) = gram_matrix(space(Vf))

@doc raw"""
    det(Vf::QuadSpaceWithIsom) -> QQFieldElem

Given a quadratic space with isometry $(V, f)$, return the determinant
of the underlying space $V$.

See [`det(::AbstractSpace)`](@ref).

# Examples
```jldoctest
julia> V = quadratic_space(QQ, 2);

julia> Vf = quadratic_space_with_isometry(V; neg=true);

julia> is_one(det(Vf))
true
```
"""
det(Vf::QuadSpaceWithIsom) = det(space(Vf))

@doc raw"""
    discriminant(Vf::QuadSpaceWithIsom) -> QQFieldElem

Given a quadratic space with isometry $(V, f)$, return the discriminant
of the underlying space $V$.

See [`discriminant(::AbstractSpace)`](@ref).

# Examples
```jldoctest
julia> V = quadratic_space(QQ, 2);

julia> Vf = quadratic_space_with_isometry(V; neg=true);

julia> discriminant(Vf)
-1
```
"""
discriminant(Vf::QuadSpaceWithIsom) = discriminant(space(Vf))

@doc raw"""
    is_positive_definite(Vf::QuadSpaceWithIsom) -> Bool

Given a quadratic space with isometry $(V, f)$, return whether the underlying
space $V$ is positive definite.

See [`is_positive_definite(::AbstractSpace)`](@ref).

# Examples
```jldoctest
julia> V = quadratic_space(QQ, 2);

julia> Vf = quadratic_space_with_isometry(V; neg=true);

julia> is_positive_definite(Vf)
true
```
"""
is_positive_definite(Vf::QuadSpaceWithIsom) = is_positive_definite(space(Vf))

@doc raw"""
    is_negative_definite(Vf::QuadSpaceWithIsom) -> Bool

Given a quadratic space with isometry $(V, f)$, return whether the underlying
space $V$ is negative definite.

See [`is_negative_definite(::AbstractSpace)`](@ref).

# Examples
```jldoctest
julia> V = quadratic_space(QQ, 2);

julia> Vf = quadratic_space_with_isometry(V; neg=true);

julia> is_negative_definite(Vf)
false
```
"""
is_negative_definite(Vf::QuadSpaceWithIsom) = is_negative_definite(space(Vf))

@doc raw"""
    is_definite(Vf::QuadSpaceWithIsom) -> Bool

Given a quadratic space with isometry $(V, f)$, return whether the underlying
space $V$ is definite.

See [`is_definite(::AbstractSpace)`](@ref).

# Examples
```jldoctest
julia> V = quadratic_space(QQ, 2);

julia> Vf = quadratic_space_with_isometry(V; neg=true);

julia> is_definite(Vf)
true
```
"""
is_definite(Vf::QuadSpaceWithIsom) = is_definite(space(Vf))

@doc raw"""
    diagonal(Vf::QuadSpaceWithIsom) -> Vector{QQFieldElem}

Given a quadratic space with isometry $(V, f)$, return the diagonal of the
underlying space $V$.

See [`diagonal(::AbstractSpace)`](@ref).

# Examples
```jldoctest
julia> V = quadratic_space(QQ, 2);

julia> Vf = quadratic_space_with_isometry(V; neg=true);

julia> diagonal(Vf)
2-element Vector{QQFieldElem}:
 1
 1
```
"""
diagonal(Vf::QuadSpaceWithIsom) = diagonal(space(Vf))

@doc raw"""
    signature_tuple(Vf::QuadSpaceWithIsom) -> Tuple{Int, Int, Int}

Given a quadratic space with isometry $(V, f)$, return the signature
tuple of the underlying space $V$.

# Examples
```jldoctest
julia> V = quadratic_space(QQ, 2);

julia> Vf = quadratic_space_with_isometry(V; neg=true);

julia> signature_tuple(Vf)
(2, 0, 0)
```
"""
signature_tuple(Vf::QuadSpaceWithIsom) = signature_tuple(space(Vf))

###############################################################################
#
#  Constructors
#
###############################################################################

@doc raw"""
    quadratic_space_with_isometry(
      V:QuadSpace,
      f::QQMatrix;
      check::Bool=false
    ) -> QuadSpaceWithIsom

Given a quadratic space $V$ and a matrix $f$, if $f$ defines an isometry of $V$
of order $n$ (possibly infinite), return the corresponding quadratic space with
isometry pair $(V, f)$.

# Examples
```jldoctest
julia> V = quadratic_space(QQ, QQ[ 2 -1;
                                  -1  2])
Quadratic space of dimension 2
  over rational field
with gram matrix
[ 2   -1]
[-1    2]

julia> f = matrix(QQ, 2, 2, [1  1;
                             0 -1])
[1    1]
[0   -1]

julia> Vf = quadratic_space_with_isometry(V, f)
Quadratic space of dimension 2
  with isometry of finite order 2
  given by
  [1    1]
  [0   -1]
```
"""
function quadratic_space_with_isometry(
    V::Hecke.QuadSpace,
    f::QQMatrix;
    check::Bool=true
  )
  if rank(V) == 0
    return QuadSpaceWithIsom(V, zero_matrix(QQ, 0, 0), -1)
  end

  if check
    @req !is_zero(f) "f must be non-zero"
    @req det(f) != 0 "Matrix must be invertible"
    @req f*gram_matrix(V)*transpose(f) == gram_matrix(V) "Matrix does not define an isometry of the given quadratic space"
  end

  n = multiplicative_order(f)
  return QuadSpaceWithIsom(V, f, n)
end

@doc raw"""
    quadratic_space_with_isometry(
      V::QuadSpace;
      neg::Bool=false
    ) -> QuadSpaceWithIsom

Given a quadratic space $V$, return the quadratic space with isometry pair
$(V, f)$ where $f$ is represented by the identity matrix.

If `neg` is set to `true`, then the isometry $f$ is negative the identity on
$V$.

# Examples
```jldoctest
julia> V = quadratic_space(QQ, QQ[ 2 -1;
                                  -1  2])
Quadratic space of dimension 2
  over rational field
with gram matrix
[ 2   -1]
[-1    2]

julia> Vf = quadratic_space_with_isometry(V)
Quadratic space of dimension 2
  with isometry of finite order 1
  given by
  [1   0]
  [0   1]
```
"""
function quadratic_space_with_isometry(V::Hecke.QuadSpace; neg::Bool=false)
  f = identity_matrix(QQ, dim(V))
  f = neg ? -f : f
  return quadratic_space_with_isometry(V, f; check=false)
end

###############################################################################
#
#  Operations on quadratic space with isometry
#
###############################################################################

@doc raw"""
    rescale(Vf::QuadSpaceWithIsom, a::RationalUnion)

Given a quadratic space with isometry $(V, f)$, return the pair $(V^a, f$)
where $V^a$ is the same space as $V$ with the associated quadratic form rescaled
by $a$.

# Examples
```jldoctest
julia> V = quadratic_space(QQ, QQ[ 2 -1;
                                  -1  2])
Quadratic space of dimension 2
  over rational field
with gram matrix
[ 2   -1]
[-1    2]

julia> Vf = quadratic_space_with_isometry(V)
Quadratic space of dimension 2
  with isometry of finite order 1
  given by
  [1   0]
  [0   1]

julia> Vf2 = rescale(Vf, 1//2)
Quadratic space of dimension 2
  with isometry of finite order 1
  given by
  [1   0]
  [0   1]

julia> space(Vf2)
Quadratic space of dimension 2
  over rational field
with gram matrix
[    1   -1//2]
[-1//2       1]
```
"""
function rescale(Vf::QuadSpaceWithIsom, a::RationalUnion)
  return quadratic_space_with_isometry(rescale(space(Vf), a), isometry(Vf); check=false)
end

@doc raw"""
    ^(Vf::QuadSpaceWithIsom, n::Int) -> QuadSpaceWithIsom

Given a quadratic space with isometry $(V, f)$ and an integer $n$, return the
pair $(V, f^n)$.

# Examples
```jldoctest
julia> V = quadratic_space(QQ, QQ[ 2 -1;
                                  -1  2])
Quadratic space of dimension 2
  over rational field
with gram matrix
[ 2   -1]
[-1    2]

julia> f = matrix(QQ, 2, 2, [1  1;
                             0 -1])
[1    1]
[0   -1]

julia> Vf = quadratic_space_with_isometry(V, f)
Quadratic space of dimension 2
  with isometry of finite order 2
  given by
  [1    1]
  [0   -1]

julia> Vf^2
Quadratic space of dimension 2
  with isometry of finite order 1
  given by
  [1   0]
  [0   1]
```
"""
function Base.:^(Vf::QuadSpaceWithIsom, n::Int)
  return quadratic_space_with_isometry(space(Vf), isometry(Vf)^n; check=false)
end

@doc raw"""
    direct_sum(
      x::Union{Vector{QuadSpaceWithIsom}, Vararg{QuadSpaceWithIsom}}
    ) -> QuadSpaceWithIsom, Vector{AbstractSpaceMor}, Vector{AbstractSpaceMor}

Given a finite collection of quadratic spaces with isometries
$(V_1, f_1), \ldots, (V_n, f_n)$, return the quadratic space with isometry
$(V, f)$ together with the embeddings of space $V_i \to V$ and the projections,
of $\mathbb{Q}$-vector spaces, $V\to V_i$.
Here $V$ is the direct sum of spaces $V := V_1 \oplus \ldots \oplus V_n$ and
$f$ is the isometry of $V$ induced by the diagonal actions of the $f_i$'s.

# Examples
```jldoctest
julia> V1 = quadratic_space(QQ, QQ[2 5;
                                   5 6])
Quadratic space of dimension 2
  over rational field
with gram matrix
[2   5]
[5   6]

julia> Vf1 = quadratic_space_with_isometry(V1; neg=true)
Quadratic space of dimension 2
  with isometry of finite order 2
  given by
  [-1    0]
  [ 0   -1]

julia> V2 = quadratic_space(QQ, QQ[ 2 -1;
                                   -1  2])
Quadratic space of dimension 2
  over rational field
with gram matrix
[ 2   -1]
[-1    2]

julia> f = matrix(QQ, 2, 2, [1  1;
                             0 -1])
[1    1]
[0   -1]

julia> Vf2 = quadratic_space_with_isometry(V2, f)
Quadratic space of dimension 2
  with isometry of finite order 2
  given by
  [1    1]
  [0   -1]

julia> Vf3, _, _ = direct_sum(Vf1, Vf2);

julia> Vf3
Quadratic space of dimension 4
  with isometry of finite order 2
  given by
  [-1    0   0    0]
  [ 0   -1   0    0]
  [ 0    0   1    1]
  [ 0    0   0   -1]

julia> space(Vf3)
Quadratic space of dimension 4
  over rational field
with gram matrix
[2   5    0    0]
[5   6    0    0]
[0   0    2   -1]
[0   0   -1    2]
```
"""
function direct_sum(x::Vector{T}) where T <: QuadSpaceWithIsom
  V, inj, proj = Hecke._biproduct(space.(x))
  f = block_diagonal_matrix(isometry.(x))
  return quadratic_space_with_isometry(V, f; check=false), inj
end

direct_sum(x::Vararg{QuadSpaceWithIsom}) = direct_sum(collect(x))

###############################################################################
#
#  Equality and hash
#
###############################################################################

function Base.:(==)(V1::QuadSpaceWithIsom, V2::QuadSpaceWithIsom)
  space(V1) == space(V2) || return false
  return isometry(V1) == isometry(V2)
end

function Base.hash(V::QuadSpaceWithIsom, u::UInt)
  u = Base.hash(space(V), u)
  return Base.hash(isometry(V), u)
end

###############################################################################
#
#  Spinor norm
#
###############################################################################

@doc raw"""
    rational_spinor_norm(Vf::QuadSpaceWithIsom; b::Int=-1) -> QQFieldElem

Given a rational quadratic space with isometry $(V, b, f)$, return the rational
spinor norm of $f$.

If $\Phi$ is the form on $V$, then the spinor norm is computed with respect to
$b\Phi$.
"""
function rational_spinor_norm(Vf::QuadSpaceWithIsom; b::Int=-1)
  @req dim(Vf) > 0 "V must have positive dimension"
  D, U = Hecke._gram_schmidt(gram_matrix(Vf), QQ)
  fD = U*isometry(Vf)*inv(U)
  return spin(b*D, fD)
end

###############################################################################
#
#  Useful
#
###############################################################################

function to_oscar(io::IO, Vf::QuadSpaceWithIsom)
  V = space(Vf)
  f = isometry(Vf)
  println(io, "G = matrix(QQ, $(dim(V)), $(dim(V)), ", gram_matrix(V), ");")
  println(io, "V = quadratic_space(QQ, G);")
  println(io, "f = matrix(QQ, $(dim(V)), $(dim(V)), ", f, ");")
  println(io, "Vf = quadratic_space_with_isometry(V, f);")
end

to_oscar(Vf::QuadSpaceWithIsom) = to_oscar(stdout, Vf)
