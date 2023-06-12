
###############################################################################
#
#  Accesors
#
###############################################################################

@doc raw"""
    space(Vf::QuadSpaceWithIsom) -> QuadSpace

Given a quadratic space with isometry $(V, f)$, return the underlying space `V`.
"""
space(Vf::QuadSpaceWithIsom) = Vf.V

@doc raw"""
    isometry(Vf::QuadSpaceWithIsom) -> QQMatrix

Given a quadratic space with isometry $(V, f)$, return the underlying isometry
`f`.
"""
isometry(Vf::QuadSpaceWithIsom) = Vf.f

@doc raw"""
    order_of_isometry(Vf::QuadSpaceWithIsom) -> IntExt

Given a quadratic space with isometry $(V, f)$, return the order of the
underlying isometry `f`
"""
order_of_isometry(Vf::QuadSpaceWithIsom) = Vf.n

###############################################################################
#
#  Attributes
#
###############################################################################

@doc raw"""
    rank(Vf::QuadSpaceWithIsom) -> Integer

Given a quadratic space with isometry $(V, f)$, return the rank of the underlying
space `V`.
"""
rank(Vf::QuadSpaceWithIsom) = rank(space(Vf))::Integer

@doc raw"""
    dim(Vf::QuadSpaceWithIsom) -> Integer

Given a quadratic space with isometry $(V, f)$, return the dimension of the
underlying space of `V`
"""
dim(Vf::QuadSpaceWithIsom) = dim(space(Vf))::Integer

@doc raw"""
    charpoly(Vf::QuadSpaceWithIsom) -> QQPolyRingElem

Given a quadratic space with isometry $(V, f)$, return the characteristic
polynomial of the underlying isometry `f`
"""
charpoly(Vf::QuadSpaceWithIsom) = charpoly(isometry(Vf))::QQPolyRingElem

@doc raw"""
    minpoly(Vf::QuadSpaceWithIsom) -> QQPolyRingElem

Given a quadratic space with isometry $(V, f)$, return the minimal
polynomial of the underlying isometry `f`.
"""
minpoly(Vf) = minpoly(isometry(Vf))::QQPolyRingElem

@doc raw"""
    gram_matrix(Vf::QuadSpaceWithIsom) -> QQMatrix

Given a quadratic space with isometry $(V, f)$, return the Gram matrix
of the underlying space `V` with respect to its standard basis.
"""
gram_matrix(Vf::QuadSpaceWithIsom) = gram_matrix(space(Vf))::QQMatrix

@doc raw"""
    det(Vf::QuadSpaceWithIsom) -> QQFieldElem

Given a quadratic space with isometry $(V, f)$, return the determinant
of the underlying space `V`.
"""
det(Vf::QuadSpaceWithIsom) = det(space(Vf))::QQFieldElem

@doc raw"""
    discriminant(Vf::QuadSpaceWithIsom) -> QQFieldElem

Given a quadratic space with isometry $(V, f)$, return the discriminant
of the underlying space `V`.
"""
discriminant(Vf::QuadSpaceWithIsom) = discriminant(space(Vf))::QQFieldElem

@doc raw"""
    is_positive_definite(Vf::QuadSpaceWithIsom) -> Bool

Given a quadratic space with isometry $(V, f)$, return whether the underlying
space `V` is positive definite.
"""
is_positive_definite(Vf::QuadSpaceWithIsom) = is_positive_definite(space(Vf))::Bool

@doc raw"""
    is_negative_definite(Vf::QuadSpaceWithIsom) -> Bool

Given a quadratic space with isometry $(V, f)$, return whether the underlying
space `V` is negative definite.
"""
is_negative_definite(Vf::QuadSpaceWithIsom) = is_negative_definite(space(Vf))::Bool

@doc raw"""
    is_definite(Vf::QuadSpaceWithIsom) -> Bool

Given a quadratic space with isometry $(V, f)$, return whether the underlying
space `V` is definite.
"""
is_definite(Vf::QuadSpaceWithIsom) = is_definite(space(Vf))::Bool

@doc raw"""
    diagonal(Vf::QuadSpaceWithIsom) -> Vector{QQFieldElem}

Given a quadratic space with isometry $(V, f)$, return the diagonal of the
underlying space `V`, that is a list of rational numbers $a_1, \ldots a_n$
such that `V` is isometric to the space whose Gram matrix is diagonal with
entries $a_1,\ldots, a_n$.
"""
diagonal(Vf::QuadSpaceWithIsom) = diagonal(space(Vf))::Vector{QQFieldElem}

@doc raw"""
    signature_tuple(Vf::QuadSpaceWithIsom) -> Tuple{Int, Int, Int}

Given a quadratic space with isometry $(V, f)$, return the signature
tuple of the underlying space `V`.
"""
signature_tuple(Vf::QuadSpaceWithIsom) = signature_tuple(space(Vf))::Tuple{Int, Int, Int}

###############################################################################
#
#  Constructors
#
###############################################################################

@doc raw"""
    quadratic_space_with_isometry(V:QuadSpace, f::QQMatrix; check::Bool = false)
                                                            -> QuadSpaceWithIsom

Given a quadratic space `V` and a matrix `f`, if `f` defines an isometry of `V`
of order `n` (possibly infinite), return the corresponding quadratic space with
isometry pair $(V, f)$.
"""
function quadratic_space_with_isometry(V::Hecke.QuadSpace, f::QQMatrix;
                                                     check::Bool = true)
  if rank(V) == 0
    return QuadSpaceWithIsom(V, zero_matrix(QQ, 0, 0), -1)
  end

  if check
    @req det(f) != 0 "Matrix must be invertible"
    @req f*gram_matrix(V)*transpose(f) == gram_matrix(V) "Matrix does not define an isometry of the given quadratic space"
  end

  n = multiplicative_order(f)
  return QuadSpaceWithIsom(V, f, n)
end

@doc raw"""
    quadratic_space_with_isometry(V::QuadSpace; neg::Bool = false) -> QuadSpaceWithIsom

Given a quadratic space `V`, return the quadratic space with isometry pair $(V, f)$
where `f` is represented by the identity matrix.

If `neg` is set to true, then the isometry `f` is negative the identity on `V`.
"""
function quadratic_space_with_isometry(V::Hecke.QuadSpace; neg::Bool = false)
  f = identity_matrix(QQ, dim(V))
  f = neg ? -f : f
  return quadratic_space_with_isometry(V, f, check=false)
end

###############################################################################
#
#  Operations on quadratic space with isometry
#
###############################################################################

@doc raw"""
    rescale(Vf::QuadSpaceWithIsom, a::RationalUnion)

Given a quadratic space with isometry $(V, f)$, return the pair $(V^a, f$) where
$V^a$ is the same space as `V` with the associated quadratic form rescaled by `a`.
"""
function rescale(Vf::QuadSpaceWithIsom, a::Hecke.RationalUnion)
  return quadratic_space_with_isometry(rescale(space(Vf), a), isometry(Vf), check = false)
end

@doc raw"""
    direct_sum(x::Vector{QuadSpaceWithIsom}) -> QuadSpaceWithIsom, Vector{AbstractSpaceMor}
    direct_sum(x::Vararg{QuadSpaceWithIsom}) -> QuadSpaceWithIsom, Vector{AbstractSpaceMor}

Given a collection of quadratic spaces with isometries $(V_1, f_1) \ldots, (V_n, f_n)$,
return the quadratic space with isometry $(V, f)$ together with the injections
$V_i \to V$, where `V` is the direct sum $V := V_1 \oplus \ldots \oplus V_n$ and
`f` is the isometry of `V` induced by the diagonal actions of the $f_i$'s.

For objects of type `QuadSpaceWithIsom`, finite direct sums and finite direct products
agree and they are therefore called biproducts.
If one wants to obtain $(V, f)$ as a direct product with the projections $V \to V_i$,
one should call `direct_product(x)`.
If one wants to obtain $(V, f)$ as a biproduct with the injections $V_i \to V$ and
the projections $V \to V_i$, one should call `biproduct(x)`.
"""

function direct_sum(x::Vector{T}) where T <: QuadSpaceWithIsom
  V, inj = direct_sum(space.(x))
  f = block_diagonal_matrix(isometry.(x))
  return quadratic_space_with_isometry(V, f, check=false), inj
end

direct_sum(x::Vararg{QuadSpaceWithIsom}) = direct_sum(collect(x))

@doc raw"""
    direct_product(x::Vector{QuadSpaceWithIsom}) -> QuadSpaceWithIsom, Vector{AbstractSpaceMor}
    direct_product(x::Vararg{QuadSpaceWithIsom}) -> QuadSpaceWithIsom, Vector{AbstractSpaceMor}

Given a collection of quadratic spaces with isometries $(V_1, f_1), \ldots, (V_n, f_n)$,
return the quadratic space with isometry $(V, f)$ together with the projections
$V \to V_i$, where `V` is the direct product $V := V_1 \times \ldots \times V_n$ and
`f` is the isometry of `V` induced by the diagonal actions of the $f_i$'s.

For objects of type `QuadSpaceWithIsom`, finite direct sums and finite direct products
agree and they are therefore called biproducts.
If one wants to obtain $(V, f)$ as a direct sum with the injections $V_i \to V$,
one should call `direct_sum(x)`.
If one wants to obtain $(V, f)$ as a biproduct with the injections $V_i \to V$ and
the projections $V \to V_i$, one should call `biproduct(x)`.
"""

function direct_product(x::Vector{T}) where T <: QuadSpaceWithIsom
  V, proj = direct_product(space.(x))
  f = block_diagonal_matrix(isometry.(x))
  return quadratic_space_with_isometry(V, f, check=false), proj
end

direct_product(x::Vararg{QuadSpaceWithIsom}) = direct_product(collect(x))

@doc raw"""
    biproduct(x::Vector{QuadSpaceWithIsom}) -> QuadSpaceWithIsom, Vector{AbstractSpaceMor}, Vector{AbstractSpaceMor}
    biproduct(x::Vararg{QuadSpaceWithIsom}) -> QuadSpaceWithIsom, Vector{AbstractSpaceMor}, Vector{AbstractSpaceMor}

Given a collection of quadratic spaces with isometries $(V_1, f_1), \ldots, (V_n, f_n)$,
return the quadratic space with isometry $(V, f)$ together with the injections
$V_i \to V$ and the projections $V \to V_i$, where `V` is the biproduct
$V := V_1 \oplus \ldots \oplus V_n$ and `f` is the isometry of `V` induced by the
diagonal actions of the $f_i$'s.

For objects of type `QuadSpaceWithIsom`, finite direct sums and finite direct products
agree and they are therefore called biproducts.
If one wants to obtain $(V, f)$ as a direct sum with the injections $V_i \to V$,
one should call `direct_sum(x)`.
If one wants to obtain $(V, f)$ as a direct product with the projections $V \to V_i$,
one should call `direct_product(x)`.
"""

function biproduct(x::Vector{T}) where T <: QuadSpaceWithIsom
  V, inj, proj = biproduct(space.(x))
  f = block_diagonal_matrix(isometry.(x))
  return quadratic_space_with_isometry(V, f, check=false), inj, proj
end

biproduct(x::Vararg{QuadSpaceWithIsom}) = biproduct(collect(x))

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

