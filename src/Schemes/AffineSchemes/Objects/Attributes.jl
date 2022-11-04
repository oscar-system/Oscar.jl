export OO, ambient_ring, base_ring, dim, codim, name

export defining_ideal, strict_modulus

export ring_type, base_ring_type, base_ring_elem_type, poly_type, ring_type



########################################################################
# (1) Attributes of AbsSpec
########################################################################

# Here is the inferface for AbsSpec

@Markdown.doc """
    OO(X::AbsSpec)

On an affine scheme ``X = Spec(R)`` this returns the ring ``R``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> OO(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field
```
"""
function OO(X::AbsSpec{BRT, RT}) where {BRT, RT} 
  OO(underlying_scheme(X))::RT
end


@Markdown.doc """
    ambient_ring(X::AbsSpec)

On an affine scheme ``X = Spec(R)`` over ``𝕜`` this returns a 
polynomial ring ``P = 𝕜[x₁,…,xₙ]`` with natural coercion 
``P → R`` and the property that for every other (commutative) 
ring ``S`` and any homomorphism ``φ : R → S`` there is a morphism 
``ψ : P → S`` factoring through ``φ`` and such that ``φ`` 
is uniquely determined by ``ψ``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> ambient_ring(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field
```
"""
function ambient_ring(X::AbsSpec)
  return ambient_ring(underlying_scheme(X))::MPolyRing
end


@Markdown.doc """
    base_ring(X::AbsSpec)

On an affine scheme ``X/𝕜`` over ``𝕜`` this returns the ring ``𝕜``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> base_ring(X)
Rational Field
```
"""
function base_ring(X::AbsSpec{BRT, RT}) where {BRT, RT}
  return base_ring(underlying_scheme(X))::BRT
end


@Markdown.doc """
    dim(X::AbsSpec)

This method returns the dimension of an affine scheme ``X = Spec(R)``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> dim(X)
3
```
"""
@attr function dim(X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing})
  return dim(saturated_ideal(modulus(OO(X))))
end

@attr function dim(X::AbsSpec{<:Ring, <:MPolyLocalizedRing})
  # the following line is supposed to refer the problem to the
  # algebra side of the problem
  return dim(ideal(ambient_ring(X), [zero(ambient_ring(X))]))
end

@attr function dim(X::AbsSpec{<:Ring, <:MPolyRing})
  return dim(ideal(ambient_ring(X), [zero(ambient_ring(X))]))
end

@attr function dim(X::AbsSpec{<:Ring, <:MPolyQuo})
  return dim(modulus(OO(X)))
end


@Markdown.doc """
    codim(X::Spec)

In Oscar, we can compute for an affine scheme ``X``
the ring ``R = ambient_ring(X)``. This allows to consider
an embedding of ``X`` into ``Spec(R)``. This method returns
the codimension of the image of this embedding in ``Spec(R)``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> codim(X)
0

julia> R = OO(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> (x1,x2,x3) = gens(R)
3-element Vector{fmpq_mpoly}:
 x1
 x2
 x3

julia> Y = subscheme(X, x1)
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1)

julia> codim(Y)
1
```
"""
@attr function codim(X::Spec)
  return dim(ideal(ambient_ring(X), [zero(ambient_ring(X))])) - dim(X)
end


@doc Markdown.doc"""
    name(X::Spec)

Returns the current name of an affine scheme.
This name can be specified via `set_name!`.

# Examples
```jldoctest
julia> X = affine_space(QQ, 3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> name(X)
"unnamed affine variety"

julia> set_name!(X, "affine 3-dimensional space")

julia> name(X)
"affine 3-dimensional space"
```
"""
@attr String function name(X::Spec)
  return "unnamed affine variety"
end


function set_name!(X::AbsSpec, name::String)
  return set_attribute!(X, :name, name)
end


########################################################################
# (2) Further attributes for AbsSpec
########################################################################

# TODO: Needed?

@Markdown.doc """
    defining_ideal(X::AbsSpec{<:Any, <:MPolyRing})

This method return the defining ideal of an affine scheme ``X``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> R = OO(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> (x1,x2,x3) = gens(R)
3-element Vector{fmpq_mpoly}:
 x1
 x2
 x3

julia> Y = subscheme(X,ideal(R,[x1*x2]))
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1*x2)

julia> defining_ideal(Y)
ideal(x1*x2)
```
"""
@attr defining_ideal(X::AbsSpec{<:Any, <:MPolyRing}) = ideal(OO(X), [zero(OO(X))])
defining_ideal(X::AbsSpec{<:Any, <:MPolyQuo}) = modulus(OO(X))
@attr defining_ideal(X::AbsSpec{<:Any, <:MPolyLocalizedRing}) = ideal(OO(X), [zero(OO(X))])
defining_ideal(X::AbsSpec{<:Any, <:MPolyQuoLocalizedRing}) = modulus(OO(X))


@Markdown.doc """
    strict_modulus(X::AbsSpec)

This method return the strict modulus of an affine scheme ``X``.
This is the ideal ``I`` in the `ambient_ring` of ``X`` consisting 
of elements ``f ∈ R`` such that their restriction to ``X`` vanishes.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> R = OO(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> (x1,x2,x3) = gens(R)
3-element Vector{fmpq_mpoly}:
 x1
 x2
 x3

julia> Y = subscheme(X,ideal(R,[x1*x2]))
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1*x2)

julia> strict_modulus(Y)
ideal(x1*x2)
```
"""
strict_modulus(X::AbsSpec) = saturated_ideal(modulus(OO(X)))



########################################################################
# (3) Further attributes for Spec
########################################################################

strict_modulus(X::Spec) = saturated_ideal(modulus(OO(X)))



########################################################################
# (4) Implementation of the AbsSpec interface for the basic Spec
########################################################################

OO(X::Spec) = X.OO
base_ring(X::Spec) = X.kk
ambient_ring(X::Spec{<:Any, <:MPolyRing}) = OO(X)
ambient_ring(X::Spec{<:Any, <:MPolyQuo}) = base_ring(OO(X))
ambient_ring(X::Spec{<:Any, <:MPolyLocalizedRing}) = base_ring(OO(X))
ambient_ring(X::Spec{<:Any, <:MPolyQuoLocalizedRing}) = base_ring(OO(X))
ambient_ring(X::Spec{T, T}) where {T<:Field} = base_ring(X)



########################################################################
# (5) Type getters
########################################################################

# TODO: Needed?

ring_type(::Type{SpecType}) where {BRT, RT, SpecType<:AbsSpec{BRT, RT}} = RT
ring_type(X::AbsSpec) = ring_type(typeof(X))

base_ring_type(::Type{SpecType}) where {BRT, RT, SpecType<:AbsSpec{BRT, RT}} = BRT
base_ring_type(X::AbsSpec) = base_ring_type(typeof(X))
base_ring_elem_type(::Type{SpecType}) where {BRT, RT, SpecType<:AbsSpec{BRT, RT}} = elem_type(BRT)
base_ring_elem_type(X::AbsSpec) = base_ring_elem_type(typeof(X))

poly_type(::Type{SpecType}) where {BRT, RT<:MPolyRing, SpecType<:AbsSpec{BRT, RT}} = elem_type(RT)
poly_type(::Type{SpecType}) where {BRT, T, RT<:MPolyQuo{T}, SpecType<:AbsSpec{BRT, RT}} = T
poly_type(::Type{SpecType}) where {BRT, T, RT<:MPolyLocalizedRing{<:Any, <:Any, <:Any, T}, SpecType<:AbsSpec{BRT, RT}} = T
poly_type(::Type{SpecType}) where {BRT, T, RT<:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, T}, SpecType<:AbsSpec{BRT, RT}} = T
poly_type(X::AbsSpec) = poly_type(typeof(X))

ring_type(::Type{Spec{BRT, RT}}) where {BRT, RT} = RT
ring_type(X::Spec) = ring_type(typeof(X))
base_ring_type(::Type{Spec{BRT, RT}}) where {BRT, RT} = BRT
base_ring_type(X::Spec) = base_ring_type(typeof(X))

