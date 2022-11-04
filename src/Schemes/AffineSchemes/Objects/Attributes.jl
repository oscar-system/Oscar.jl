export OO, coordinate_ring, base_ring, dim, codim, name

export ambient_affine_space, ambient_coordinate_ring, ambient_coordinates, ambient_embedding

export ring_type, base_ring_type, base_ring_elem_type, poly_type, ring_type



########################################################################
# (1) Attributes of AbsSpec
########################################################################

# Here is the inferface for AbsSpec

@Markdown.doc """
    coordinate_ring(X::AbsSpec)

On an affine scheme ``X = Spec(R)`` this returns the ring ``R``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> coordinate_ring(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field
```
We allow the shortcut `OO`
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> OO(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field
"""
coordinate_ring(X::AbsSpec) = OO(X)

@Markdown.doc """
    OO(X::AbsSpec)

On an affine scheme ``X = Spec(R)`` this returns the ring ``R``.
"""
function OO(X::AbsSpec{BRT, RT}) where {BRT, RT}
  OO(underlying_scheme(X))::RT
end


@Markdown.doc """
    ambient_affine_space(X::AbsSpec)

Return the ambient affine space of ``X``.

# Examples
```jldoctest
julia> X = affine_space(QQ, [:x,:y])
Spec of Multivariate Polynomial Ring in x, y over Rational Field

julia> ambient_affine_space(X) === X
true

julia> (x, y) = coordinates(X);

julia> Y = subscheme(X, [x])
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x, y)

julia> X === ambient_affine_space(Y)
true

julia> Z = subscheme(Y, y)
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x, y)

julia> ambient_affine_space(Z) === X
true

julia> V = hypersurface_complement(Y, y)

julia> ambient_affine_space(U) === X
true
```

We can create ``X``, ``Y`` and ``Z`` also by first constructing the corresponding
coordinate rings. The subset relations are inferred from the coordinate rings.
More precisely, for a polynomial ring ``P`` an ideal ``I < P `` and a multiplicatively closed subset
``U`` of ``P`` let ``R`` be one of ``P``, ``U^{-1}P``, ``P/I`` or ``U^{-1}(P/I)``.
In each case the ambient affine space is given by ``Spec(P)``.

# Examples
```jldoctest
julia> P, (x, y) = PolynomialRing(QQ, [:x, :y])

julia> X = Spec(P)

julia> I = ideal(P, x)
ideal(x)

julia> RmodI, quotient_map = quo(P, I);

julia> Y = Spec(RmodI)
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x)

julia> ambient_affine_space(Y) == X
true

julia> J = ideal(RmodI, y);

julia> RmodJ, quotient_map2 = quo(RmodI, J);

julia> Z = Spec(RmodJ)
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x, y)

julia> ambient_space(Z) == X
true

julia> U = powers_of_element(y)
powers of fmpq_mpoly[y]

julia> URmodI, _ = localization(RmodI, U);

julia> V = Spec(URmodI)
Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x) at the multiplicative set powers of fmpq_mpoly[y]

julia> ambient_affine_space(V) == X
true
```
"""
function ambient_affine_space(X::AbsSpec)
  error("$X does not have an ambient affine space")
end

function ambient_affine_space(X::AbsSpec{BRT, RT}) where {BRT, RT<:MPolyRing}
  return X
end

@attr function ambient_affine_space(X::AbsSpec{BRT,RT}) where {BRT, RT <: Union{MPolyQuo,MPolyLocalizedRing,MPolyQuoLocalizedRing}}
  return Spec(ambient_coordinate_ring(X))
end

@Markdown.doc """
    ambient_embedding(X::AbsSpec)

Return the embedding of ``X`` in its ambient affine space.

# Examples
```jldoctest
julia> X = affine_space(QQ, [:x,:y])
Spec of Multivariate Polynomial Ring in x, y over Rational Field

julia> (x, y) = coordinates(X);

julia> Y = subscheme(X, [x]);

julia> inc = ambient_embedding(Y)

julia> inc == inclusion_morphism(Y, X)
true
```
"""
ambient_embedding(X::AbsSpec) = inclusion_morphism(X, ambient_affine_space(X), check=false)

@Markdown.doc """
    ambient_coordinate_ring(X::AbsSpec)

Return the coordinate ring of the ambient affine space of ``X``.

See also [`ambient_affine_space](@ref).

# Examples
```jldoctest
julia> X = affine_space(QQ, [:x,:y])
Spec of Multivariate Polynomial Ring in x, y over Rational Field

julia> (x,y) = coordinates(X);

julia> Y = subscheme(X, [x])
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x, y)

julia> ambient_coordinate_ring(Y)
Multivariate Polynomial Ring in x, y over Rational Field
```
"""
function ambient_coordinate_ring(X::AbsSpec)
  return ambient_coordinate_ring(underlying_scheme(X))::MPolyRing
end

@Markdown.doc """
    ambient_coordinates(X::AbsSpec)

Return the coordinate functions of the ambient affine space of ``X``.

See also [`ambient_affine_space](@ref).

# Examples
```jldoctest
julia> X = affine_space(QQ, [:x,:y])
Spec of Multivariate Polynomial Ring in x, y over Rational Field

julia> (x,y) = coordinates(X);

julia> Y = subscheme(X, [x])
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x, y)

julia> (x,y) == ambient_coordinates(Y)
true
```
"""
ambient_coordinates(X::AbsSpec) = gens(ambient_coordinate_ring(X))

@Markdown.doc """
    coordinates(X::AbsSpec)

Return the coordinate functions of ``X`` as elements of its coordinate ring.

For ``X`` a subscheme of an ambient affine space, the coordinate functions are induced
by the ambient affine space.

# Examples
```jldoctest
julia> X = affine_space(QQ, [:x,:y])
Spec of Multivariate Polynomial Ring in x, y over Rational Field

julia> (x, y) = coordinates(X);
2-element Vector{fmpq_mpoly}:
 x
 y

julia> Y = subscheme(X, [x])
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x, y)

julia> (xY, yY) = coordinates(Y)
2-element Vector{MPolyQuoElem{fmpq_mpoly}}:
 x
 y

julia> parent(xY) == coordinate_ring(Y)
true
```
"""
coordinates(X::AbsSpec) = gens(OO(X))

ambient_coordinates(X::AbsSpec) = gens(ambient_coordinate_ring(X))

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

Return the dimension the affine scheme ``X = Spec(R)``.

By definition, this is the Krull dimension of ``R``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> dim(X)
3

julia> Y = affine_space(ZZ, 2)
Spec of Multivariate Polynomial Ring in x1, x2 over Integer Ring

julia> dim(Y) # one dimension comes from ZZ and one from x1
2
```
"""
@attr function dim(X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing})
  return dim(saturated_ideal(modulus(OO(X))))
end

@attr function dim(X::AbsSpec{<:Ring, <:MPolyLocalizedRing})
  # the following line is supposed to refer the problem to the
  # algebra side of the problem
  return dim(ideal(ambient_coordinate_ring(X), [zero(ambient_coordinate_ring(X))]))
end

@attr function dim(X::AbsSpec{<:Ring, <:MPolyRing})
  return dim(ideal(ambient_coordinate_ring(X), [zero(ambient_coordinate_ring(X))]))
end

@attr function dim(X::AbsSpec{<:Ring, <:MPolyQuo})
  return dim(modulus(OO(X)))
end


@Markdown.doc """
    codim(X::Spec)

Return the codimension of ``X`` in its ambient affine space.

Throws and error if ``X`` does not have an ambient affine space.

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
  return dim(ideal(ambient_coordinate_ring(X), [zero(ambient_coordinate_ring(X))])) - dim(X)
end


@doc Markdown.doc"""
    name(X::Spec)

Return the current name of an affine scheme.

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
    ambient_closure_ideal(X::AbsSpec{<:Any, <:MPolyRing})

Return the defining ideal of the closure of ``X`` in its ambient affine space.

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

julia> Y = subscheme(X, ideal(R, [x1*x2]))
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1*x2)

julia> I = Oscar.ambient_closure_ideal(Y)
ideal(x1*x2)

julia> base_ring(I) == OO(Y)
false

julia> base_ring(I) == R
true
```
"""
@attr ambient_closure_ideal(X::AbsSpec{<:Any, <:MPolyRing}) = ideal(OO(X), [zero(OO(X))])
ambient_closure_ideal(X::AbsSpec{<:Any, <:MPolyQuo}) = modulus(OO(X))
@attr ambient_closure_ideal(X::AbsSpec{<:Any, <:MPolyLocalizedRing}) = ideal(ambient_coordinate_ring(X), [zero(ambient_coordinate_ring(X))])
ambient_closure_ideal(X::AbsSpec{<:Any, <:MPolyQuoLocalizedRing}) = saturated_ideal(modulus(OO(X)))


########################################################################
# (3) Implementation of the AbsSpec interface for the basic Spec
########################################################################

OO(X::Spec) = X.OO
base_ring(X::Spec) = X.kk
ambient_coordinate_ring(X::Spec{<:Any, <:MPolyRing}) = OO(X)
ambient_coordinate_ring(X::Spec{<:Any, <:MPolyQuo}) = base_ring(OO(X))
ambient_coordinate_ring(X::Spec{<:Any, <:MPolyLocalizedRing}) = base_ring(OO(X))
ambient_coordinate_ring(X::Spec{<:Any, <:MPolyQuoLocalizedRing}) = base_ring(OO(X))
ambient_coordinate_ring(X::Spec{T, T}) where {T<:Field} = base_ring(X)



########################################################################
# (4) Type getters
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

