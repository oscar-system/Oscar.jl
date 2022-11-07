export affine_space
export closure
export hypersurface_complement
export standard_spec
export subscheme


########################################################
# (1) Generic constructors
########################################################

@doc Markdown.doc"""
    Spec(R::MPolyRing, I::MPolyIdeal)

Constructs the affine scheme of the ideal ``I`` in the ring ``R``.
This is the spectrum of the quotient ring ``R/I``.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);

julia> I = ideal(R, [x]);

julia> Spec(R, I)
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x)
```
"""
Spec(R::MPolyRing, I::MPolyIdeal) = Spec(quo(R, I)[1])


@doc Markdown.doc"""
    Spec(R::MPolyRing, U::AbsMPolyMultSet)

Given a polynomial ring ``R``, we can localize that polynomial
ring at a multiplicatively closed subset ``U`` of ``R``. The spectrum
of the localized ring $U^{-1} R$ is computed by this method.

# Examples
Set ``C`` to be the positive orthant in two dimensions.
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);

julia> I = ideal(R, [x]);

julia> U = complement_of_ideal(I);

julia> Spec(R, U)
Spec of localization of Multivariate Polynomial Ring in x, y over Rational Field at the complement of ideal(x)
```
"""
Spec(R::MPolyRing, U::AbsMPolyMultSet) = Spec(Localization(R, U)[1])


@doc Markdown.doc"""
    Spec(R::MPolyRing, I::MPolyIdeal, U::AbsMPolyMultSet)

We allow to combine quotients and localizations at the same time.
That is, consider a polynomial ring ``R``, an ideal ``I`` of ``R`` and
a multiplicatively closed subset ``U`` of ``R``. The spectrum of the
localized ring $U^{-1} (R/I)$ is computed by this method.

# Examples
Set ``C`` to be the positive orthant in two dimensions.
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);

julia> I = ideal(R, [x]);

julia> U = complement_of_ideal(ideal(R, [y]));

julia> Spec(R, I, U)
Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x) at the multiplicative set complement of ideal(y)
```
"""
Spec(R::MPolyRing, I::MPolyIdeal, U::AbsMPolyMultSet) = Spec(MPolyQuoLocalizedRing(R, I, U))



########################################################
# (2) Copy constructors
########################################################

@doc Markdown.doc"""
    Spec(X::Spec)

For convenience, an affine spectrum can be passed to `Spec`
to create a new spectrum. This can be particularly useful
when in need to copy an affine spectrum.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);

julia> I = ideal(R, [x]);

julia> X = Spec(R, I)
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x)

julia> Y = Spec(X)
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x)
```
"""
Spec(X::Spec) = Spec(OO(X))
Base.deepcopy_internal(X::Spec, dict::IdDict) = Spec(deepcopy_internal(OO(X), dict))



########################################################
# (3) Affine n-dimensional space
########################################################

@doc Markdown.doc"""
    affine_space(kk::BRT, n::Int; variable_name="x") where {BRT<:Ring}

The ``n``-dimensional affine space over a ring ``kk`` is created
by this method. By default, the variable names are chosen as $x_1$, $x_2$
and so on. This choice can be overwritten with a third optional argument.

# Examples
```jldoctest
julia> affine_space(QQ, 5)
Spec of Multivariate Polynomial Ring in x1, x2, x3, x4, x5 over Rational Field

julia> affine_space(QQ,5,variable_name="y")
Spec of Multivariate Polynomial Ring in y1, y2, y3, y4, y5 over Rational Field
```
"""
function affine_space(kk::BRT, n::Int; variable_name="x") where {BRT<:Ring}
  R, _ = PolynomialRing(kk, [variable_name * "$i" for i in 1:n])
  return Spec(R)
end


@doc Markdown.doc"""
    affine_space(kk::BRT, var_symbols::Vector{Symbol}) where {BRT<:Ring}

Creates the ``n``-dimensional affine space over a ring ``kk``,
but allows more flexibility in the choice of variable names.
The following example demonstrates this.

# Examples
```jldoctest
julia> affine_space(QQ,[:y1,:z2,:a])
Spec of Multivariate Polynomial Ring in y1, z2, a over Rational Field
```
"""
function affine_space(kk::BRT, var_symbols::Vector{Symbol}) where {BRT<:Ring}
  R, _ = PolynomialRing(kk, var_symbols)
  return Spec(R)
end


########################################################
# (4) StdSpec (needed?)
########################################################

@doc Markdown.doc"""
    standard_spec(X::AbsSpec)

For an affine spectrum with coordinate ring of type `MPolyRing`, 
`MPolyQuo`, or `MPolyLocalizedRing`, this returns the canonical 
transform to a `Spec` of an `MPolyQuoLocalizedRing`. 

# Examples
```jldoctest
julia> standard_spec(affine_space(QQ,5))
Spec of Localization of Quotient of Multivariate Polynomial Ring in x1, x2, x3, x4, x5 over Rational Field by ideal(0) at the multiplicative set powers of fmpq_mpoly[1]
```
"""
function standard_spec(X::AbsSpec)
  error("not implemented for input of type $(typeof(X))")
end

standard_spec(X::AbsSpec{<:Any, <:MPolyRing}) = Spec(MPolyQuoLocalizedRing(OO(X), ideal(OO(X), [zero(OO(X))]), units_of(OO(X))))


#@doc Markdown.doc"""
#    standard_spec(X::AbsSpec{<:Any, <:MPolyQuo})
#
#For an affine spectrum whose coordinate ring is the
#quotient of a polynomial ring, this method computes
#the standard spectrum.
#
## Examples
#```jldoctest
#julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);
#
#julia> I = ideal(R, [x]);
#
#julia> X = Spec(R, I)
#Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x)
#
#julia> standard_spec(X)
#Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x) at the multiplicative set powers of fmpq_mpoly[1]
#```
#"""
function standard_spec(X::AbsSpec{<:Any, <:MPolyQuo})
  A = OO(X)
  R = base_ring(A)
  return Spec(MPolyQuoLocalizedRing(R, modulus(A), units_of(R)))
end


#@doc Markdown.doc"""
#    standard_spec(X::AbsSpec{<:Any, <:MPolyLocalizedRing})
#
#For an affine spectrum whose coordinate ring is the
#quotient of a polynomial ring, this method computes
#the standard spectrum.
#
## Examples
#```jldoctest
#julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);
#
#julia> I = ideal(R, [x]);
#
#julia> U = complement_of_ideal(I);
#
#julia> X = Spec(R, U)
#Spec of localization of Multivariate Polynomial Ring in x, y over Rational Field at the complement of ideal(x)
#
#julia> standard_spec(X)
#Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(0) at the multiplicative set complement of ideal(x)
#```
#"""
standard_spec(X::AbsSpec{<:Any, <:MPolyLocalizedRing}) = Spec(MPolyQuoLocalizedRing(ambient_ring(X), ideal(ambient_ring(X), [zero(ambient_ring(X))]), inverted_set(OO(X))))


#@doc Markdown.doc"""
#    standard_spec(X::AbsSpec{<:Any, <:MPolyQuoLocalizedRing})
#
#For an affine spectrum whose coordinate ring is the
#quotient of a polynomial ring, this method computes
#the standard spectrum.
#
## Examples
#```jldoctest
#julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);
#
#julia> I = ideal(R, [x]);
#
#julia> U = complement_of_ideal(ideal(R, [y]));
#
#julia> X = Spec(R, I, U)
#Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x) at the multiplicative set complement of ideal(y)
#
#julia> standard_spec(X)
#Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x) at the multiplicative set complement of ideal(y)
#```
#"""
standard_spec(X::AbsSpec{<:Any, <:MPolyQuoLocalizedRing}) = Spec(OO(X))



########################################################
# (5) Closed subschemes
########################################################

@doc Markdown.doc"""
    subscheme(X::AbsSpec, f::Vector{<:RingElem})

For an affine spectrum ``X`` and elements ``f_1``, ``f_2``,
etc. of the coordinate ring of ``X``, this method computes
the subscheme ``V(f_1, f_2, \dots)`` of ``X``.

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

julia> subscheme(X,x1)
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1)

julia> subscheme(X,[x1,x2])
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1, x2)
```
"""
subscheme(X::AbsSpec, f::Vector{<:RingElem}) = subscheme(X, ideal(OO(X), f))
function subscheme(X::Spec, f::Vector{<:RingElem})
  all(x->(parent(x) == OO(X)), f) || return subscheme(X, OO(X).(f))
  return subscheme(X, ideal(OO(X), f))
end
subscheme(X::AbsSpec, f::RingElem) = subscheme(X, ideal(OO(X), [f]))


@Markdown.doc """
    subscheme(X::AbsSpec, I::Ideal)

For a scheme ``X = Spec(R)`` and an ideal ``I ⊂ 𝒪(X)``
this returns the closed subscheme defined by ``I``.

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

julia> subscheme(X,ideal(R,[x1*x2]))
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1*x2)
```
"""
function subscheme(X::AbsSpec, I::Ideal)
  base_ring(I) == OO(X) || return subscheme(X, ideal(OO(X), OO(X).(gens(I)))) # this will throw if coercion is not possible
  return Spec(quo(OO(X), I)[1])
end



########################################################
# (6) Open subschemes
########################################################

@Markdown.doc """
    hypersurface_complement(X::AbsSpec, f::RingElem)

For a scheme ``X = Spec(R)`` and an element ``f ∈ R``
this returns the open subscheme ``U = Spec(R[f⁻¹]) = X ∖ V(f)``
defined by the complement of the vanishing locus of ``f``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> R = OO(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> (x1, x2, x3) = gens(R)
3-element Vector{fmpq_mpoly}:
 x1
 x2
 x3

julia> hypersurface_complement(X, x1)
Spec of localization of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field at the powers of fmpq_mpoly[x1]
```
"""
function hypersurface_complement(X::AbsSpec, f::RingElem)
  return hypersurface_complement(underlying_scheme(X), f)::AbsSpec
end

function hypersurface_complement(X::SpecType, f::RingElem) where {SpecType<:AbsSpec{<:Any, <:MPolyQuoLocalizedRing}}
  parent(f) == OO(X) || return hypersurface_complement(X, OO(X)(f))
  h = lifted_numerator(f)
  U = MPolyPowersOfElement(h)
  W, _ = Localization(OO(X), U)
  return Spec(W)
end

function hypersurface_complement(X::SpecType, f::RingElem) where {SpecType<:AbsSpec{<:Any, <:MPolyLocalizedRing}}
  parent(f) == OO(X) || return hypersurface_complement(X, OO(X)(f))
  h = numerator(f)
  U = MPolyPowersOfElement(h)
  W, _ = Localization(OO(X), U)
  return Spec(W)
end

function hypersurface_complement(X::SpecType, f::RingElem) where {SpecType<:AbsSpec{<:Any, <:MPolyRing}}
  parent(f) == OO(X) || return hypersurface_complement(X, OO(X)(f))
  U = MPolyPowersOfElement(f)
  W, _ = Localization(OO(X), U)
  return Spec(W)
end

function hypersurface_complement(X::SpecType, f::RingElem) where {SpecType<:AbsSpec{<:Any, <:MPolyQuo}}
  parent(f) == OO(X) || return hypersurface_complement(X, OO(X)(f))
  U = MPolyPowersOfElement(lift(f))
  W, _ = Localization(OO(X), U)
  return Spec(W)
end


@Markdown.doc """
    hypersurface_complement(X::AbsSpec, f::Vector{<:RingElem})

For a scheme ``X = Spec(R)`` and elements ``f₁, f₂, ... ∈ R``
this returns the open subscheme ``U = Spec(R[f₁⁻¹,f₂⁻¹, ...]) = X ∖ V(f₁⋅f₂⋅…)``
defined by the complement of the vanishing locus of the product ``f₁⋅f₂⋅…``.

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

julia> hypersurface_complement(X,[x1,x2])
Spec of localization of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field at the powers of fmpq_mpoly[x1, x2]
```
"""
function hypersurface_complement(X::AbsSpec, f::Vector{<:RingElem})
  return hypersurface_complement(underlying_scheme(X), f)::AbsSpec
end

function hypersurface_complement(X::SpecType, f::Vector{<:RingElem}) where {SpecType<:AbsSpec{<:Any, <:MPolyQuoLocalizedRing}}
  all(x->(parent(x) == OO(X)), f) || return hypersurface_complement(X, OO(X).(f))
  h = lifted_numerator.(f)
  U = MPolyPowersOfElement(ambient_ring(X), h)
  W, _ = Localization(OO(X), U)
  return Spec(W)
end

function hypersurface_complement(X::SpecType, f::Vector{<:RingElem}) where {SpecType<:AbsSpec{<:Any, <:MPolyLocalizedRing}}
  all(x->(parent(x) == OO(X)), f) || return hypersurface_complement(X, OO(X).(f))
  h = numerator.(f)
  U = MPolyPowersOfElement(ambient_ring(X), h)
  W, _ = Localization(OO(X), U)
  return Spec(W)
end

function hypersurface_complement(X::SpecType, f::Vector{<:RingElem}) where {SpecType<:AbsSpec{<:Any, <:MPolyRing}}
  all(x->(parent(x) == OO(X)), f) || return hypersurface_complement(X, OO(X).(f))
  U = MPolyPowersOfElement(ambient_ring(X), f)
  W, _ = Localization(OO(X), U)
  return Spec(W)
end

function hypersurface_complement(X::SpecType, f::Vector{<:RingElem}) where {SpecType<:AbsSpec{<:Any, <:MPolyQuo}}
  all(x->(parent(x) == OO(X)), f) || return hypersurface_complement(X, OO(X).(f))
  U = MPolyPowersOfElement(ambient_ring(X), lift.(f))
  W, _ = Localization(OO(X), U)
  return Spec(W)
end



########################################################
# (7) Intersections of (affine) schemes
########################################################


# 7.1 Intersection methods involving an empty scheme
Base.intersect(E::EmptyScheme{BRT}, X::Scheme{BRT}) where {BRT<:Ring} = E
Base.intersect(X::Scheme{BRT}, E::EmptyScheme{BRT}) where {BRT<:Ring} = E
Base.intersect(X::EmptyScheme{BRT}, E::EmptyScheme{BRT}) where {BRT<:Ring} = E


# 7.2 Intersection methods not involving empty schemes

### For Specs of MPolyRings
# TODO  intersect X,Y for X<Y should return a copy of X with === ambient_rings
# Spec(X) does not apply for instance to principal open subsets hence a change
# is necessary
@Markdown.doc """
    Base.intersect(X::AbsSpec, Y::AbsSpec)

This method computes the intersection to two affine
schemes that reside in the same ambient affine space.

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

julia> Y1 = subscheme(X,[x1])
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1)

julia> Y2 = subscheme(X,[x2])
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x2)

julia> intersect(Y1, Y2)
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1, x2)
```
"""
function Base.intersect(X::AbsSpec{BRT, <:Ring}, Y::AbsSpec{BRT, <:Ring}) where {BRT<:Ring}
  error("method not implemted for arguments of type $(typeof(X)) and $(typeof(Y))")
end

function Base.intersect(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyRing}
  ) where {BRT<:Ring}
  R = OO(X)
  R == OO(Y) || error("schemes can not be compared")
  return Spec(X)
end


function Base.intersect(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyQuo}
  ) where {BRT<:Ring}
  R = OO(X)
  R === ambient_ring(Y) || error("schemes can not be compared")
  return Spec(Y)
end


function Base.intersect(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyLocalizedRing}
  ) where {BRT<:Ring}
  R = OO(X)
  R === ambient_ring(Y) || error("schemes can not be compared")
  return Spec(Y)
end


function Base.intersect(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing}
  ) where {BRT<:Ring}
  R = OO(X)
  R === ambient_ring(Y) || error("schemes can not be compared")
  return Spec(Y)
end


function Base.intersect(
    Y::AbsSpec{BRT, <:Ring},
    X::AbsSpec{BRT, <:MPolyRing}
  ) where {BRT<:Ring}
  return intersect(X, Y)
end


### For Specs of MPolyQuos
function Base.intersect(
    X::AbsSpec{BRT, <:MPolyQuo},
    Y::AbsSpec{BRT, <:MPolyQuo}
  ) where {BRT<:Ring}
  R = ambient_ring(X)
  R === ambient_ring(Y) || error("schemes can not be compared")
  return Spec(quo(R, modulus(OO(X)) + modulus(OO(Y)))[1])
end


function Base.intersect(
    X::AbsSpec{BRT, <:MPolyQuo},
    Y::AbsSpec{BRT, <:MPolyLocalizedRing}
  ) where {BRT<:Ring}
  R = ambient_ring(X)
  R === ambient_ring(Y) || error("schemes can not be compared")
  return Spec(quo(OO(Y), OO(Y)(modulus(OO(X))))[1])
end


function Base.intersect(
    X::AbsSpec{BRT, <:MPolyQuo},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing}
  ) where {BRT<:Ring}
  R = ambient_ring(X)
  R === ambient_ring(Y) || error("schemes can not be compared")
  return Spec(quo(OO(Y), OO(Y)(modulus(OO(X))))[1])
end


function Base.intersect(
    Y::AbsSpec{BRT, <:Ring},
    X::AbsSpec{BRT, <:MPolyQuo}
  ) where {BRT<:Ring}
  return intersect(X, Y)
end


### For Specs of MPolyLocalizedRings
function Base.intersect(
    X::AbsSpec{BRT, <:MPolyLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyLocalizedRing}
  ) where {BRT<:Ring}
  R = ambient_ring(X)
  R === ambient_ring(Y) || error("schemes can not be compared")
  return Spec(Localization(R, inverted_set(OO(X)) * inverted_set(OO(Y)))[1])
end


function Base.intersect(
    X::AbsSpec{BRT, <:MPolyLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing}
  ) where {BRT<:Ring}
  R = ambient_ring(X)
  R === ambient_ring(Y) || error("schemes can not be compared")
  return Spec(R, modulus(quotient_ring(OO(Y))), inverted_set(OO(X))*inverted_set(OO(Y)))
end


function Base.intersect(
    Y::AbsSpec{BRT, <:Ring},
    X::AbsSpec{BRT, <:MPolyLocalizedRing}
  ) where {BRT<:Ring}
  return intersect(X, Y)
end


### For Specs of MPolyQuoLocalizedRings
function Base.intersect(
    X::AbsSpec{BRT, <:MPolyQuoLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing}
  ) where {BRT<:Ring}
  R = ambient_ring(X)
  R === ambient_ring(Y) || error("schemes can not be compared")
#  Q, _ = quo(R, modulus(quotient_ring(OO(X))) + modulus(quotient_ring(OO(Y))))
  return Spec(R, modulus(quotient_ring(OO(X))) + modulus(quotient_ring(OO(Y))), 
              inverted_set(OO(X)) * inverted_set(OO(Y)))
end



########################################################
# (8) Closure of affine scheme in another affine scheme
########################################################

#TODO: Add more cross-type methods as needed.
@Markdown.doc """
    closure(X::AbsSpec, Y::AbsSpec) 

Returns the closure of ``X`` in ``Y``.

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

julia> H = subscheme(X,ideal(R,[x1]))
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1)

julia> closure(H, X)
Spec of Localization of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1) at the multiplicative set powers of fmpq_mpoly[1]
```
"""
function closure(X::AbsSpec, Y::AbsSpec) 
  return closure(standard_spec(X), standard_spec(Y))
end

function closure(
    X::Spec{BRT, RT}, 
    Y::Spec{BRT, RT}
  ) where {BRT, RT<:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                          <:MPolyPowersOfElement}}
  issubset(X, Y) || error("the first argument is not a subset of the second")
  is_closed_embedding(X, Y) && return X
  W, _ = Localization(inverted_set(OO(X))*inverted_set(OO(Y)))
  I = ideal(W, W.(gens(modulus(OO(X)))))
  Isat = saturated_ideal(I)
  R = ambient_ring(Y)
  return Spec(MPolyQuoLocalizedRing(R, Isat, inverted_set(OO(Y))))
end

