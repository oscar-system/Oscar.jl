

########################################################
# (1) Generic constructors
########################################################

@doc raw"""
    affine_scheme -> AffineScheme
    affine_scheme(::Ring)
    affine_scheme(::Ring, ::Ideal)
    affine_scheme(::Ideal)

Return the affine scheme defined by the input.
"""
affine_scheme


spec(I::Ideal) = affine_scheme(base_ring(I),I)
affine_scheme(I::Ideal) = affine_scheme(base_ring(I),I)

@doc raw"""
    spec(R::Ring) -> AffineScheme

Return the spectrum of the given ring `R` as an affine scheme.
"""
spec(R::Ring) = AffineScheme(R)
affine_scheme(R::Ring) = AffineScheme(R)

spec(kk::Ring, R::Ring) = AffineScheme(kk, R)
affine_scheme(kk::Ring, R::Ring) = AffineScheme(kk, R)

@doc raw"""
    spec(R::MPolyRing, I::MPolyIdeal)

Construct the affine scheme of the ideal ``I`` in the ring ``R``.
This is the spectrum of the quotient ring ``R/I``.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> I = ideal(R, [x]);

julia> spec(R, I)
Spectrum
  of quotient
    of multivariate polynomial ring in 2 variables x, y
      over rational field
    by ideal (x)
```
"""
spec(R::MPolyRing, I::MPolyIdeal) = AffineScheme(quo(R, I)[1])
affine_scheme(R::MPolyRing, I::MPolyIdeal) = AffineScheme(quo(R, I)[1])


@doc raw"""
    spec(R::MPolyRing, U::AbsMPolyMultSet)

Given a polynomial ring ``R``, we can localize that polynomial
ring at a multiplicatively closed subset ``U`` of ``R``. The spectrum
of the localized ring $U^{-1} R$ is computed by this method.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> I = ideal(R, [x]);

julia> U = complement_of_prime_ideal(I);

julia> spec(R, U)
Spectrum
  of localization
    of multivariate polynomial ring in 2 variables x, y
      over rational field
    at complement of prime ideal (x)
```
"""
spec(R::MPolyRing, U::AbsMPolyMultSet) = AffineScheme(localization(R, U)[1])
affine_scheme(R::MPolyRing, U::AbsMPolyMultSet) = AffineScheme(localization(R, U)[1])


@doc raw"""
    spec(R::MPolyRing, I::MPolyIdeal, U::AbsMPolyMultSet)

We allow to combine quotients and localizations at the same time.
That is, consider a polynomial ring ``R``, an ideal ``I`` of ``R`` and
a multiplicatively closed subset ``U`` of ``R``. The spectrum of the
localized ring $U^{-1} (R/I)$ is computed by this method.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> I = ideal(R, [x]);

julia> U = complement_of_prime_ideal(ideal(R, [y]));

julia> spec(R, I, U)
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 2 variables x, y
        over rational field
      by ideal (x)
    at complement of prime ideal (y)
```
"""
spec(R::MPolyRing, I::MPolyIdeal, U::AbsMPolyMultSet) = AffineScheme(MPolyQuoLocRing(R, I, U))
affine_scheme(R::MPolyRing, I::MPolyIdeal, U::AbsMPolyMultSet) = AffineScheme(MPolyQuoLocRing(R, I, U))



########################################################
# (2) Copy constructors
########################################################
#TODO: Do we need this? It is quite unusual.

@doc raw"""
    spec(X::AffineScheme)

For convenience, an affine spectrum can be passed to `AffineScheme`
to create a new spectrum. This can be particularly useful
when in need to copy an affine spectrum.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> I = ideal(R, [x]);

julia> X = spec(R, I)
Spectrum
  of quotient
    of multivariate polynomial ring in 2 variables x, y
      over rational field
    by ideal (x)

julia> Y = spec(X)
Spectrum
  of quotient
    of multivariate polynomial ring in 2 variables x, y
      over rational field
    by ideal (x)
```
"""
spec(X::AffineScheme) = AffineScheme(OO(X))
affine_scheme(X::AffineScheme) = AffineScheme(OO(X))

Base.deepcopy_internal(X::AffineScheme, dict::IdDict) = AffineScheme(deepcopy_internal(OO(X), dict))



########################################################
# (3) Affine n-dimensional space
########################################################

@doc raw"""
    affine_space(kk::BRT, n::Int; variable_name::VarName="x#") where {BRT<:Ring}

The ``n``-dimensional affine space over a ring ``kk`` is created
by this method. By default, the variable names are chosen as `x1`, `x2`
and so on. This choice can be overwritten with a third optional argument.

# Examples
```jldoctest
julia> affine_space(QQ, 5)
Affine space of dimension 5
  over rational field
with coordinates [x1, x2, x3, x4, x5]

julia> affine_space(QQ,5,variable_name="y#")
Affine space of dimension 5
  over rational field
with coordinates [y1, y2, y3, y4, y5]
```
"""
function affine_space(kk::BRT, n::Int; variable_name::VarName="x#") where {BRT<:Ring}
  R, _ = polynomial_ring(kk, variable_name => 1:n; cached=false)
  return spec(R)
end


@doc raw"""
    affine_space(kk::BRT, var_names::AbstractVector{<:VarName}) where {BRT<:Ring}

Create the ``n``-dimensional affine space over a ring ``kk``,
but allows more flexibility in the choice of variable names.
The following example demonstrates this.

# Examples
```jldoctest
julia> affine_space(QQ, [:x, :y, :z])
Affine space of dimension 3
  over rational field
with coordinates [x, y, z]

julia> affine_space(QQ, ['x', 'y', 'z'])
Affine space of dimension 3
  over rational field
with coordinates [x, y, z]

julia> affine_space(QQ, ["x", "y", "z"])
Affine space of dimension 3
  over rational field
with coordinates [x, y, z]
```
"""
function affine_space(kk::BRT, var_names::AbstractVector{<:VarName}) where {BRT<:Ring}
  R, _ = polynomial_ring(kk, var_names; cached=false)
  return spec(R)
end

function affine_space(kk::BRT, n::Int; variable_name::VarName="x#") where {BRT<:Field}
  R, _ = polynomial_ring(kk, variable_name => 1:n; cached=false)
  return variety(spec(R), check=false)
end

function affine_space(kk::BRT, var_names::AbstractVector{<:VarName}) where {BRT<:Field}
  R, _ = polynomial_ring(kk, var_names; cached=false)
  return variety(spec(R), check=false)
end

########################################################
# (4) StdAffineScheme (needed?)
# Calling this in a constructor should be avoided
# since we want to support all types of affine schemes.
# But StdAffineScheme can be useful when implementing some
# comparison or properties.
########################################################

@doc raw"""
    standard_spec(X::AbsAffineScheme)

For an affine scheme with coordinate ring of type `MPolyRing`,
`MPolyQuoRing`, or `MPolyLocRing`, return the canonical
transform to an `AffineScheme` of an `MPolyQuoLocRing`.

# Examples
```jldoctest
julia> Oscar.standard_spec(affine_space(QQ,5))
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 5 variables x1, x2, x3, x4, x5
        over rational field
      by ideal (0)
    at products of (1)

julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> I = ideal(R, [x]);

julia> X = spec(R, I)
Spectrum
  of quotient
    of multivariate polynomial ring in 2 variables x, y
      over rational field
    by ideal (x)

julia> Oscar.standard_spec(X)
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 2 variables x, y
        over rational field
      by ideal (x)
    at products of (1)

julia> I = ideal(R, [x]);

julia> U = complement_of_prime_ideal(I);

julia> X = spec(R, U)
Spectrum
  of localization
    of multivariate polynomial ring in 2 variables x, y
      over rational field
    at complement of prime ideal (x)

julia> Oscar.standard_spec(X)
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 2 variables x, y
        over rational field
      by ideal (0)
    at complement of prime ideal (x)
```
"""
function standard_spec(X::AbsAffineScheme)
  error("not implemented for input of type $(typeof(X))")
end

standard_spec(X::AbsAffineScheme{<:Any, <:MPolyRing}) = spec(MPolyQuoLocRing(OO(X), ideal(OO(X), [zero(OO(X))]), units_of(OO(X))))


# documented above
function standard_spec(X::AbsAffineScheme{<:Any, <:MPolyQuoRing})
  A = OO(X)
  R = base_ring(A)
  return AffineScheme(MPolyQuoLocRing(R, modulus(A), units_of(R)))
end


# documented above
standard_spec(X::AbsAffineScheme{<:Any, <:MPolyLocRing}) = AffineScheme(MPolyQuoLocRing(ambient_coordinate_ring(X), ideal(ambient_coordinate_ring(X), [zero(ambient_coordinate_ring(X))]), inverted_set(OO(X))))

#documented above
standard_spec(X::AbsAffineScheme{<:Any, <:MPolyQuoLocRing}) = AffineScheme(OO(X))



########################################################
# (5) Closed subschemes
########################################################

@doc raw"""
    subscheme(X::AbsAffineScheme, f::Vector{<:RingElem})

For an affine spectrum ``X`` and elements ``f_1``, ``f_2``,
etc. of the coordinate ring of ``X``, this method computes
the subscheme ``V(f_1, f_2, \dots)`` of ``X``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates [x1, x2, x3]

julia> R = OO(X)
Multivariate polynomial ring in 3 variables x1, x2, x3
  over rational field

julia> (x1,x2,x3) = gens(R)
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> subscheme(X,x1)
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal (x1)

julia> subscheme(X,[x1,x2])
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal (x1, x2)
```
"""
subscheme(X::AbsAffineScheme, f::Vector{<:RingElem}) = subscheme(X, ideal(OO(X), f))
function subscheme(X::AffineScheme, f::Vector{<:RingElem})
  all(x->(parent(x) == OO(X)), f) || return subscheme(X, OO(X).(f))
  return subscheme(X, ideal(OO(X), f))
end
subscheme(X::AbsAffineScheme, f::RingElem) = subscheme(X, ideal(OO(X), [f]))


@doc raw"""
    subscheme(X::AbsAffineScheme, I::Ideal)

For a scheme ``X = Spec(R)`` and an ideal ``I ⊂ 𝒪(X)``,
return the closed subscheme defined by ``I``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates [x1, x2, x3]

julia> R = OO(X)
Multivariate polynomial ring in 3 variables x1, x2, x3
  over rational field

julia> (x1,x2,x3) = gens(R)
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> subscheme(X,ideal(R,[x1*x2]))
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal (x1*x2)
```
"""
function subscheme(X::AbsAffineScheme, I::Ideal)
  base_ring(I) == OO(X) || return subscheme(X, ideal(OO(X), OO(X).(gens(I)))) # this will throw if coercion is not possible
  Y =  spec(quo(OO(X), I)[1])
  set_attribute!(Y, :ambient_space, ambient_space(X))
  return Y
end

function sub(X::AbsAffineScheme, a)
  (a isa RingElem && parent(a) === OO(X)) || return sub(X, OO(X)(a))
  return sub(X, ideal(OO(X), a))
end

function sub(X::AbsAffineScheme, a::Vector)
  all(a->a isa RingElem && parent(a) === OO(X), a) || return sub(X, OO(X).(a))
  return sub(X, ideal(OO(X), a))
end


########################################################
# (6) Open subschemes
########################################################

@doc raw"""
    hypersurface_complement(X::AbsAffineScheme, f::RingElem)

For a scheme ``X = Spec(R)`` and an element ``f ∈ R``,
return the open subscheme ``U = Spec(R[f⁻¹]) = X ∖ V(f)``
defined by the complement of the vanishing locus of ``f``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates [x1, x2, x3]

julia> R = OO(X)
Multivariate polynomial ring in 3 variables x1, x2, x3
  over rational field

julia> (x1, x2, x3) = gens(R)
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> hypersurface_complement(X, x1)
Spectrum
  of localization
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    at products of (x1)
```
"""
function hypersurface_complement(X::AbsAffineScheme, f::RingElem)
  return hypersurface_complement(underlying_scheme(X), f)::AbsAffineScheme
end

function hypersurface_complement(X::AffineSchemeType, f::RingElem) where {AffineSchemeType<:AbsAffineScheme{<:Any, <:MPolyQuoLocRing}}
  parent(f) == OO(X) || return hypersurface_complement(X, OO(X)(f))
  h = lifted_numerator(f)
  U = MPolyPowersOfElement(h)
  #simplify!(U)
  simplify_light!(U)
  W, _ = localization(OO(X), U)
  Y = spec(W)
  set_attribute!(Y, :ambient_space, ambient_space(X))
  return Y
end

function hypersurface_complement(X::AffineSchemeType, f::RingElem) where {AffineSchemeType<:AbsAffineScheme{<:Any, <:MPolyLocRing}}
  parent(f) == OO(X) || return hypersurface_complement(X, OO(X)(f))
  h = numerator(f)
  U = MPolyPowersOfElement(h)
  #simplify!(U)
  simplify_light!(U)
  W, _ = localization(OO(X), U)
  Y = spec(W)
  set_attribute!(Y, :ambient_space, ambient_space(X))
  return Y
end

function hypersurface_complement(X::AffineSchemeType, f::RingElem) where {AffineSchemeType<:AbsAffineScheme{<:Any, <:MPolyRing}}
  parent(f) == OO(X) || return hypersurface_complement(X, OO(X)(f))
  U = MPolyPowersOfElement(f)
  #simplify!(U)
  simplify_light!(U)
  W, _ = localization(OO(X), U)
  Y = spec(W)
  set_attribute!(Y, :ambient_space, ambient_space(X))
  return Y
end

function hypersurface_complement(X::AffineSchemeType, f::RingElem) where {AffineSchemeType<:AbsAffineScheme{<:Any, <:MPolyQuoRing}}
  parent(f) == OO(X) || return hypersurface_complement(X, OO(X)(f))
  U = MPolyPowersOfElement(lift(f))
  #simplify!(U)
  simplify_light!(U)
  W, _ = localization(OO(X), U)
  Y = spec(W)
  set_attribute!(Y, :ambient_space, ambient_space(X))
  return Y
end


@doc raw"""
    hypersurface_complement(X::AbsAffineScheme, f::Vector{<:RingElem})

For a scheme ``X = Spec(R)`` and elements ``f₁, f₂, ... ∈ R``,
return the open subscheme ``U = Spec(R[f₁⁻¹,f₂⁻¹, ...]) = X ∖ V(f₁⋅f₂⋅…)``
defined by the complement of the vanishing locus of the product ``f₁⋅f₂⋅…``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates [x1, x2, x3]

julia> R = OO(X)
Multivariate polynomial ring in 3 variables x1, x2, x3
  over rational field

julia> (x1,x2,x3) = gens(R)
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> hypersurface_complement(X,[x1,x2])
Spectrum
  of localization
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    at products of (x1, x2)
```
"""
function hypersurface_complement(X::AbsAffineScheme, f::Vector{<:RingElem})
  return hypersurface_complement(underlying_scheme(X), f)::AbsAffineScheme
end

function hypersurface_complement(X::AffineSchemeType, f::Vector{<:RingElem}) where {AffineSchemeType<:AbsAffineScheme{<:Any, <:MPolyQuoLocRing}}
  all(x->(parent(x) == OO(X)), f) || return hypersurface_complement(X, OO(X).(f))
  h = lifted_numerator.(f)
  U = MPolyPowersOfElement(ambient_coordinate_ring(X), h)
  #simplify!(U)
  simplify_light!(U)
  W, _ = localization(OO(X), U)
  Y = spec(W)
  set_attribute!(Y, :ambient_space, ambient_space(X))
  return Y
end

function hypersurface_complement(X::AffineSchemeType, f::Vector{<:RingElem}) where {AffineSchemeType<:AbsAffineScheme{<:Any, <:MPolyLocRing}}
  all(x->(parent(x) == OO(X)), f) || return hypersurface_complement(X, OO(X).(f))
  h = numerator.(f)
  U = MPolyPowersOfElement(ambient_coordinate_ring(X), h)
  #simplify!(U)
  simplify_light!(U)
  W, _ = localization(OO(X), U)
  Y = spec(W)
  set_attribute!(Y, :ambient_space, ambient_space(X))
  return Y
end

function hypersurface_complement(X::AffineSchemeType, f::Vector{<:RingElem}) where {AffineSchemeType<:AbsAffineScheme{<:Any, <:MPolyRing}}
  all(x->(parent(x) == OO(X)), f) || return hypersurface_complement(X, OO(X).(f))
  U = MPolyPowersOfElement(ambient_coordinate_ring(X), f)
  #simplify!(U)
  simplify_light!(U)
  W, _ = localization(OO(X), U)
  Y = spec(W)
  set_attribute!(Y, :ambient_space, ambient_space(X))
  return Y
end

function hypersurface_complement(X::AffineSchemeType, f::Vector{<:RingElem}) where {AffineSchemeType<:AbsAffineScheme{<:Any, <:MPolyQuoRing}}
  all(x->(parent(x) == OO(X)), f) || return hypersurface_complement(X, OO(X).(f))
  U = MPolyPowersOfElement(ambient_coordinate_ring(X), lift.(f))
  #simplify!(U)
  simplify_light!(U)
  W, _ = localization(OO(X), U)
  Y = spec(W)
  set_attribute!(Y, :ambient_space, ambient_space(X))
  return Y
end



########################################################
# (7) Intersections of (affine) schemes
########################################################


# 7.1 Intersection methods involving an empty scheme
Base.intersect(E::EmptyScheme{BRT}, X::Scheme{BRT}) where {BRT<:Ring} = E
Base.intersect(X::Scheme{BRT}, E::EmptyScheme{BRT}) where {BRT<:Ring} = E
Base.intersect(X::EmptyScheme{BRT}, E::EmptyScheme{BRT}) where {BRT<:Ring} = E


# 7.2 Intersection methods not involving empty schemes

### For AffineSchemes of MPolyRings
# TODO  intersect X,Y for X<Y should return a copy of X with === ambient_coordinate_rings
# spec(X) does not apply for instance to principal open subsets hence a change
# is necessary
@doc raw"""
    Base.intersect(X::AbsAffineScheme, Y::AbsAffineScheme)

This method computes the intersection to two affine
schemes that reside in the same ambient affine space.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates [x1, x2, x3]

julia> R = OO(X)
Multivariate polynomial ring in 3 variables x1, x2, x3
  over rational field

julia> (x1,x2,x3) = gens(R)
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> Y1 = subscheme(X,[x1])
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal (x1)

julia> Y2 = subscheme(X,[x2])
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal (x2)

julia> intersect(Y1, Y2)
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal (x1, x2)
```
"""
function Base.intersect(X::AbsAffineScheme{BRT, <:Ring}, Y::AbsAffineScheme{BRT, <:Ring}) where {BRT<:Ring}
  error("method not implemented for arguments of type $(typeof(X)) and $(typeof(Y))")
end

function Base.intersect(
    X::AbsAffineScheme{BRT, <:MPolyRing},
    Y::AbsAffineScheme{BRT, <:MPolyRing}
  ) where {BRT<:Ring}
  R = OO(X)
  R == OO(Y) || error("schemes can not be compared")
  return spec(X)
end


function Base.intersect(
    X::AbsAffineScheme{BRT, <:MPolyRing},
    Y::AbsAffineScheme{BRT, <:MPolyQuoRing}
  ) where {BRT<:Ring}
  R = OO(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  return spec(Y)
end


function Base.intersect(
    X::AbsAffineScheme{BRT, <:MPolyRing},
    Y::AbsAffineScheme{BRT, <:MPolyLocRing}
  ) where {BRT<:Ring}
  R = OO(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  return spec(Y)
end


function Base.intersect(
    X::AbsAffineScheme{BRT, <:MPolyRing},
    Y::AbsAffineScheme{BRT, <:MPolyQuoLocRing}
  ) where {BRT<:Ring}
  R = OO(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  return spec(Y)
end


function Base.intersect(
    Y::AbsAffineScheme{BRT, <:Ring},
    X::AbsAffineScheme{BRT, <:MPolyRing}
  ) where {BRT<:Ring}
  return intersect(X, Y)
end


### For AffineSchemes of MPolyQuos
function Base.intersect(
    X::AbsAffineScheme{BRT, <:MPolyQuoRing},
    Y::AbsAffineScheme{BRT, <:MPolyQuoRing}
  ) where {BRT<:Ring}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  return spec(quo(R, modulus(OO(X)) + modulus(OO(Y)))[1])
end


function Base.intersect(
    X::AbsAffineScheme{BRT, <:MPolyQuoRing},
    Y::AbsAffineScheme{BRT, <:MPolyLocRing}
  ) where {BRT<:Ring}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  return spec(quo(OO(Y), OO(Y)(modulus(OO(X))))[1])
end


function Base.intersect(
    X::AbsAffineScheme{BRT, <:MPolyQuoRing},
    Y::AbsAffineScheme{BRT, <:MPolyQuoLocRing}
  ) where {BRT<:Ring}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  return spec(quo(OO(Y), OO(Y)(modulus(OO(X))))[1])
end


function Base.intersect(
    Y::AbsAffineScheme{BRT, <:Ring},
    X::AbsAffineScheme{BRT, <:MPolyQuoRing}
  ) where {BRT<:Ring}
  return intersect(X, Y)
end


### For AffineSchemes of MPolyLocalizedRings
function Base.intersect(
    X::AbsAffineScheme{BRT, <:MPolyLocRing},
    Y::AbsAffineScheme{BRT, <:MPolyLocRing}
  ) where {BRT<:Ring}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  return spec(localization(R, inverted_set(OO(X)) * inverted_set(OO(Y)))[1])
end


function Base.intersect(
    X::AbsAffineScheme{BRT, <:MPolyLocRing},
    Y::AbsAffineScheme{BRT, <:MPolyQuoLocRing}
  ) where {BRT<:Ring}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  return spec(R, modulus(underlying_quotient(OO(Y))), inverted_set(OO(X))*inverted_set(OO(Y)))
end


function Base.intersect(
    Y::AbsAffineScheme{BRT, <:Ring},
    X::AbsAffineScheme{BRT, <:MPolyLocRing}
  ) where {BRT<:Ring}
  return intersect(X, Y)
end


### For AffineSchemes of MPolyQuoLocalizedRings
function Base.intersect(
    X::AbsAffineScheme{BRT, <:MPolyQuoLocRing},
    Y::AbsAffineScheme{BRT, <:MPolyQuoLocRing}
  ) where {BRT<:Ring}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
#  Q, _ = quo(R, modulus(underlying_quotient(OO(X))) + modulus(underlying_quotient(OO(Y))))
  return spec(R, modulus(underlying_quotient(OO(X))) + modulus(underlying_quotient(OO(Y))),
              inverted_set(OO(X)) * inverted_set(OO(Y)))
end



########################################################
# (8) Closure of affine scheme in another affine scheme
########################################################

#TODO: Add more cross-type methods as needed.
@doc raw"""
    closure(X::AbsAffineScheme, Y::AbsAffineScheme)

Return the closure of ``X`` in ``Y``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates [x1, x2, x3]

julia> R = OO(X)
Multivariate polynomial ring in 3 variables x1, x2, x3
  over rational field

julia> (x1,x2,x3) = gens(R)
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> H = subscheme(X,ideal(R,[x1]))
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal (x1)

julia> closure(H, X)
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal (x1)
```
"""
function closure(X::AbsAffineScheme, Y::AbsAffineScheme, check= true)
  error("not implemented")
end

function closure(
    X::AbsAffineScheme{BRT, <:Union{MPolyQuoRing,MPolyRing}},
    Y::AbsAffineScheme{BRT, <:MPolyAnyRing};
    check::Bool=true
  ) where {BRT}
  @check is_subscheme(X, Y) "the first argument is not a subset of the second"
  return X
end

function closure(
    X::AbsAffineScheme{BRT, <:MPolyLocRing},
    Y::AbsAffineScheme{BRT, <:MPolyAnyRing};
    check::Bool=true
  ) where {BRT}
  @check is_subscheme(X, Y) "the first argument is not a subset of the second"
  return Y
end

function closure(
    X::AbsAffineScheme{BRT, <:MPolyLocRing},
    Y::AbsAffineScheme{BRT, <:MPolyLocRing};
    check::Bool=true
  ) where {BRT}
  @check is_subscheme(X, Y) "the first argument is not a subset of the second"
  return Y
end


function closure(
    X::AbsAffineScheme{BRT, <:MPolyQuoLocRing},
    Y::AbsAffineScheme{BRT, <:Union{MPolyRing,MPolyQuoRing}};
    check::Bool=true
  ) where {BRT}
  @check is_subscheme(X, Y) "the first argument is not a subset of the second"
  I = saturated_ideal(defining_ideal(X))
  return spec(base_ring(I),I)
end

function closure(
    X::AbsAffineScheme{BRT, <:MPolyQuoLocRing},
    Y::AbsAffineScheme{BRT, <:MPolyLocRing};
    check::Bool=true
  ) where {BRT}
  @check is_subscheme(X, Y) "the first argument is not a subset of the second"
  I = saturated_ideal(defining_ideal(X))
  R = base_ring(I)
  return spec(MPolyQuoLocRing(R, I, inverted_set(Y)))
end

function closure(
    X::AbsAffineScheme{BRT, RT},
    Y::AbsAffineScheme{BRT, RT};
    check::Bool=true
  ) where {BRT, RT<:MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any,
                                    <:MPolyPowersOfElement}}
  @check is_subscheme(X, Y) "the first argument is not a subset of the second"
  #is_closed_embedding(X, Y) && return X
  W, _ = localization(inverted_set(OO(X))*inverted_set(OO(Y)))
  I = ideal(W, W.(gens(modulus(OO(X)))))
  Isat = saturated_ideal(I)
  R = ambient_coordinate_ring(Y)
  return spec(MPolyQuoLocRing(R, Isat, inverted_set(OO(Y))))
end

@doc raw"""
    closure(X::AbsAffineScheme) -> AbsAffineScheme

Return the closure of `X` in its ambient affine space.
"""
closure(X::AbsAffineScheme) = closure(X, ambient_space(X), check= true)



######################################################################
# Unions
######################################################################

function union(X::AbsAffineScheme{BRT,RT}, Y::AbsAffineScheme{BRT,RT}) where {BRT, RT<:MPolyQuoRing}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  IX = modulus(OO(X))
  IY = modulus(OO(Y))
  return spec(R, intersect(IX, IY))
end

function union(X::AbsAffineScheme{BRT,<:MPolyQuoRing}, Y::AbsAffineScheme{BRT,<:MPolyRing}) where {BRT}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  return X
end

function union(X::AbsAffineScheme{BRT,<:MPolyRing}, Y::AbsAffineScheme{BRT,<:MPolyQuoRing}) where {BRT}
  return union(Y,X)
end
