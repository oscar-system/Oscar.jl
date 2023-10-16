

########################################################
# (1) Generic constructors
########################################################

@doc raw"""
    Spec(R::MPolyRing, I::MPolyIdeal)

Constructs the affine scheme of the ideal ``I`` in the ring ``R``.
This is the spectrum of the quotient ring ``R/I``.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> I = ideal(R, [x]);

julia> Spec(R, I)
Spectrum
  of quotient
    of multivariate polynomial ring in 2 variables x, y
      over rational field
    by ideal(x)
```
"""
Spec(R::MPolyRing, I::MPolyIdeal) = Spec(quo(R, I)[1])


@doc raw"""
    Spec(R::MPolyRing, U::AbsMPolyMultSet)

Given a polynomial ring ``R``, we can localize that polynomial
ring at a multiplicatively closed subset ``U`` of ``R``. The spectrum
of the localized ring $U^{-1} R$ is computed by this method.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> I = ideal(R, [x]);

julia> U = complement_of_prime_ideal(I);

julia> Spec(R, U)
Spectrum
  of localization
    of multivariate polynomial ring in 2 variables x, y
      over rational field
    at complement of prime ideal(x)
```
"""
Spec(R::MPolyRing, U::AbsMPolyMultSet) = Spec(Localization(R, U)[1])


@doc raw"""
    Spec(R::MPolyRing, I::MPolyIdeal, U::AbsMPolyMultSet)

We allow to combine quotients and localizations at the same time.
That is, consider a polynomial ring ``R``, an ideal ``I`` of ``R`` and
a multiplicatively closed subset ``U`` of ``R``. The spectrum of the
localized ring $U^{-1} (R/I)$ is computed by this method.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> I = ideal(R, [x]);

julia> U = complement_of_prime_ideal(ideal(R, [y]));

julia> Spec(R, I, U)
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 2 variables x, y
        over rational field
      by ideal(x)
    at complement of prime ideal(y)
```
"""
Spec(R::MPolyRing, I::MPolyIdeal, U::AbsMPolyMultSet) = Spec(MPolyQuoLocRing(R, I, U))



########################################################
# (2) Copy constructors
########################################################
#TODO: Do we need this? It is quite unusual.

@doc raw"""
    Spec(X::Spec)

For convenience, an affine spectrum can be passed to `Spec`
to create a new spectrum. This can be particularly useful
when in need to copy an affine spectrum.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> I = ideal(R, [x]);

julia> X = Spec(R, I)
Spectrum
  of quotient
    of multivariate polynomial ring in 2 variables x, y
      over rational field
    by ideal(x)

julia> Y = Spec(X)
Spectrum
  of quotient
    of multivariate polynomial ring in 2 variables x, y
      over rational field
    by ideal(x)
```
"""
Spec(X::Spec) = Spec(OO(X))

Base.deepcopy_internal(X::Spec, dict::IdDict) = Spec(deepcopy_internal(OO(X), dict))



########################################################
# (3) Affine n-dimensional space
########################################################

@doc raw"""
    affine_space(kk::BRT, n::Int; variable_name="x") where {BRT<:Ring}

The ``n``-dimensional affine space over a ring ``kk`` is created
by this method. By default, the variable names are chosen as $x_1$, $x_2$
and so on. This choice can be overwritten with a third optional argument.

# Examples
```jldoctest
julia> affine_space(QQ, 5)
Affine space of dimension 5
  over rational field
with coordinates [x1, x2, x3, x4, x5]

julia> affine_space(QQ,5,variable_name="y")
Affine space of dimension 5
  over rational field
with coordinates [y1, y2, y3, y4, y5]
```
"""
function affine_space(kk::BRT, n::Int; variable_name="x") where {BRT<:Ring}
  R, _ = polynomial_ring(kk, [variable_name * "$i" for i in 1:n])
  return Spec(R)
end


@doc raw"""
    affine_space(kk::BRT, var_symbols::Vector{Symbol}) where {BRT<:Ring}

Creates the ``n``-dimensional affine space over a ring ``kk``,
but allows more flexibility in the choice of variable names.
The following example demonstrates this.

# Examples
```jldoctest
julia> affine_space(QQ,[:y1,:z2,:a])
Affine space of dimension 3
  over rational field
with coordinates [y1, z2, a]
```
"""
function affine_space(kk::BRT, var_symbols::Vector{Symbol}) where {BRT<:Ring}
  R, _ = polynomial_ring(kk, var_symbols)
  return variety(Spec(R), check=false)
end

function affine_space(kk::BRT, n::Int; variable_name="x") where {BRT<:Field}
  R, _ = polynomial_ring(kk, [variable_name * "$i" for i in 1:n])
  return variety(Spec(R), check=false)
end

function affine_space(kk::BRT, var_symbols::Vector{Symbol}) where {BRT<:Field}
  R, _ = polynomial_ring(kk, var_symbols)
  return variety(Spec(R), check=false)
end

########################################################
# (4) StdSpec (needed?)
# Calling this in a constructor should be avoided
# since we want to support all types of affine schemes.
# But StdSpec can be useful when implementing some
# comparison or properties.
########################################################

@doc raw"""
    standard_spec(X::AbsSpec)

For an affine spectrum with coordinate ring of type `MPolyRing`, 
`MPolyQuoRing`, or `MPolyLocRing`, return the canonical
transform to a `Spec` of an `MPolyQuoLocRing`. 

# Examples
```jldoctest
julia> standard_spec(affine_space(QQ,5))
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 5 variables x1, x2, x3, x4, x5
        over rational field
      by ideal(0)
    at products of (1)

julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> I = ideal(R, [x]);

julia> X = Spec(R, I)
Spectrum
  of quotient
    of multivariate polynomial ring in 2 variables x, y
      over rational field
    by ideal(x)

julia> standard_spec(X)
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 2 variables x, y
        over rational field
      by ideal(x)
    at products of (1)

julia> I = ideal(R, [x]);

julia> U = complement_of_prime_ideal(I);

julia> X = Spec(R, U)
Spectrum
  of localization
    of multivariate polynomial ring in 2 variables x, y
      over rational field
    at complement of prime ideal(x)

julia> standard_spec(X)
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 2 variables x, y
        over rational field
      by ideal(0)
    at complement of prime ideal(x)
```
"""
function standard_spec(X::AbsSpec)
  error("not implemented for input of type $(typeof(X))")
end

standard_spec(X::AbsSpec{<:Any, <:MPolyRing}) = Spec(MPolyQuoLocRing(OO(X), ideal(OO(X), [zero(OO(X))]), units_of(OO(X))))


# documented above
function standard_spec(X::AbsSpec{<:Any, <:MPolyQuoRing})
  A = OO(X)
  R = base_ring(A)
  return Spec(MPolyQuoLocRing(R, modulus(A), units_of(R)))
end


# documented above
standard_spec(X::AbsSpec{<:Any, <:MPolyLocRing}) = Spec(MPolyQuoLocRing(ambient_coordinate_ring(X), ideal(ambient_coordinate_ring(X), [zero(ambient_coordinate_ring(X))]), inverted_set(OO(X))))

#documented above
standard_spec(X::AbsSpec{<:Any, <:MPolyQuoLocRing}) = Spec(OO(X))



########################################################
# (5) Closed subschemes
########################################################

@doc raw"""
    subscheme(X::AbsSpec, f::Vector{<:RingElem})

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
    by ideal(x1)

julia> subscheme(X,[x1,x2])
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal(x1, x2)
```
"""
subscheme(X::AbsSpec, f::Vector{<:RingElem}) = subscheme(X, ideal(OO(X), f))
function subscheme(X::Spec, f::Vector{<:RingElem})
  all(x->(parent(x) == OO(X)), f) || return subscheme(X, OO(X).(f))
  return subscheme(X, ideal(OO(X), f))
end
subscheme(X::AbsSpec, f::RingElem) = subscheme(X, ideal(OO(X), [f]))


@doc raw"""
    subscheme(X::AbsSpec, I::Ideal)

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
    by ideal(x1*x2)
```
"""
function subscheme(X::AbsSpec, I::Ideal)
  base_ring(I) == OO(X) || return subscheme(X, ideal(OO(X), OO(X).(gens(I)))) # this will throw if coercion is not possible
  Y =  Spec(quo(OO(X), I)[1])
  set_attribute!(Y, :ambient_space, ambient_space(X))
  return Y
end



########################################################
# (6) Open subschemes
########################################################

@doc raw"""
    hypersurface_complement(X::AbsSpec, f::RingElem)

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
function hypersurface_complement(X::AbsSpec, f::RingElem)
  return hypersurface_complement(underlying_scheme(X), f)::AbsSpec
end

function hypersurface_complement(X::SpecType, f::RingElem) where {SpecType<:AbsSpec{<:Any, <:MPolyQuoLocRing}}
  parent(f) == OO(X) || return hypersurface_complement(X, OO(X)(f))
  h = lifted_numerator(f)
  U = MPolyPowersOfElement(h)
  simplify!(U)
  W, _ = Localization(OO(X), U)
  Y = Spec(W)
  set_attribute!(Y, :ambient_space, ambient_space(X))
  return Y
end

function hypersurface_complement(X::SpecType, f::RingElem) where {SpecType<:AbsSpec{<:Any, <:MPolyLocRing}}
  parent(f) == OO(X) || return hypersurface_complement(X, OO(X)(f))
  h = numerator(f)
  U = MPolyPowersOfElement(h)
  simplify!(U)
  W, _ = Localization(OO(X), U)
  Y = Spec(W)
  set_attribute!(Y, :ambient_space, ambient_space(X))
  return Y
end

function hypersurface_complement(X::SpecType, f::RingElem) where {SpecType<:AbsSpec{<:Any, <:MPolyRing}}
  parent(f) == OO(X) || return hypersurface_complement(X, OO(X)(f))
  U = MPolyPowersOfElement(f)
  simplify!(U)
  W, _ = Localization(OO(X), U)
  Y = Spec(W)
  set_attribute!(Y, :ambient_space, ambient_space(X))
  return Y
end

function hypersurface_complement(X::SpecType, f::RingElem) where {SpecType<:AbsSpec{<:Any, <:MPolyQuoRing}}
  parent(f) == OO(X) || return hypersurface_complement(X, OO(X)(f))
  U = MPolyPowersOfElement(lift(f))
  simplify!(U)
  W, _ = Localization(OO(X), U)
  Y = Spec(W)
  set_attribute!(Y, :ambient_space, ambient_space(X))
  return Y
end


@doc raw"""
    hypersurface_complement(X::AbsSpec, f::Vector{<:RingElem})

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
    at products of (x1,x2)
```
"""
function hypersurface_complement(X::AbsSpec, f::Vector{<:RingElem})
  return hypersurface_complement(underlying_scheme(X), f)::AbsSpec
end

function hypersurface_complement(X::SpecType, f::Vector{<:RingElem}) where {SpecType<:AbsSpec{<:Any, <:MPolyQuoLocRing}}
  all(x->(parent(x) == OO(X)), f) || return hypersurface_complement(X, OO(X).(f))
  h = lifted_numerator.(f)
  U = MPolyPowersOfElement(ambient_coordinate_ring(X), h)
  simplify!(U)
  W, _ = Localization(OO(X), U)
  Y = Spec(W)
  set_attribute!(Y, :ambient_space, ambient_space(X))
  return Y
end

function hypersurface_complement(X::SpecType, f::Vector{<:RingElem}) where {SpecType<:AbsSpec{<:Any, <:MPolyLocRing}}
  all(x->(parent(x) == OO(X)), f) || return hypersurface_complement(X, OO(X).(f))
  h = numerator.(f)
  U = MPolyPowersOfElement(ambient_coordinate_ring(X), h)
  simplify!(U)
  W, _ = Localization(OO(X), U)
  Y = Spec(W)
  set_attribute!(Y, :ambient_space, ambient_space(X))
  return Y
end

function hypersurface_complement(X::SpecType, f::Vector{<:RingElem}) where {SpecType<:AbsSpec{<:Any, <:MPolyRing}}
  all(x->(parent(x) == OO(X)), f) || return hypersurface_complement(X, OO(X).(f))
  U = MPolyPowersOfElement(ambient_coordinate_ring(X), f)
  simplify!(U)
  W, _ = Localization(OO(X), U)
  Y = Spec(W)
  set_attribute!(Y, :ambient_space, ambient_space(X))
  return Y
end

function hypersurface_complement(X::SpecType, f::Vector{<:RingElem}) where {SpecType<:AbsSpec{<:Any, <:MPolyQuoRing}}
  all(x->(parent(x) == OO(X)), f) || return hypersurface_complement(X, OO(X).(f))
  U = MPolyPowersOfElement(ambient_coordinate_ring(X), lift.(f))
  simplify!(U)
  W, _ = Localization(OO(X), U)
  Y = Spec(W)
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

### For Specs of MPolyRings
# TODO  intersect X,Y for X<Y should return a copy of X with === ambient_coordinate_rings
# Spec(X) does not apply for instance to principal open subsets hence a change
# is necessary
@doc raw"""
    Base.intersect(X::AbsSpec, Y::AbsSpec)

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
    by ideal(x1)

julia> Y2 = subscheme(X,[x2])
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal(x2)

julia> intersect(Y1, Y2)
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal(x1, x2)
```
"""
function Base.intersect(X::AbsSpec{BRT, <:Ring}, Y::AbsSpec{BRT, <:Ring}) where {BRT<:Ring}
  error("method not implemented for arguments of type $(typeof(X)) and $(typeof(Y))")
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
    Y::AbsSpec{BRT, <:MPolyQuoRing}
  ) where {BRT<:Ring}
  R = OO(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  return Spec(Y)
end


function Base.intersect(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyLocRing}
  ) where {BRT<:Ring}
  R = OO(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  return Spec(Y)
end


function Base.intersect(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocRing}
  ) where {BRT<:Ring}
  R = OO(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
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
    X::AbsSpec{BRT, <:MPolyQuoRing},
    Y::AbsSpec{BRT, <:MPolyQuoRing}
  ) where {BRT<:Ring}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  return Spec(quo(R, modulus(OO(X)) + modulus(OO(Y)))[1])
end


function Base.intersect(
    X::AbsSpec{BRT, <:MPolyQuoRing},
    Y::AbsSpec{BRT, <:MPolyLocRing}
  ) where {BRT<:Ring}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  return Spec(quo(OO(Y), OO(Y)(modulus(OO(X))))[1])
end


function Base.intersect(
    X::AbsSpec{BRT, <:MPolyQuoRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocRing}
  ) where {BRT<:Ring}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  return Spec(quo(OO(Y), OO(Y)(modulus(OO(X))))[1])
end


function Base.intersect(
    Y::AbsSpec{BRT, <:Ring},
    X::AbsSpec{BRT, <:MPolyQuoRing}
  ) where {BRT<:Ring}
  return intersect(X, Y)
end


### For Specs of MPolyLocalizedRings
function Base.intersect(
    X::AbsSpec{BRT, <:MPolyLocRing},
    Y::AbsSpec{BRT, <:MPolyLocRing}
  ) where {BRT<:Ring}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  return Spec(Localization(R, inverted_set(OO(X)) * inverted_set(OO(Y)))[1])
end


function Base.intersect(
    X::AbsSpec{BRT, <:MPolyLocRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocRing}
  ) where {BRT<:Ring}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
  return Spec(R, modulus(underlying_quotient(OO(Y))), inverted_set(OO(X))*inverted_set(OO(Y)))
end


function Base.intersect(
    Y::AbsSpec{BRT, <:Ring},
    X::AbsSpec{BRT, <:MPolyLocRing}
  ) where {BRT<:Ring}
  return intersect(X, Y)
end


### For Specs of MPolyQuoLocalizedRings
function Base.intersect(
    X::AbsSpec{BRT, <:MPolyQuoLocRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocRing}
  ) where {BRT<:Ring}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || error("schemes can not be compared")
#  Q, _ = quo(R, modulus(underlying_quotient(OO(X))) + modulus(underlying_quotient(OO(Y))))
  return Spec(R, modulus(underlying_quotient(OO(X))) + modulus(underlying_quotient(OO(Y))), 
              inverted_set(OO(X)) * inverted_set(OO(Y)))
end



########################################################
# (8) Closure of affine scheme in another affine scheme
########################################################

#TODO: Add more cross-type methods as needed.
@doc raw"""
    closure(X::AbsSpec, Y::AbsSpec) 

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
    by ideal(x1)

julia> closure(H, X)
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal(x1)
```
"""
function closure(X::AbsSpec, Y::AbsSpec, check= true)
  error("not implemented")
end

function closure(
    X::AbsSpec{BRT, <:Union{MPolyQuoRing,MPolyRing}},
    Y::AbsSpec{BRT, <:MPolyAnyRing};
    check::Bool=true
  ) where {BRT}
  @check issubset(X, Y) "the first argument is not a subset of the second"
  return X
end

function closure(
    X::AbsSpec{BRT, <:MPolyLocRing},
    Y::AbsSpec{BRT, <:MPolyAnyRing};
    check::Bool=true
  ) where {BRT}
  @check issubset(X, Y) "the first argument is not a subset of the second"
  return Y
end

function closure(
    X::AbsSpec{BRT, <:MPolyLocRing},
    Y::AbsSpec{BRT, <:MPolyLocRing};
    check::Bool=true
  ) where {BRT}
  @check issubset(X, Y) "the first argument is not a subset of the second"
  return Y
end


function closure(
    X::AbsSpec{BRT, <:MPolyQuoLocRing},
    Y::AbsSpec{BRT, <:Union{MPolyRing,MPolyQuoRing}};
    check::Bool=true
  ) where {BRT}
  @check issubset(X, Y) "the first argument is not a subset of the second"
  I = ambient_closure_ideal(X)
  return Spec(base_ring(I),I)
end

function closure(
    X::AbsSpec{BRT, <:MPolyQuoLocRing},
    Y::AbsSpec{BRT, <:MPolyLocRing};
    check::Bool=true
  ) where {BRT}
  @check issubset(X, Y) "the first argument is not a subset of the second"
  I = ambient_closure_ideal(X)
  R = base_ring(I)
  return Spec(MPolyQuoLocRing(R, I, inverted_set(Y)))
end

function closure(
    X::AbsSpec{BRT, RT},
    Y::AbsSpec{BRT, RT};
    check::Bool=true
  ) where {BRT, RT<:MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any,
                                    <:MPolyPowersOfElement}}
  @check issubset(X, Y) "the first argument is not a subset of the second"
  #is_closed_embedding(X, Y) && return X
  W, _ = Localization(inverted_set(OO(X))*inverted_set(OO(Y)))
  I = ideal(W, W.(gens(modulus(OO(X)))))
  Isat = saturated_ideal(I)
  R = ambient_coordinate_ring(Y)
  return Spec(MPolyQuoLocRing(R, Isat, inverted_set(OO(Y))))
end

@doc raw"""
    closure(X::AbsSpec) -> AbsSpec

Return the closure of `X` in its ambient affine space.
"""
closure(X::AbsSpec) = closure(X, ambient_space(X), check= true)
