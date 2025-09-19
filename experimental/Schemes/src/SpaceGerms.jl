
import AbstractAlgebra: Ring

@doc raw"""
    AbsSpaceGerm{BaseRingType<:Ring,RingType<:Ring}

A space germ, i.e. a ringed space ``(X,O_{(X,x)})`` with local ring ``O_{(X,x)}`` of type `RingType` over a coefficient field ``k`` of type `BaseRingType`.
"""
abstract type AbsSpaceGerm{BaseRingType<:Ring, RingType<:Ring} <: AbsAffineScheme{BaseRingType, RingType} end


####################################################################
## two short hand definitions for internal use only
####################################################################
const GermAtClosedPoint = AffineScheme{<:Field,
                         <:AbsLocalizedRing{<:Ring, <:Any,
                                            <:MPolyComplementOfKPointIdeal}
                        }
const GermAtGeometricPoint = AffineScheme{<:Field,
                            <:AbsLocalizedRing{<:Ring, <:Any,
                                               <:MPolyComplementOfPrimeIdeal}
                           }

###################################################################
## start of the definition of space germ functionality
###################################################################

@doc raw"""
    SpaceGerm{BaseRingType, RingType, AffineSchemeType}
A space germ ``(X,O_{(X,x)}``, i.e. a ringed space with underlying scheme ``X`` of type AffineSchemeType and local ring ``O_{(X,x)}`` of type `RingType` over some base ring ``k`` of type `BaseRingType`.
"""
@attributes mutable struct SpaceGerm{
                BaseRingType<:Ring,
                RingType<:Ring,
                AffineSchemeType<:AffineScheme} <: AbsSpaceGerm{BaseRingType, RingType}
  X::AffineSchemeType

  function SpaceGerm(X::GermAtClosedPoint)
    return new{typeof(base_ring(X)), typeof(OO(X)), typeof(X)}(X)
  end

## the following one is currently unused..
  function SpaceGerm(X::GermAtGeometricPoint)
    return new{typeof(base_ring(X)), typeof(OO(X)), typeof(X)}(X)
  end
end

@doc raw"""
    HypersurfaceGerm{BaseRingType, RingType, AffineSchemeType}

A hypersurface germ ``(X,O_{(X,x)}``, i.e. a ringed space with underlying scheme ``X`` of type `AffineSchemeType` and local ring ``O_{(X,x)}`` of type `RingType` over some base ring ``k`` of type `BaseRingType`.
"""
@attributes mutable struct HypersurfaceGerm{
                 BaseRingType<:Ring,
                 RingType<:Ring,
                 AffineSchemeType<:AffineScheme} <: AbsSpaceGerm{BaseRingType, RingType}
  X::AffineSchemeType
  f::RingElem

  function HypersurfaceGerm(X::GermAtClosedPoint,f::MPolyLocRingElem; check::Bool=true)
    base_ring(modulus(OO(X))) == parent(f) || error("baserings do not match")
    @check begin
      (ideal(parent(f),[f]) == modulus(OO(X))) || error("given f does not define given X")
    end
    return new{typeof(base_ring(X)), typeof(OO(X)), typeof(X)}(X,f)
  end

## the following is currently unused....
## as no backend for groebner computations is currently available in this case
  function HypersurfaceGerm(X::GermAtGeometricPoint, f::MPolyLocRingElem; check::Bool=true)
    base_ring(modulus(OO(X))) == parent(f) || error("baserings do not match")
    @check begin
      (ideal(parent(f),[f]) == modulus(OO(X))) || error("given f does not define given X")
    end
    return new{typeof(base_ring(X)), typeof(OO(X)), typeof(X)}(X,f)
  end
end

@doc raw"""
    CompleteIntersectionGerm{BaseRingType, RingType, AffineSchemeType}

A complete intersection germ ``(X,O_{(X,x)}``, i.e. a ringed space with underlying scheme ``X`` of type AffineSchemeType and local ring ``O_{(X,x)}`` of type `RingType` over some base ring ``k`` of type `BaseRingType`.
"""
@attributes mutable struct CompleteIntersectionGerm{BaseRingType<:Ring, RingType<:Ring, AffineSchemeType<:AffineScheme} <: AbsSpaceGerm{BaseRingType, RingType}
  X::AffineSchemeType
  v::Vector{<:RingElem}

  function CompleteIntersectionGerm(X::GermAtClosedPoint, v::Vector{T}; check::Bool=true) where T<:MPolyLocRingElem
    R = base_ring(modulus(OO(X)))
    all(x->parent(x) == R, v) || error("base_rings do not coincide")
    @check begin
      length(v) == krull_dim(R) - dim(X) || error("not a complete intersection")
      modulus(OO(X)) == ideal(R,v) || error("given tuple does not generate modulus")
    end
    return new{typeof(base_ring(X)), typeof(OO(X)), typeof(X)}(X,v)
  end

## the following one is currently unused...
## as no backend for groebner computations is currently available in this case
  function CompleteIntersectionGerm(X::GermAtGeometricPoint, v::Vector{MPolyLocRingElem}; check::Bool=true)
    R = base_ring(OO(X))
    all(x->parent(x) == R, v) || error("base_rings do not coincide")
    @check begin
      length(v) == krull_dim(R) - dim(X) || error("not a complete intersection")
      modulus(OO(X)) == ideal(R,v) || error("given tuple does not generate modulus")
    end
    return new{typeof(base_ring(X)), typeof(OO(X)), typeof(X)}(X,v)
  end
end

##############################################################################
## Some more shorthand notation
##############################################################################
const AnySpaceGerm = Union{SpaceGerm, HypersurfaceGerm, CompleteIntersectionGerm}
const AnySpaceGermClosedPoint = Union{SpaceGerm{<:Ring,<:Ring,<:GermAtClosedPoint},
                                HypersurfaceGerm{<:Ring,<:Ring,<:GermAtClosedPoint},
                                CompleteIntersectionGerm{<:Ring,<:Ring,<:GermAtClosedPoint}}
const AnySpaceGermGeometricPoint = Union{SpaceGerm{<:Ring,<:Ring,<:GermAtGeometricPoint},
                                HypersurfaceGerm{<:Ring,<:Ring,<:GermAtGeometricPoint},
                                CompleteIntersectionGerm{<:Ring,<:Ring,<:GermAtGeometricPoint}}

##############################################################################
### Getter functions
##############################################################################

function underlying_scheme(X::AnySpaceGerm)
  return X.X
end

@doc raw"""
    representative(X::AnySpaceGerm)
    representative(X::SpaceGerm{<:Ring, <:MPolyLocRing})

Return a representative `Y` of a space germ `(X,p)` at a point `p`.

More precisely, let `(X,p)` be given by `Spec U^{-1}(R /I)`, where `R` is a polynomial
ring, `I` an ideal of it and `U` the complement of the maximal ideal corresponding
to `p. Then the representative `Y = Spec R/I` is returned.

# Examples
```jldoctest
julia> R, (x,y,z) = QQ[:x, :y, :z];

julia> I = ideal(R, [(x-1)*(x^2 - y^2 + z^2)]);

julia> X = spec(R, I);

julia> XS = SpaceGerm(X,[0,0,0])
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 3 variables x, y, z
        over rational field
      by ideal (x^3 - x^2 - x*y^2 + x*z^2 + y^2 - z^2)
    at complement of maximal ideal of point (0, 0, 0)

julia> representative(XS)
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x, y, z
      over rational field
    by ideal (x^3 - x^2 - x*y^2 + x*z^2 + y^2 - z^2)

julia> representative(XS) == X
true

julia> L, phi = localization(R,complement_of_point_ideal(R,[0,0,0]));

julia> IL = phi(I);

julia> Z = germ_at_point(quo(L,IL)[1])[1];

julia> representative(Z)
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x, y, z
      over rational field
    by ideal (x^3 - x^2 - x*y^2 + x*z^2 + y^2 - z^2)

```
"""
@attr AffineScheme function representative(X::AnySpaceGerm)
  return AffineScheme(underlying_quotient(OO(X)))
end

@attr AffineScheme function representative(X::SpaceGerm{<:Ring, <:MPolyLocRing})
  R = ambient_coordinate_ring(X)
  return AffineScheme(R)
end

@doc raw"""
    point(X::AnySpaceGermClosedPoint)
    point(X::AnySpaceGermGeometricPoint)

Return the point `p` of a germ `(X,p)`, where p is specified
- as its point_coordinates in the first case
- as the respective prime ideal of `p` in the second case

# Examples
```jldoctest
julia> R, (x,y,z) = QQ[:x, :y, :z];

julia> I = ideal(R, [(x-1)*(x^2 - y^2 + z^2)]);

julia> X = spec(R, I);

julia> XS = SpaceGerm(X,[0,0,0]);

julia> point(XS)
3-element Vector{QQFieldElem}:
 0
 0
 0

```
"""
function point(X::AnySpaceGermClosedPoint)
  return point_coordinates(inverted_set(OO(X)))
end

## TODO: move to higher level, i not already available
function point(X::AnySpaceGermGeometricPoint)
  return prime_ideal(inverted_set(OO(X)))
end

coordinate_ring(X::AbsSpaceGerm) = OO(X)

function defining_ideal(X::AnySpaceGerm)
    return modulus(OO(X))
end

function defining_ideal(X::SpaceGerm{<:Ring,<:MPolyLocRing})
    return ideal(OO(X),[])
end

@doc raw"""
    ambient_germ(X::AbsSpaceGerm)

Return the ambient germ of a given germ `(X,p)`.

More precisely, let `(X,p)` be given by `Spec U^{-1}(R /I)`, where `R` is a polynomial
ring, `I` an ideal of it and `U` the complement of the maximal ideal corresponding
to `p. Then the ambient germ `Spec U^{-1}R` is returned.

# Examples
```jldoctest
julia> R, (x,y,z) = QQ[:x, :y, :z];

julia> I = ideal(R, [(x-1)*(x^2 - y^2 + z^2)]);

julia> X = spec(R, I);

julia> XS = SpaceGerm(X,[0,0,0]);

julia> ambient_germ(XS)
Spectrum
  of localization
    of multivariate polynomial ring in 3 variables x, y, z
      over rational field
    at complement of maximal ideal of point (0, 0, 0)

```
"""
@attr SpaceGerm function ambient_germ(X::AnySpaceGerm)
    Y,_ = germ_at_point(localized_ring(OO(X)))
    return Y
end

@attr SpaceGerm function ambient_germ(X::SpaceGerm{<:Ring,<:MPolyLocRing})
    return X
end

@doc raw"""
    defining_ring_element(X::HypersurfaceGerm)
    defining_ring_elements(X::CompleteIntersectionGerm)

Return the (fixed) defining element(s) of the ideal of `X` in the ring of the ambient germ of `X`. Note that the return value is not an element of a polynomial ring, but of a localization of a polynomial ring the complement of a maximal ideal. (Hence each such element has a numerator and a denominator.)

Caution: This command is not exported and is only provided for convenience in programming.
# Examples:
```jldoctest
julia> R, (x,y,z) = QQ[:x, :y, :z];

julia> I = ideal(R, [(x-1)*(x^2 - y^2 + z^2)]);

julia> X = spec(R, I);

julia> XS = HypersurfaceGerm(X,[0,0,0])
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 3 variables x, y, z
        over rational field
      by ideal (x^3 - x^2 - x*y^2 + x*z^2 + y^2 - z^2)
    at complement of maximal ideal of point (0, 0, 0)

julia> defining_ring_element(XS)
-x^3 + x^2 + x*y^2 - x*z^2 - y^2 + z^2

```
"""
defining_ring_element(X::HypersurfaceGerm) = X.f::elem_type(localized_ring_type(ring_type(X)))
defining_ring_elements(X::CompleteIntersectionGerm) = X.v::Vector{elem_type(localized_ring_type(ring_type(X)))}

#######################################################################################
### constructors
#######################################################################################
@doc raw"""
    SpaceGerm(X::AbsAffineScheme, a::Vector{T}) where T<:Union{Integer, FieldElem}
    SpaceGerm(X::AbsAffineScheme, I:Ideal)
    SpaceGerm(X::AffineScheme(LocalRing))
    SpaceGerm(A::LocalRing)

Return the space germ `(X,p)` arising from the given representative `X` or the given
`X = Spec(A)` for a local ring `A`, where the point `p` may be specified in several
equivalent ways:
- by its coordinates `a` in the ambient_space of `X` or
- by a maximal ideal `I`in the coordinate ring of `X` or
- by a maximal ideal `I` in the ambient_coordinate_ring of `X`
- by the maximal ideal of the local ring `A`

!!!note
    Only `LocalRing`s localized at rational points over the coefficient field are currently fully supported.

# Examples
```jldoctest
julia> X = affine_space(QQ,3);

julia> R = coordinate_ring(X);

julia> (x,y,z) = gens(R);

julia> XL = SpaceGerm(X,ideal(R,[x,y,z]))
Spectrum
  of localization
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    at complement of maximal ideal of point (0, 0, 0)

julia> RL = coordinate_ring(XL);

julia> J = ideal(RL,[x^2-y^2+z^2])
Ideal generated by
  x1^2 - x2^2 + x3^2

julia> Q,_ = quo(RL,J);

julia> SpaceGerm(Q)
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 3 variables x1, x2, x3
        over rational field
      by ideal (x1^2 - x2^2 + x3^2)
    at complement of maximal ideal of point (0, 0, 0)

```
"""
function SpaceGerm(X::AbsAffineScheme, a::Vector{T}) where T<:Union{Integer, FieldElem}
  R = ambient_coordinate_ring(X)
  kk = coefficient_ring(R)
  b = [kk.(v) for v in a]  ## throws an error, if vector entries are not compatible
  U = MPolyComplementOfKPointIdeal(R,b)
  Y = spec(localization(OO(X), U)[1])
  Z = SpaceGerm(Y)
  set_attribute!(Z,:representative,X)
  return SpaceGerm(Y)
end

function SpaceGerm(X::AbsAffineScheme, I::MPolyIdeal)
  R = base_ring(I)
  R === ambient_coordinate_ring(X) || error("rings are not compatible")
  a = rational_point_coordinates(I)
  Y = SpaceGerm(X,a)
  return Y
end

function SpaceGerm(X::AbsAffineScheme, I::Ideal)
  A = base_ring(I)
  A === OO(X) || error("rings are incompatible")
  J = saturated_ideal(I)
  a = rational_point_coordinates(J)
  return SpaceGerm(X,a)
end

@doc raw"""
    SpaceGerm(X::AbsAffineScheme, p::AbsAffineRationalPoint)
    SpaceGerm(p::AbsAffineRationalPoint)

Return a space germ `(X,p)` for a given `X` and a rational point `p` on some affine scheme `Y`. If no `X` is specified, `Y` is used in the place of `Y`.
"""
SpaceGerm(p::AbsAffineRationalPoint) = SpaceGerm(codomain(p), coordinates(p))

function SpaceGerm(X::AbsAffineScheme, p::AbsAffineRationalPoint)
  ambient_space(X) == ambient_space(codomain(p)) || error("ambient spaces do not match")
  return SpaceGerm(X,coordinates(p))
end

@doc raw"""
    HypersurfaceGerm(X::AbsAffineScheme, a::Vector{T}) where T<:Union{Integer, FieldElem}
    HypersurfaceGerm(X::AbsAffineScheme, I:Ideal)
    HypersurfaceGerm(A::LocalRing)

Check that `X` (or `Spec(A)` respectively) represents a hypersurface germ at the given
point `p` and returns the hypersurface germ `(X,p)` from `X` in the affirmative case, where `p`
may be specified in several equivalent ways:
- by its coordinates `a` in the ambient_space of `X` or
- by a maximal ideal `I` in the coordinate ring of `X`  or
- by a maximal ideal `I` in the ambient_coordinate_ring of `X`
- by the maximal ideal of the local ring `A`

    HypersurfaceGerm(X::AffineScheme(LocalRing),f::MPolyLocRingElem)

This variant allows explicit specification of the generator for the hypersurface. The given `f` is checked to generate to modulus of OO(X) or A respectively. In the affirmative case, the given generator will subsequently be used by all methods explicitly accessing a generator.

!!!note
    Only `LocalRing`s localized at rational points over the coefficient field are currently fully supported.

# Examples
```jldoctest
julia> X = affine_space(QQ,3);

julia> R = coordinate_ring(X);

julia> (x,y,z) = gens(R);

julia> Q,_ = quo(R,ideal(R,[x^2+y^2+z^2]));

julia> HypersurfaceGerm(spec(Q),[0,0,0])
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 3 variables x1, x2, x3
        over rational field
      by ideal (x1^2 + x2^2 + x3^2)
    at complement of maximal ideal of point (0, 0, 0)

```
"""
function HypersurfaceGerm(X::AbsAffineScheme, a::Vector{T}) where T<:Union{Integer, FieldElem}
  R = ambient_coordinate_ring(X)
  kk = coefficient_ring(R)
  b = [kk.(v) for v in a]  ## throws an error, if vector entries are not compatible
  U = MPolyComplementOfKPointIdeal(R,b)
  LX,_ = localization(OO(X), U)
  mingens = minimal_generating_set(modulus(LX))
  length(mingens) == 1 || error("not a hypersurface")
  f = mingens[1]
  Y = HypersurfaceGerm(spec(LX),f)
  set_attribute!(Y,:representative,X)
  return Y
end

function HypersurfaceGerm(X::AbsAffineScheme, I::MPolyIdeal)
  R = base_ring(I)
  R === ambient_coordinate_ring(X) || error("rings are not compatible")
  a = rational_point_coordinates(I)
  Y = HypersurfaceGerm(X,a)
  return Y
end

function HypersurfaceGerm(X::AbsAffineScheme, I::Ideal)
  A = base_ring(I)
  A === OO(X) || error("rings are incompatible")
  J = saturated_ideal(I)
  a = rational_point_coordinates(J)
  return HypersurfaceGerm(X,a)
end

@doc raw"""
    HypersurfaceGerm(X::AbsAffineScheme, p::AbsAffineRationalPoint)
    HypersurfaceGerm(p::AbsAffineRationalPoint)

Return a hypersurface germ `(X,p)` for a given `X` and a rational point `p` on some scheme `Y`. If no `X` is specified, `Y` is used in its place.
"""
HypersurfaceGerm(p::AbsAffineRationalPoint) = HypersurfaceGerm(codomain(p), coordinates(p))

function HypersurfaceGerm(X::AbsAffineScheme, p::AbsAffineRationalPoint)
  ambient_space(X) == ambient_space(codomain(p)) || error("ambient spaces do not match")
  return HypersurfaceGerm(X,coordinates(p))
end

@doc raw"""
    CompleteIntersectionGerm(X::AbsAffineScheme, a::Vector{T}) where T<:Union{Integer, FieldElem}
    CompleteIntersectionGerm(X::AbsAffineScheme, I:Ideal)

Check that `X` represents a complete intersection germ at the given point `p` and returns
the complete intersection germ `(X,p)` from `X` in the affirmative case, where `p` may
be specified in several equivalent ways:
- by its coordinates `a` in the ambient_space of `X` or
- by a maximal ideal in the coordinate ring of `X` or
- by a maximal ideal in the ambient_coordinate_ring of `X`

# Examples
```jldoctest
julia> X = affine_space(QQ,3);

julia> R = coordinate_ring(X);

julia> (x,y,z) = gens(R);

julia> Q,_ = quo(R,ideal(R,[x^2+y^2+z^2,x*y]));

julia> CompleteIntersectionGerm(spec(Q),[0,0,0])
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 3 variables x1, x2, x3
        over rational field
      by ideal (x1^2 + x2^2 + x3^2, x1*x2)
    at complement of maximal ideal of point (0, 0, 0)

```
"""
function CompleteIntersectionGerm(X::AbsAffineScheme, a::Vector{T}) where T<:Union{Integer, FieldElem}
  R = ambient_coordinate_ring(X)
  kk = coefficient_ring(R)
  b = [kk.(v) for v in a]  ## throws an error, if vector entries are not compatible
  U = MPolyComplementOfKPointIdeal(R,b)
  L,_ = localization(OO(X), U)
  mingens = minimal_generating_set(modulus(L))
  length(mingens) == krull_dim(R) - dim(L) || error("not a complete intersection")
  w = mingens
  Y = CompleteIntersectionGerm(spec(L),w)
  set_attribute!(Y,:representative,X)
  return Y
end

function CompleteIntersectionGerm(X::AbsAffineScheme, I::MPolyIdeal)
  R = base_ring(I)
  R === ambient_coordinate_ring(X) || error("rings are not compatible")
  a = rational_point_coordinates(I)
  Y = CompleteIntersectionGerm(X,a)
  return Y
end

function ComleteIntersectionGerm(X::AbsAffineScheme, I::Ideal)
  A = base_ring(I)
  A === OO(X) || error("rings are incompatible")
  J = saturated_ideal(I)
  a = rational_point_coordinates(J)
  return CompleteIntersectionGerm(X,a)
end

@doc raw"""
    CompleteIntersectionGerm(X::AbsAffineScheme, p::AbsAffineRationalPoint)
    CompleteIntersectionGerm(p::AbsAffineRationalPoint)

Return a complete intersection germ `(X,p)` for a given `X`and a rational point `p` on some affine scheme `Y`, provided that $X$ is locally a complete intersection in some neighborhood of `p`. If no `X` is specified, `Y` is used in its place.
"""
CompleteIntersectionGerm(p::AbsAffineRationalPoint) = CompleteIntersectionGerm(codomain(p), coordinates(p))

function CompleteIntersectionGerm(X::AbsAffineScheme, p::AbsAffineRationalPoint)
  ambient_space(X) == ambient_space(codomain(p)) || error("ambient spaces do not match")
  return CompleteIntersectionGerm(X,coordinates(p))
end

@doc raw"""
    germ_at_point(X::AbsAffineScheme, I::Union{Ideal,Vector})
    germ_at_point(X::AbsAffineScheme, I::Union{Ideal,Vector})
    germ_at_point(A::Union{MPolyRing,MPolyQuoRing},
              I::Union{Ideal,Vector})
    germ_at_point(A::LocalRing, I::Union{Ideal,Vector})

Return a SpaceGerm `(X,p)` and the corresponding inclusion morphism of spectra
arising from the given representative `X` or the given
`X = Spec(A)` for a local ring `A`, where the point `p` may be specified in several
equivalent ways:
- by its coordinates `a` in the ambient_space of `X` or
- by a maximal ideal `I`in the coordinate ring of `X` or
- by a maximal ideal `I` in the ambient_coordinate_ring of `X`
- by the maximal ideal of the local ring `A`

!!!note
    Only `LocalRing`s localized at rational points over the coefficient field are currently fully supported.

# Examples

```jldoctest
julia> X = affine_space(QQ,3);

julia> R = coordinate_ring(X);

julia> (x,y,z) = gens(R);

julia> Q,_ = quo(R,ideal(R,[x^2+y^2+z^2]));

julia> Y,phi = germ_at_point(Q,[0,0,0]);

julia> Y
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 3 variables x1, x2, x3
        over rational field
      by ideal (x1^2 + x2^2 + x3^2)
    at complement of maximal ideal of point (0, 0, 0)

julia> phi
Affine scheme morphism
  from [x1, x2, x3]  Spec of localization of Q at complement of maximal ideal
  to   [x1, x2, x3]  scheme(x1^2 + x2^2 + x3^2)
given by the pullback function
  x1 -> x1
  x2 -> x2
  x3 -> x3

```
"""
function germ_at_point(X::AbsAffineScheme, I::Union{Ideal,Vector})
  Y = SpaceGerm(X, I)
  restr_map = morphism(Y, X, hom(OO(X), OO(Y), gens(OO(Y)), check=false), check=false)
  return Y, restr_map
end

germ_at_point(A::Union{MPolyRing,MPolyQuoRing},
              I::Union{Ideal,Vector}) = germ_at_point(spec(A),I)

@doc raw"""
    germ_at_point(X::AbsAffineScheme, p::AbsAffineRationalPoint)
    germ_at_point(p::AbsAffineRationalPoint)

Return a space germ `(X,p)` and the corresponding inclusion morphism of spectra arising
from the representative `X` for a given `X` and a rational point `p` on some affine scheme `Y`. If no `X` is specified, `Y` is used in its place..

"""
germ_at_point(p::AbsAffineRationalPoint) = germ_at_point(codomain(p), coordinates(p))

function germ_at_point(X::AbsAffineScheme, p::AbsAffineRationalPoint)
  ambient_space(X) == ambient_space(codomain(p)) || error("ambient spaces do not match")
  return germ_at_point(X,coordinates(p))
end

@doc raw"""
    hypersurface_germ(X::AbsAffineScheme, I::Union{Ideal,Vector})
    hypersurface_germ(X::AbsAffineScheme, I::Union{Ideal,Vector})
    hypersurface_germ(A::Union{MPolyRing,MPolyQuoRing},
              I::Union{Ideal,Vector})
    hypersurface_germ(A::LocalRing, I::Union{Ideal,Vector})

Return a HypersurfaceGerm `(X,p)` and the corresponding inclusion morphism of spectra
arising from the given representative `X` or the given
`X = Spec(A)` for a local ring `A`, where the point `p` may be specified in several
equivalent ways:
- by its coordinates `a` in the ambient_space of `X` or
- by a maximal ideal `I`in the coordinate ring of `X` or
- by a maximal ideal `I` in the ambient_coordinate_ring of `X`
- by the maximal ideal of the local ring `A`

!!!note
    Only `LocalRing`s localized at rational points over the coefficient field are currently fully supported.

!!!note
    If the defining ideal of `(X,p)` is not principal, an error exception occurs.

# Examples
```jldoctest
julia> X = affine_space(QQ,3);

julia> R = coordinate_ring(X);

julia> (x,y,z) = gens(R);

julia> Q,_ = quo(R,ideal(R,[x^2+y^2+z^2]));

julia> Y,phi = hypersurface_germ(Q,[0,0,0]);

julia> Y
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 3 variables x1, x2, x3
        over rational field
      by ideal (x1^2 + x2^2 + x3^2)
    at complement of maximal ideal of point (0, 0, 0)

julia> phi
Affine scheme morphism
  from [x1, x2, x3]  Spec of localization of Q at complement of maximal ideal
  to   [x1, x2, x3]  scheme(x1^2 + x2^2 + x3^2)
given by the pullback function
  x1 -> x1
  x2 -> x2
  x3 -> x3

```
"""
function hypersurface_germ(X::AbsAffineScheme, I::Union{Ideal,Vector})
  Y = HypersurfaceGerm(X,I)
  restr_map = morphism(Y, X, hom(OO(X), OO(Y), gens(OO(Y)), check=false), check=false)
  return Y, restr_map
end

hypersurface_germ(A::Union{MPolyRing,MPolyQuoRing},
                  I::Union{Ideal,Vector}) = hypersurface_germ(spec(A),I)

@doc raw"""
    hypersurface_germ(X::AbsAffineScheme, p::AbsAffineRationalPoint)
    hypersurface_germ(p::AbsAffineRationalPoint)

Return a hypersurface germ `(X,p)` and the corresponding inclusion morphism of spectra for a given `X` and a rational point `p` on some affine scheme `Y`. If no `X` is specified, `Y` is used in its place.

!!!note
    If the defining ideal of `(X,p)` is not principal. an error exception occurs.
"""
hypersurface_germ(p::AbsAffineRationalPoint) = hypersurface_germ(codomain(p), coordinates(p))

function hypersurface_germ(X::AbsAffineScheme, p::AbsAffineRationalPoint)
  ambient_space(X) == ambient_space(codomain(p)) || error("ambient spaces do not match")
  return hypersurface_germ(X,coordinates(p))
end

@doc raw"""
    complete_intersection_germ(X::AbsAffineScheme, I::Union{Ideal,Vector})
    complete_intersection_germ(X::AbsAffineScheme, I::Union{Ideal,Vector})
    complete_intersection_germ(A::Union{MPolyRing,MPolyQuoRing},
              I::Union{Ideal,Vector})
    complete_intersection_germ(A::LocalRing, I::Union{Ideal,Vector})

Return a CompleteIntersectionGerm `(X,p)` and the corresponding inclusion morphism of spectra arising from the given representative `X` or the given `X = Spec(A)` for a local ring `A`, where the point `p` may be specified in several
equivalent ways:
- by its coordinates `a` in the ambient_space of `X` or
- by a maximal ideal `I`in the coordinate ring of `X` or
- by a maximal ideal `I` in the ambient_coordinate_ring of `X`
- by the maximal ideal of the local ring `A`

!!!note: Only `LocalRing`s localized at rational points over the coefficient field are currently fully supported.

!!!note: If the defining ideal of `(X,p)` is not principal, an error exception occurs.

# Examples
```jldoctest
julia> X = affine_space(QQ,3);

julia> R = coordinate_ring(X);

julia> (x,y,z) = gens(R);

julia> Q,_ = quo(R,ideal(R,[x^2+y^2+z^2,x*y]));

julia> Y,phi = complete_intersection_germ(Q,[0,0,0]);

julia> Y
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 3 variables x1, x2, x3
        over rational field
      by ideal (x1^2 + x2^2 + x3^2, x1*x2)
    at complement of maximal ideal of point (0, 0, 0)

julia> phi
Affine scheme morphism
  from [x1, x2, x3]  Spec of localization of Q at complement of maximal ideal
  to   [x1, x2, x3]  scheme(x1^2 + x2^2 + x3^2, x1*x2)
given by the pullback function
  x1 -> x1
  x2 -> x2
  x3 -> x3

```
"""
function complete_intersection_germ(X::AbsAffineScheme, I::Union{Ideal,Vector})
  Y = CompleteIntersectionGerm(X,I)
  restr_map = morphism(Y, X, hom(OO(X), OO(Y), gens(OO(Y)), check=false), check=false)
  return Y, restr_map
end

complete_intersection_germ(A::Union{MPolyRing,MPolyQuoRing},
                  I::Union{Ideal,Vector}) = complete_intersection_germ(spec(A),I)

@doc raw"""
    complete_intersection_germ(X::AbsAffineScheme, p::AbsAffineRationalPoint)
    complete_intersection_germ(p::AbsAffineRationalPoint)

Return a complete intersection germ `(X,p)` and the corresponding inclusion morphism of spectra for a given `X` and a rational point `p` on some affine scheme `Y`. If no `X` is specified, `Y` is used in its place.

!!!note
    If the defining ideal of `(X,p)` does not describe a complete intersection. an error exception occurs.
"""
complete_intersection_germ(p::AbsAffineRationalPoint) = hypersurface_germ(codomain(p), coordinates(p))

function complete_intersection_germ(X::AbsAffineScheme, p::AbsAffineRationalPoint)
  ambient_space(X) == ambient_space(codomain(p)) || error("ambient spaces do not match")
  return complete_intersection_germ(X,coordinates(p))
end

#########################################################################################
## for convenience of users thinking in terms of local rings
#########################################################################################

const LocalRing = Union{
                  MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any,
                                        <:MPolyComplementOfKPointIdeal},
                  MPolyLocRing{<:Any, <:Any, <:Any, <:Any,
                                     <:MPolyComplementOfKPointIdeal},
                  MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any,
                                        <:MPolyComplementOfPrimeIdeal},
                  MPolyLocRing{<:Any, <:Any, <:Any, <:Any,
                                     <:MPolyComplementOfPrimeIdeal}
                 }


function SpaceGerm(A::LocalRing)
  return SpaceGerm(spec(A))
end

function HypersurfaceGerm(A::LocalRing)
  I = modulus(A)
  v  = minimal_generating_set(I)
  length(v) == 1 || error("not a hypersurface germ")
  return HypersurfaceGerm(spec(A),v[1])
end

function CompleteIntersectionGerm(A::LocalRing)
  I = modulus(A)
  !iszero(I) || error("zero ideal not allowed for complete intersection germ")
  v = minimal_generating_set(I)
  length(v) == krull_dim(base_ring(I)) - krull_dim(A) || error("not a complete intersection germ")
  return CompleteIntersectionGerm(spec(A),v)
end

## and with identity map to keep usage consistent
function germ_at_point(A::LocalRing)
  X = SpaceGerm(A)
  restr_map = morphism(X, X, hom(OO(X), OO(X), gens(OO(X)), check=false), check=false)
  return X, restr_map
end

function hypersurface_germ(A::LocalRing)
  X = HypersurfaceGerm(A)
  restr_map = morphism(X, X, hom(OO(X), OO(X), gens(OO(X)), check=false), check=false)
  return X, restr_map
end

function complete_intersection_germ(A::LocalRing)
  X = CompleteIntersectionGerm(A)
  restr_map = morphism(X, X, hom(OO(X), OO(X), gens(OO(X)), check=false), check=false)
  return X, restr_map
end

### basic functionality for space germs

##############################################################################
# note: ==, issubset are basically inherited from AffineScheme
#       intersect uses explicit fallback to AffineScheme and adjusted return types
#       union uses explicit fallback and adjusted return types
##############################################################################
#@doc raw"""
#    is_subset(X::AbsSpaceGerm,Y::AbsSpaceGerm)
#
#Return whether `X` is a subset of `Y` as space germs
#"""
function is_subset(X::AbsSpaceGerm{<:Any, <:MPolyQuoLocRing}, Y::AbsSpaceGerm{<:Any, <:MPolyQuoLocRing})
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  point(X) == point(Y) || return false
  IY=ideal(localized_ring(OO(X)), gens(modulus(underlying_quotient(OO(Y)))))
  return issubset(IY,modulus(OO(X)))
end

function is_subset(X::AbsSpaceGerm{<:Any, <:MPolyLocRing}, Y::AbsSpaceGerm{<:Any, <:MPolyQuoLocRing})
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  point(X) == point(Y) || return false
  return iszero(modulus(OO(Y)))
end

function is_subset(X::AbsSpaceGerm{<:Any, <:MPolyLocRing}, Y::AbsSpaceGerm{<:Any, <:MPolyLocRing})
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  return point(X) == point(Y)
end

function is_subset(X::AbsSpaceGerm{<:Any, <:MPolyQuoLocRing}, Y::AbsSpaceGerm{<:Any, <:MPolyLocRing})
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  return point(X) == point(Y)
end

## note: intersection of hypersurfaces is hardly ever a hypersurface
##       intersection of complete intersections need not be complete intersection
##       hence return type always SpaceGerm
@doc raw"""
    intersect(X::AbsSpaceGerm,Y::AbsSpaceGerm) --> SpaceGerm

Return the intersection of `X` and `Y`, provided they are germs at the same point.
"""
function Base.intersect(X::AbsSpaceGerm, Y::AbsSpaceGerm)
  point(X) == point(Y) || error("not the same point of the germ")
  Z = intersect(underlying_scheme(X),underlying_scheme(Y))
  return SpaceGerm(Z)
end

@doc raw"""
    union(X::AbsSpaceGerm,Y::AbsSpaceGerm) --> SpaceGerm
    union(X::HypersurfaceGerm, Y:: HypersurfaceGerm) --> HypersurfaceGerm

Return the union of `X`and `Y`. If `X`and `Y` happen to be HypersurfaceGerms, so is the result.

# Examples
```jldoctest
julia> X = affine_space(QQ,3);

julia> R = coordinate_ring(X);

julia> (x,y,z) = gens(R);

julia> Q1,_ = quo(R,ideal(R,[x*y]));

julia> Q2,_ = quo(R,ideal(R,[x*(x-y)]));

julia> X1 = HypersurfaceGerm(spec(Q1),[0,0,0]);

julia> X2 = HypersurfaceGerm(spec(Q2),[0,0,0]);

julia> union(X1,X2)
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 3 variables x1, x2, x3
        over rational field
      by ideal (x1^2*x2 - x1*x2^2)
    at complement of maximal ideal of point (0, 0, 0)

```
"""
function Base.union(X::AbsSpaceGerm, Y::AbsSpaceGerm)
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || error("not subgerms of a common space germ")
  point(X) == point(Y) || error("not the same point of the germ")
  # comparison of points implicitly also checks that localization was performed at points
  # otherwise 'point' is not implemented
  I = intersect(modulus(underlying_quotient(OO(X))),modulus(underlying_quotient(OO(Y))))
  Z,_ = germ_at_point(MPolyQuoLocRing(R, I ,inverted_set(OO(X))))
  return Z
end

## note: union of hypersurface germs is again a hypersurface germ
##       union of complete intersection germs need not even be equidimensional
function Base.union(X::HypersurfaceGerm, Y::HypersurfaceGerm)
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || error("not subgerms of a common space germ")
  point(X) == point(Y) || error("not the same point of the germ")
  # comparison of points implicitly also checks that localization was performed at points
  # otherwise 'point' is not implemented
  f_new = numerator(defining_ring_element(X))*numerator(defining_ring_element(Y))
  f_new = radical(ideal(R,[f_new]))[1]
  Y = HypersurfaceGerm(spec(quo(R,ideal(R,f_new))[1]), point(X))
  Y.f = f_new
  return Y
end

##############################################################################
# note: singular_locus, is_smooth and is_regular are inherited from AffineScheme
##############################################################################

# We want the singular locus of an `AbsSpaceGerm` to be a `SpaceGerm` again and
# not a plain `AffineScheme`.
@doc raw"""
    singular_locus(X::AbsSpaceGerm)  --> SpaceGerm, ClosedEmbedding

Return the space germ (Y,p) for a given germ (X,p) and the closed embedding of (Y,p) into (X,p), where Y is the singular locus of X.

# Examples
```jldoctest
julia> X = affine_space(QQ,3);

julia> R = coordinate_ring(X);

julia> (x,y,z) = gens(R);

julia> Q,_ = quo(R,ideal(R,[x^2+y^2]));

julia> Y = SpaceGerm(spec(Q),[0,0,0]);

julia> XSL, incSL = singular_locus(Y);

julia> XSL
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 3 variables x1, x2, x3
        over rational field
      by ideal (x1^2 + x2^2, x2, x1)
    at complement of maximal ideal of point (0, 0, 0)

julia> incSL
Affine scheme morphism
  from [x1, x2, x3]  Spec of localization of quotient of multivariate polynomial ring at complement of maximal ideal
  to   [x1, x2, x3]  Spec of localization of Q at complement of maximal ideal
given by the pullback function
  x1 -> x1
  x2 -> x2
  x3 -> x3

```
"""
function singular_locus(X::AbsSpaceGerm)
  S, inc = singular_locus(underlying_scheme(X))
  Sgerm = SpaceGerm(S)
  return Sgerm, ClosedEmbedding(morphism(Sgerm, X, pullback(inc), check=false), image_ideal(inc), check=false)
end

## note: subgerms of hypersurface and complete intersection germs are simply space germs
@doc raw"""
    subgerm(X::AbsSpaceGerm, I::Ideal) --> SpaceGerm

Return the space germ (Y,p) of (X,p) defined by the ideal I in the local ring of X at p

!!! note: (Y,p) is of type SpaceGerm, even if (X,p) is a HypersurfaceGerm or a CompleteIntersectionGerm.
"""
function subgerm(X::AbsSpaceGerm, I::Ideal)
  base_ring(I) === OO(X) || error("ideal does not belong to the correct ring")
  Y = subscheme(underlying_scheme(X), I)
  return SpaceGerm(Y)
end

subscheme(X::AbsSpaceGerm, I::Ideal) = subgerm(X::AbsSpaceGerm, I::Ideal)

@doc raw"""
    is_isolated_singularity(X::AbsSpaceGerm)

Return whether `(X,p)` has at most an isolated singularity.
```jldoctest
julia> X = affine_space(QQ,3);

julia> R = coordinate_ring(X);

julia> (x,y,z) = gens(R);

julia> Q,_ = quo(R,ideal(R,[x^2+y^2]));

julia> Y = SpaceGerm(spec(Q),[0,0,0]);

julia> is_isolated_singularity(Y)
false

```
"""
@attr Bool function is_isolated_singularity(X::AbsSpaceGerm)
  return dim(singular_locus(X)[1]) < 1
end

##############################################################################
# milnor_number, milnor_algebra for IHS
# milnor_number for ICIS
# -- beyond this we do no longer have a bouquet of spheres of same dimension
##############################################################################
@doc raw"""
    milnor_algebra(X::HypersurfaceGerm)

Return the local Milnor algebra of `(X,p)` at p
"""
function milnor_algebra(X::HypersurfaceGerm)
  R = localized_ring(OO(X))
  ## milnor number independent of choice of representative
  ## hence choose a polynomial representative for easier computation
  f_poly = numerator(defining_ring_element(X))
  I = ideal(R, R.([derivative(f_poly, i) for i=1:nvars(base_ring(R))]))
  return quo(R,I)[1]
end

@doc raw"""
    milnor_number(X::HypersurfaceGerm)
    milnor_number(X::CompleteIntersectionGerm)

Return the local Milnor number of `(X,p)` at p
"""
function milnor_number(X::HypersurfaceGerm)
  return vector_space_dim(milnor_algebra(X))
end

function milnor_number(X::CompleteIntersectionGerm)
  R = localized_ring(OO(X))
  ## milnor number independent of choice of representative
  ## hence choose polynomial representatives for easier computation
  v = [numerator(a) for a in defining_ring_elements(X)]
  w = typeof(v[1])[]   ## already used entries of v
  dims = 0             ## for building up the alternating sum
  sign = 1             ##     and the sign

  ## alternating sum of the Le-Greuel formula, one summand per while-loop pass
  while !is_empty(v)
    found = false      ## keep track of 'colength condition satisfied?'

    ## run through the potential choices, until colength condition satisfied
    for f in v
      ## note: although we have already moved to polynomial data, we need to
      ##       hand localized ring to helper to ensure shift to origin
      ##       before computations
      dtemp = _icis_milnor_helper(R,w,f)   ## colength; -1 if condition violated
      if dtemp > -1
        dims = dims + sign * dtemp         ## contribute to alternating sum
        sign = -sign
        push!(w,f)                          ## put f in 'used' list
        deleteat!(v, findfirst(==(f), v)) ## remove f from 'unused' list
        found = true
        break
      end
    end
    found == true || error("retry/general linear combinations not implemented yet")
  end
  if dims > 0
    return dims
  end
  return -dims
end

## internal helper for the summands in the Le-Greuel formula -- local case
function _icis_milnor_helper(L::MPolyLocRing, v::Vector,f::RingElem)
  R = parent(f)
  R == base_ring(L) || error("base_rings do not match")
  all(a -> parent(a)==R,v) || error("base_rings do not match")

  ## establish the shift to origin to allow computations w.r.t. local ordering
  shift,back_shift = base_ring_shifts(L)
  w = [shift(a) for a in v]
  g = shift(f)

  ## compute the appropriate summand in the Le-Greuel formula
  ## for the (length(w))-th contribution with specified choice
  I = ideal(R,w)
  push!(w,g)
  n = nvars(R)
  JM = matrix(R,n, length(w),
              [derivative(h,i) for i=1:n for h in w])
  mo = minors(JM,length(w))
  J = I + ideal(R,mo)
  !isone(J) || return 0
  o = negdegrevlex(gens(R))
  LJ = leading_ideal(J;ordering=o)

  ## we might have a violated colength condition (i.e. dim(LJ)>0)
  krull_dim(quo(R,LJ)[1]) == 0 || return (-1)
  return vector_space_dim(LJ)
end

@doc raw"""
    milnor_algebra(X::AffineScheme{<:Field,<:MPolyQuoRing})

Return the global milnor algebra of the affine hypersurface `X`. If `X` is not a hypersurface, an error occurs.
"""
function milnor_algebra(X::AffineScheme{<:Field,<:MPolyQuoRing})
  R = base_ring(OO(X))
  ngens(modulus(OO(X))) == 1 || error("not a hypersurface (or unnecessary generators in specified generating set)")
  v = gen(modulus(OO(X)),1)
  I = ideal(R,R.([derivative(v,i) for i in 1:nvars(R)]))
  return quo(R,I)[1]
end

@doc raw"""
    milnor_number(X::AffineScheme{<:Field,<:MPolyQuoRing})

Return the global milnor number of the affine hypersurface or complete intersection `X`. If `X` is neither of the two, an error occurs.
"""
function milnor_number(X::AffineScheme{<:Field,<:MPolyQuoRing})
  R = base_ring(OO(X))
  v = gens(modulus(OO(X)))
  if length(v) == 1
    return vector_space_dim(milnor_algebra(X))
  end
  length(v) == krull_dim(R) - dim(X) || error("not a complete intersection (or unnecessary generators in specified generating set)")
  w = typeof(v[1])[]   ## already used entries of v
  dims = 0             ## for building up the alternating sum
  sign = 1
  while !is_empty(v)
    found = false      ## colength condition satisfied?
    for f in v
    ## note: in the global case, we do not need to shift
    ##       hence helper with completely different signature
      dtemp = _icis_milnor_helper(w,f)   ## colength; -1 if condition violated
      if dtemp > -1
        dims = dims + sign * dtemp         ## alternating sum
        sign = -sign
        push!(w,f)                          ## put f in 'used' list
        deleteat!(v, findfirst(==(f), v)) ## remove f from 'unused' list
        found = true
        break
      end
    end
    found == true || error("retry/general linear combinations not implemented yet")
  end
  if dims > 0
     return dims
  end
  return -dims
end

## internal helper for the summands in the Le-Greuel formula -- global case
function _icis_milnor_helper(v::Vector,f::MPolyRingElem)
  R = parent(f)
  all(a-> parent(a) == R,v) || error("base rings do not match")

  ## compute the appropriate step in the Le-Greuel formula
  ## for the (length(w))-th contribution
  I = ideal(R,v)
  bla = copy(v)
  push!(bla,f)
  n = nvars(R)
  JM = matrix(R,n, length(bla),
              [derivative(h,i) for i in 1:n for h in bla])
  mo = minors(JM,length(bla))
  J = I + ideal(R,mo)
  o = degrevlex(gens(R))
  LJ = leading_ideal(J;ordering=o)

  ## we might have a violated colength condition (i.e. dim(LJ)>0)
  dim(LJ) == 0 || return -1
  return vector_space_dim(quo(R,LJ)[1])
end
