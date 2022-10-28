export affine_space
export closure
export hypersurface_complement
export standard_spec
export subscheme


########################################################
# (1) Generic constructors
########################################################

Spec(R::MPolyRing, I::MPolyIdeal) = Spec(quo(R, I)[1])


Spec(R::MPolyRing, U::AbsMPolyMultSet) = Spec(Localization(R, U)[1])


Spec(R::MPolyRing, I::MPolyIdeal, U::AbsMPolyMultSet) = Spec(MPolyQuoLocalizedRing(R, I, U))



########################################################
# (2) Copy constructors
########################################################

Spec(X::Spec) = Spec(OO(X))

### internal copy routines
Base.deepcopy_internal(X::Spec, dict::IdDict) = Spec(deepcopy_internal(OO(X), dict))



########################################################
# (3) Affine n-dimensional space
########################################################

function affine_space(kk::BRT, n::Int; variable_name="x") where {BRT<:Ring}
  R, _ = PolynomialRing(kk, [ variable_name * "$i" for i in 1:n])
  return Spec(R)
end

function affine_space(kk::BRT, var_symbols::Vector{Symbol}) where {BRT<:Ring}
  R, _ = PolynomialRing(kk, var_symbols)
  return Spec(R)
end



########################################################
# (4) StdSpec (needed?)
########################################################

standard_spec(X::AbsSpec{<:Any, <:MPolyRing}) = Spec(MPolyQuoLocalizedRing(OO(X), ideal(OO(X), [zero(OO(X))]), units_of(OO(X))))

function standard_spec(X::AbsSpec{<:Any, <:MPolyQuo})
  A = OO(X)
  R = base_ring(A)
  return Spec(MPolyQuoLocalizedRing(R, modulus(A), units_of(R)))
end

standard_spec(X::AbsSpec{<:Any, <:MPolyLocalizedRing}) = Spec(MPolyQuoLocalizedRing(ambient_ring(X), ideal(ambient_ring(X), [zero(ambient_ring(X))]), inverted_set(OO(X))))
standard_spec(X::AbsSpec{<:Any, <:MPolyQuoLocalizedRing}) = Spec(OO(X))



########################################################
# (5) Closed subschemes
########################################################

## Closed subscheme defined by elements
function subscheme(X::Spec, f::Vector{<:RingElem})
  all(x->(parent(x) == OO(X)), f) || return subscheme(X, OO(X).(f))
  return subscheme(X, ideal(OO(X), f))
end


## Closed subscheme defined by ideal
#
# Note that these default to the plain type Spec and must be overwritten 
# if something more sophisticated should be returned!
@Markdown.doc """
    subscheme(X::AbsSpec, I::Ideal)

For a scheme ``X = Spec(R)`` and an ideal ``I ⊂ 𝒪(X)`` 
this returns the closed subscheme defined by ``I``.
"""
function subscheme(X::AbsSpec, I::Ideal)
  base_ring(I) == OO(X) || return subscheme(X, ideal(OO(X), OO(X).(gens(I)))) # this will throw if coercion is not possible
  return Spec(quo(OO(X), I)[1])
end
subscheme(X::AbsSpec, f::RingElem) = subscheme(X, ideal(OO(X), [f]))
subscheme(X::AbsSpec, f::Vector{<:RingElem}) = subscheme(X, ideal(OO(X), f))



########################################################
# (6) Open subschemes
########################################################

### open subschemes defined by complements of hypersurfaces
@Markdown.doc """
    hypersurface_complement(X::AbsSpec, f::RingElem)

For a scheme ``X = Spec(R)`` and an element ``f ∈ R`` 
this returns the open subscheme ``U = Spec(R[f⁻¹]) = X ∖ V(f)`` 
defined by the complement of the vanishing 
locus of ``f``.
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

Base.intersect(E::EmptyScheme{BRT}, X::Scheme{BRT}) where {BRT} = E


Base.intersect(X::Scheme{BRT}, E::EmptyScheme{BRT}) where {BRT} = E


Base.intersect(X::EmptyScheme{BRT}, E::EmptyScheme{BRT}) where {BRT} = E


### For Specs of MPolyRings
# TODO  intersect X,Y for X<Y should return a copy of X with === ambient_rings
# Spec(X) does not apply for instance to principal open subsets hence a change
# is necessary
function Base.intersect(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyRing}
  ) where {BRT}
  R = OO(X)
  R == OO(Y) || error("schemes can not be compared")
  return Spec(X)
end


function Base.intersect(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyQuo}
  ) where {BRT}
  R = OO(X)
  R === ambient_ring(Y) || error("schemes can not be compared")
  return Spec(Y)
end


function Base.intersect(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyLocalizedRing}
  ) where {BRT}
  R = OO(X)
  R === ambient_ring(Y) || error("schemes can not be compared")
  return Spec(Y)
end


function Base.intersect(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing}
  ) where {BRT}
  R = OO(X)
  R === ambient_ring(Y) || error("schemes can not be compared")
  return Spec(Y)
end


function Base.intersect(
    Y::AbsSpec{BRT, <:Any},
    X::AbsSpec{BRT, <:MPolyRing}
  ) where {BRT}
  return intersect(X, Y)
end


### For Specs of MPolyQuos
function Base.intersect(
    X::AbsSpec{BRT, <:MPolyQuo},
    Y::AbsSpec{BRT, <:MPolyQuo}
  ) where {BRT}
  R = ambient_ring(X)
  R === ambient_ring(Y) || error("schemes can not be compared")
  return Spec(quo(R, modulus(OO(X)) + modulus(OO(Y)))[1])
end


function Base.intersect(
    X::AbsSpec{BRT, <:MPolyQuo},
    Y::AbsSpec{BRT, <:MPolyLocalizedRing}
  ) where {BRT}
  R = ambient_ring(X)
  R === ambient_ring(Y) || error("schemes can not be compared")
  return Spec(quo(OO(Y), OO(Y)(modulus(OO(X))))[1])
end


function Base.intersect(
    X::AbsSpec{BRT, <:MPolyQuo},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing}
  ) where {BRT}
  R = ambient_ring(X)
  R === ambient_ring(Y) || error("schemes can not be compared")
  return Spec(quo(OO(Y), OO(Y)(modulus(OO(X))))[1])
end


function Base.intersect(
    Y::AbsSpec{BRT, <:Any},
    X::AbsSpec{BRT, <:MPolyQuo}
  ) where {BRT}
  return intersect(X, Y)
end


### For Specs of MPolyLocalizedRings
function Base.intersect(
    X::AbsSpec{BRT, <:MPolyLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyLocalizedRing}
  ) where {BRT}
  R = ambient_ring(X)
  R === ambient_ring(Y) || error("schemes can not be compared")
  return Spec(Localization(R, inverted_set(OO(X)) * inverted_set(OO(Y)))[1])
end


function Base.intersect(
    X::AbsSpec{BRT, <:MPolyLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing}
  ) where {BRT}
  R = ambient_ring(X)
  R === ambient_ring(Y) || error("schemes can not be compared")
  return Spec(R, modulus(quotient_ring(OO(Y))), inverted_set(OO(X))*inverted_set(OO(Y)))
end


function Base.intersect(
    Y::AbsSpec{BRT, <:Any},
    X::AbsSpec{BRT, <:MPolyLocalizedRing}
  ) where {BRT}
  return intersect(X, Y)
end


### For Specs of MPolyQuoLocalizedRings
function Base.intersect(
    X::AbsSpec{BRT, <:MPolyQuoLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing}
  ) where {BRT}
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

### compute the closure of X in Y
@Markdown.doc """
    closure(X::AbsSpec, Y::AbsSpec) 

Returns the closure of ``X`` in ``Y``.
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
