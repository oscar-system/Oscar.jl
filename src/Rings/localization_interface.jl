export AbsMultSet
export AbsLocalizedRing
export original_ring, inverted_set
export reduce_fraction
export localize_at

export AbsLocalizedRingElem
export numerator, denominator, parent

export AbsLocalizedIdeal
export ideal

import AbstractAlgebra.Ring

#################################################################################
# General framework for localizations of rings; to be used with affine algebras #
#################################################################################

#################################################################################
# Multiplicatively closed sets of (commutative) rings                           #
#################################################################################

@Markdown.doc """
    AbsMultSet{RingType, RingElemType}

Abstract type for a multiplicatively closed set in a commutative (Noetherian) ring 
R of type `RingType` with elements of type `RingElemType`.
"""
abstract type AbsMultSet{RingType<:Ring, RingElemType<:RingElem} end

### required getter functions
@Markdown.doc """
    ambient_ring(S::AbsMultSet{RingType, RingElemType}) where {RingType, RingElemType}

Returns the ambient ring `R` for a multiplicatively closed set `S ⊂ R`.
"""
function ambient_ring(S::AbsMultSet{RingType, RingElemType}) where {RingType, RingElemType}
  error("method `ambient_ring` not implemented for multiplicatively closed sets of type $(typeof(S))")
end

### required functionality
@Markdown.doc """
    in(f::RingElemType, S::AbsMultSet{RingType, RingElemType}) where {RingType, RingElemType}

Returns `true` if `f` belongs to `S`; `false` otherwise.

**Note:** If this routine is not implemented, the function call will default to the 
execution of an error message. 
"""
function Base.in(f::RingElemType, S::AbsMultSet{RingType, RingElemType}) where {RingType, RingElemType}
  error("method not implemented for multiplicatively closed sets of type $(typeof(S))")
end


#################################################################################
# Localizations of (commutative) rings at multiplicatively closed sets          #
#################################################################################
@Markdown.doc """
    AbsLocalizedRing{RingType, RingElemType, MultSetType}

The localization R[S⁻¹] of a ring R of type `RingType` with elements of type `RingElemType` at a 
multiplicatively closed set S of type `MultSetType`. 

In general, the arithmetic of such a localized ring R[S⁻¹] should be implemented using fractions 
of elements in the original ring R. The methods provided for the multiplicatively closed set S 
can be used to check whether a given denominator is admissible for the specific localization. 

Depending on the actual type of R and S, further functionality can then be provided using 
different Groebner-basis driven backends. 
"""
abstract type AbsLocalizedRing{RingType, RingElemType, MultSetType} <: Ring end

### required getter functions
@Markdown.doc """
    original_ring(W::AbsLocalizedRing{RingType, RingElemType, MultSetType}) where {RingType, RingElemType, MultSetType} 

Returns the original ring R for a localized ring of the form W = R[S⁻¹].
"""
function original_ring(W::AbsLocalizedRing{RingType, RingElemType, MultSetType}) where {RingType, RingElemType, MultSetType} 
  error("`original_ring` is not implemented for localized rings of type $(typeof(W))")
end

@Markdown.doc """
    inverted_set(W::AbsLocalizedRing{RingType, RingElemType, MultSetType}) where {RingType, RingElemType, MultSetType}

Returns the set S of at which has been localized for a localized ring W = R[S⁻¹].
"""
function inverted_set(W::AbsLocalizedRing{RingType, RingElemType, MultSetType}) where {RingType, RingElemType, MultSetType}
  error("`inverted_set` is not implemented for localized rings of type $(typeof(W))")
end

### required functionality
@Markdown.doc """
    localize_at(S::AbsMultSet{RingType, RingElemType}) where {RingType, RingElemType}

Returns the localization of the `ambient_ring` of `S` at `S` itself.
"""
function localize_at(S::AbsMultSet{RingType, RingElemType}) where {RingType, RingElemType}
  error("localizations at multiplicatively closed sets of type $(typeof(S)) are not implemented")
end

@Markdown.doc """
    (W::AbsLocalizedRing{RingType, RingElemType, MultSetType})(f::AbstractAlgebra.Generic.Frac{RingElemType}) where {RingType, RingElemType, MultSetType} 

Converts a fraction f = a//b to an element of the localized ring W.
"""
function (W::AbsLocalizedRing{RingType, RingElemType, MultSetType})(f::AbstractAlgebra.Generic.Frac{RingElemType}) where {RingType, RingElemType, MultSetType} 
  error("conversion for fractions to elements of type $(typeof(W)) is not implemented")
end

### required conversions
@Markdown.doc """
    (W::AbsLocalizedRing{RingType, RingElemType, MultSetType})(a::RingElemType) where {RingType, RingElemType, MultSetType} 

Converts an element `a` to an element of `W`.
"""
function (W::AbsLocalizedRing{RingType, RingElemType, MultSetType})(a::RingElemType) where {RingType, RingElemType, MultSetType} 
  error("conversion of elements of type $(RingElemType) to elements of $(typeof(W)) is not implemented")
end

@Markdown.doc """
    (W::AbsLocalizedRing{RingType, RingElemType, MultSetType})(a::RingElemType, b::RingElemType) where {RingType, RingElemType, MultSetType} 

Converts a pair `(a, b)` to a fraction `a/b` in `W`.
"""
function (W::AbsLocalizedRing{RingType, RingElemType, MultSetType})(a::RingElemType, b::RingElemType) where {RingType, RingElemType, MultSetType} 
  error("conversion of pairs `(a, b)` of elements of type $(RingElemType) to fractions `a/b` in a ring of type $(typeof(W)) is not implemented")
end

### Other conversions for the sake of convenience
(W::AbsLocalizedRing{RingType, RingElemType, MultSetType})(a::Oscar.IntegerUnion) where {RingType, RingElemType, MultSetType} = W(original_ring(W)(a))


#################################################################################
# Elements of localized rings                                                   #
#################################################################################
@Markdown.doc """
    AbsLocalizedRingElem{RingType, RingElemType, MultSetType}

The abstract type of an element of the localization R[S⁻¹] of a commutative ring 
R of type `RingType` with elements of type `RingElemType` at a multiplicatively 
closed set S of type `MultSetType`.
"""
abstract type AbsLocalizedRingElem{RingType, RingElemType, MultSetType} end

### required getter functions 
@Markdown.doc """
numerator(f::AbsLocalizedRingElem)

Returns the numerator of `f`.
"""
function numerator(f::AbsLocalizedRingElem)
  error("`numerator` is not implemented for elements of type $(typeof(f))")
end

@Markdown.doc """
denominator(f::AbsLocalizedRingElem)

Returns the denominator of `f`.
"""
function denominator(f::AbsLocalizedRingElem)
  error("`denominator` is not implemented for elements of type $(typeof(f))")
end

@Markdown.doc """
parent(f::AbsLocalizedRingElem)

Returns the parent ring R[S⁻¹] of `f`.
"""
function parent(f::AbsLocalizedRingElem)
  error("`parent` is not implemented for the type $(typeof(f))")
end

### default functionality for printing
function Base.show(io::IO, f::AbsLocalizedRingElem)
  print(io, "($(numerator(f)))//($(denominator(f)))")
end


########################################################################
# Arithmetic; a dumb catchall implementation, NOT performant!          #
########################################################################

@Markdown.doc """
reduce_fraction(a::AbsLocalizedRingElem)

Reduce the fraction a = p/q. **Warning**: The catchall-implementation does nothing!
"""
function reduce_fraction(a::AbsLocalizedRingElem)
  return a
end

function +(a::T, b::T) where {T<:AbsLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  if denominator(a) == denominator(b) 
    return reduce_fraction((parent(a))(numerator(a) + numerator(b), denominator(a)))
  end
  return reduce_fraction((parent(a))(numerator(a)*denominator(b) + numerator(b)*denominator(a), denominator(a)*denominator(b)))
end

function -(a::T, b::T) where {T<:AbsLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  if denominator(a) == denominator(b) 
    return reduce_fraction((parent(a))(numerator(a) - numerator(b), denominator(a)))
  end
  return reduce_fraction((parent(a))(numerator(a)*denominator(b) - numerator(b)*denominator(a), denominator(a)*denominator(b)))
end

function *(a::T, b::T) where {T<:AbsLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return reduce_fraction((parent(a))(numerator(a)*numerator(b), denominator(a)*denominator(b)))
end

function Base.:(//)(a::Oscar.IntegerUnion, b::AbsLocalizedRingElem)
  error("function `//` not implemented for elements of type $(typeof(b))")
end

### Why are divisions not implemented on the default level?
# Say both a = p⋅c/q and b = f⋅c/g are elements of a localized 
# ring R[S⁻¹] such that f⋅c ∉ S, but f ∈ S. Then by cancellation 
# a/b is an admissible element of the localization. But this 
# cancellation can not be implemented on a generic level, but 
# only for the specific ring R.
function Base.:(//)(a::T, b::T) where {T<:AbsLocalizedRingElem}
  error("function `//` not implemented for elements of type $(typeof(b))")
end

function ==(a::T, b::T) where {T<:AbsLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return numerator(a)*denominator(b) == numerator(b)*denominator(a)
end

function ^(a::AbsLocalizedRingElem, i::Oscar.IntegerUnion)
  return parent(a)(numerator(a)^i, denominator(a)^i)
end

isone(a::AbsLocalizedRingElem) = (numerator(a) == denominator(a))

### Needs to be overwritten in case of zero divisors!
iszero(a::AbsLocalizedRingElem) = iszero(numerator(a))


############################################################################
# Finitely generated ideals in localized rings                             #
############################################################################

@Markdown.doc """
    AbsLocalizedIdeal{RingType, RingElemType, MultSetType} <: Ideal{RingElemType}

Abstract type for finitely generated ideals ``IS⁻¹ ⊂ R[S⁻¹]`` in localized rings. 
"""
abstract type AbsLocalizedIdeal{RingType, RingElemType, MultSetType} <: Ideal{RingElemType} end

### required getter functions
function gens(I::AbsLocalizedIdeal{RingType, RingElemType, MultSetType}) where {RingType, RingElemType, MultSetType}
  error("`gens(I)` has not been implemented for `I` of type $(typeof(I))")
end

function base_ring(I::AbsLocalizedIdeal{RingType, RingElemType, MultSetType}) where {RingType, RingElemType, MultSetType}
  error("`base_ring(I)` has not been implemented for `I` of type $(typeof(I))")
end

### required constructors
function ideal(
    W::AbsLocalizedRing{RingType, RingElemType, MultSetType}, 
    f::AbsLocalizedRingElem{RingType, RingElemType, MultSetType} 
  ) where {RingType, RingElemType, MultSetType}
  error("`ideal(W, f)` has not been implemented for `W` of type $(typeof(W)) and `f` of type $(typeof(f))")
end

function ideal(
    W::AbsLocalizedRing{RingType, RingElemType, MultSetType}, 
    v::Vector{AbsLocalizedRingElem{RingType, RingElemType, MultSetType}}
  ) where {RingType, RingElemType, MultSetType}
  error("`ideal(W, v)` has not been implemented for `W` of type $(typeof(W)) and `v` of type $(typeof(v))")
end

function ideal(
    W::AbsLocalizedRing{RingType, RingElemType, MultSetType}, 
    f::RingElemType 
  ) where {RingType, RingElemType, MultSetType}
  error("`ideal(W, f)` has not been implemented for `W` of type $(typeof(W)) and `f` of type $(typeof(f))")
end

function ideal(
    W::AbsLocalizedRing{RingType, RingElemType, MultSetType}, 
    v::Vector{RingElemType}
  ) where {RingType, RingElemType, MultSetType}
  error("`ideal(W, v)` has not been implemented for `W` of type $(typeof(W)) and `v` of type $(typeof(v))")
end

### required functionality
function Base.in(
    f::AbsLocalizedRingElem{RingType, RingElemType, MultSetType}, 
    I::AbsLocalizedIdeal{RingType, RingElemType, MultSetType}
  ) where {RingType, RingElemType, MultSetType}
  error("`in(f, I)` has not been implemented for `f` of type $(typeof(f)) and `I` of type $(typeof(I))")
end

function Base.in(
    f::RingElemType, 
    I::AbsLocalizedIdeal{RingType, RingElemType, MultSetType}
  ) where {RingType, RingElemType, MultSetType}
  error("`in(f, I)` has not been implemented for `f` of type $(typeof(f)) and `I` of type $(typeof(I))")
end


### A catchall implementation for the ideal arithmetic 
function Base.:*(I::T, J::T) where {T<:AbsLocalizedIdeal}
  W = base_ring(I) 
  W == base_ring(J) || error("the given ideals do not belong to the same ring")
  new_gens = [ f*g for f in gens(I) for g in gens(J)]
  return ideal(W, new_gens)
end

function Base.:+(I::T, J::T) where {T<:AbsLocalizedIdeal}
  W = base_ring(I) 
  W == base_ring(J) || error("the given ideals do not belong to the same ring")
  return ideal(W, vcat(gens(I), gens(J)))
end
