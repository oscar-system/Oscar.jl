export AbsMultSet
export AbsLocalizedRing
export ambient_ring, inverted_set
export reduce_fraction
export Localization

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
ambient_ring(S::AbsMultSet)

Returns the ambient ring `R` for a multiplicatively closed set `S ⊂ R`.
"""
function ambient_ring(S::AbsMultSet)
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
base_ring(W::AbsLocalizedRing)

Returns the base ring R for a localized ring of the form W = R[S⁻¹].
"""
function base_ring(W::AbsLocalizedRing)
  error("`base_ring` is not implemented for localized rings of type $(typeof(W))")
end

@Markdown.doc """
inverted_set(W::AbsLocalizedRing)

Returns the set S of at which has been localized for a localized ring W = R[S⁻¹].
"""
function inverted_set(W::AbsLocalizedRing)
  error("`inverted_set` is not implemented for localized rings of type $(typeof(W))")
end

### required functionality
@Markdown.doc """
Localization(S::AbsMultSet)

Returns the localization of the `ambient_ring` of `S` at `S` itself.
"""
function Localization(S::AbsMultSet)
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
(W::AbsLocalizedRing)(a::Oscar.IntegerUnion) = W(base_ring(W)(a))


#################################################################################
# Elements of localized rings                                                   #
#################################################################################
@Markdown.doc """
    AbsLocalizedRingElem{RingType, RingElemType, MultSetType}

The abstract type of an element of the localization R[S⁻¹] of a commutative ring 
R of type `RingType` with elements of type `RingElemType` at a multiplicatively 
closed set S of type `MultSetType`.
"""
abstract type AbsLocalizedRingElem{
    RingType <: AbstractAlgebra.Ring, 
    RingElemType <: AbstractAlgebra.RingElem, 
    MultSetType
  } <: AbstractAlgebra.RingElem end

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
  if needs_parentheses(f)
    print(io, "($(numerator(f)))//($(denominator(f)))")
  else 
    print(io, "$(numerator(f))//$(denominator(f))")
  end
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

function -(a::T) where {T<:AbsLocalizedRingElem}
  return (parent(a))(-numerator(a), denominator(a))
end

function *(a::T, b::T) where {T<:AbsLocalizedRingElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  return reduce_fraction((parent(a))(numerator(a)*numerator(b), denominator(a)*denominator(b)))
end

function *(a::RET, b::AbsLocalizedRingElem{RT, RET, MST}) where {RT, RET <: RingElem, MST}
  return reduce_fraction((parent(b))(a*numerator(b), denominator(b)))
end

function *(a::AbsLocalizedRingElem{RT, RET, MST}, b::RET) where {RT, RET <: RingElem, MST}
  return b*a
end

function Base.:(//)(a::Oscar.IntegerUnion, b::AbsLocalizedRingElem)
  error("function `//` not implemented for elements of type $(typeof(b))")
end

### Why are divisions not implemented as generic functions?
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

function ^(a::AbsLocalizedRingElem, i::fmpz)
  return parent(a)(numerator(a)^i, denominator(a)^i)
end

function ^(a::AbsLocalizedRingElem, i::Integer)
  return parent(a)(numerator(a)^i, denominator(a)^i)
end

function divexact(a::T, b::T; check::Bool=false) where {T<:AbsLocalizedRingElem} 
  error("method `divexact` not implemented for arguments of type $(typeof(a))")
end

function inv(a::AbsLocalizedRingElem) 
  return parent(a)(denominator(a), numerator(a))
end

########################################################################
# generic functions to adhere to the Oscar ring interface              #
########################################################################

isone(a::AbsLocalizedRingElem) = (numerator(a) == denominator(a))

isunit(f::AbsLocalizedRingElem) = numerator(f) in inverted_set(parent(f))

isdomain_type(T::Type{U}) where {U<:AbsLocalizedRingElem} = false # default set to false

isexact_type(T::Type{U}) where {U<:AbsLocalizedRingElem} = false # default set to false

function Base.hash(f::T, h::UInt) where {T<:AbsLocalizedRingElem} 
  r = 0x78a97cd90
  r = xor(r, hash(numerator(f), h))
  return xor(r, hash(denominator(f), h))
end

Base.deepcopy_internal(f::T, dict::IdDict) where {T<:AbsLocalizedRingElem} = parent(f)(copy(numerator(f)), copy(denominator(f)))

one(W::AbsLocalizedRing) = W(one(base_ring(W)))

zero(W::AbsLocalizedRing) = W(zero(base_ring(W)))

canonical_unit(f::LocRingElemType) where {LocRingElemType<:AbsLocalizedRingElem} = one(parent(f))

characteristic(W::AbsLocalizedRing) = characteristic(base_ring(W))

function Base.show(io::IO, W::AbsLocalizedRing) 
  print(io, "localization of ") 
  print(io, base_ring(W))
  print(io, " at the ")
  print(io, inverted_set(W))
end

needs_parentheses(f::AbsLocalizedRingElem) = true

function zero!(a::AbsLocalizedRingElem) 
  a = zero(parent(a))
  return a
end

function mul!(c::T, a::T, b::T) where {T<:AbsLocalizedRingElem} 
  c = a*b
  return c
end

function add!(c::T, a::T, b::T) where {T<:AbsLocalizedRingElem} 
  c = a+b
  return c
end

function addeq!(a::T, b::T) where {T<:AbsLocalizedRingElem}
  a = a+b
  return a
end

### promotion rules
promote_rule(::Type{AbsLocalizedRingElem{RT, RET, MST}}, ::Type{AbsLocalizedRingElem{RT, RET, MST}}) where {RT<:Ring, RET<:RingElement, MST} = AbsLocalizedRingElem{RT, RET, MST}

function promote_rule(::Type{AbsLocalizedRingElem{RT, RET, MST}}, ::Type{T}) where {RT<:Ring, RET<:RingElement, MST, T<:RingElement} 
  promote_rule(RET, T) ? AbsLocalizedRingElem{RT, RET, MST} : Union{}
end





### Needs to be overwritten in case of zero divisors!
iszero(a::AbsLocalizedRingElem) = iszero(numerator(a))


############################################################################
# Finitely generated ideals in localized rings                             #
############################################################################

@Markdown.doc """
AbsLocalizedIdeal{RingType, RingElemType, MultSetType}

Abstract type for finitely generated ideals ``I ⊂ R[S⁻¹]`` in localized rings. 
"""
abstract type AbsLocalizedIdeal{RingType, RingElemType, MultSetType} <: Ideal{RingElemType} end

### required getter functions
#Returns a Vector of generators of `I`.
function gens(I::AbsLocalizedIdeal)
  error("`gens(I)` has not been implemented for `I` of type $(typeof(I))")
end

# Returns the localized ring over which `I` is defined.
function base_ring(I::AbsLocalizedIdeal)
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
# Checks for ideal membership of `f` in `I`.
function Base.in(
    f::AbsLocalizedRingElem{RingType, RingElemType, MultSetType}, 
    I::AbsLocalizedIdeal{RingType, RingElemType, MultSetType}
  ) where {RingType, RingElemType, MultSetType}
  error("`in(f, I)` has not been implemented for `f` of type $(typeof(f)) and `I` of type $(typeof(I))")
end

# Checks for ideal membership of `f` in `I`.
function Base.in(
    f::RingElemType, 
    I::AbsLocalizedIdeal{RingType, RingElemType, MultSetType}
  ) where {RingType, RingElemType, MultSetType}
  error("`in(f, I)` has not been implemented for `f` of type $(typeof(f)) and `I` of type $(typeof(I))")
end


### A catchall implementation for the ideal arithmetic 
# Returns the product of the ideals `I` and `J`.
function Base.:*(I::T, J::T) where {T<:AbsLocalizedIdeal}
  W = base_ring(I) 
  W == base_ring(J) || error("the given ideals do not belong to the same ring")
  new_gens = [ f*g for f in gens(I) for g in gens(J)]
  return ideal(W, new_gens)
end

# Returns the sum of the ideals `I` and `J`.
function Base.:+(I::T, J::T) where {T<:AbsLocalizedIdeal}
  W = base_ring(I) 
  W == base_ring(J) || error("the given ideals do not belong to the same ring")
  return ideal(W, vcat(gens(I), gens(J)))
end


########################################################################
# Homomorphisms of localized rings                                     #
########################################################################

@Markdown.doc """
AbsLocalizedRingHom{RingType, RingElemType, DomainMultSetType, CodomainMultSetType} 

Homomorphism ``ϕ : P[U⁻¹] → Q[V⁻¹]`` of localized rings with ``P`` and 
``Q`` of type `RingType`, ``U`` of type `DomainMultSetType`, and 
``V`` of type `CodomainMultSetType`.
"""
abstract type AbsLocalizedRingHom{
    RingType, 
    RingElemType, 
    DomainMultSetType, 
    CodomainMultSetType
  } <: Map{
    AbsLocalizedRing{RingType, RingElemType, DomainMultSetType}, 
    AbsLocalizedRing{RingType, RingElemType, CodomainMultSetType}, 
    SetMap,
    AbsLocalizedRingHom
  }
end


### required getter functions
@Markdown.doc """
domain(f::AbsLocalizedRingHom) 

Returns the domain of definition of `f`.
"""
function domain(f::AbsLocalizedRingHom) 
  error("`domain(f)` not implemented for `f` of type $(typeof(f))")
end

@Markdown.doc """
codomain(f::AbsLocalizedRingHom) 

Returns the codomain of `f`.
"""
function codomain(f::AbsLocalizedRingHom) 
  error("`codomain(f)` not implemented for `f` of type $(typeof(f))")
end

### required functionality
@Markdown.doc """
(f::AbsLocalizedRingHom)(a::T) where {T<:RingElement}

Applies the map `f` to the element `a` in the domain of `f`.
"""
(f::AbsLocalizedRingHom{RT, RET, DMST, CMST})(a::AbsLocalizedRingElem{RT, RET, DMST}) where {RT, RET, DMST, CMST} = error("mapping of elements of type $(typeof(a)) not implemented")

### generic functions
(f::AbsLocalizedRingHom{RT, RET, DMST, CMST})(a::RET) where {RT, RET, DMST, CMST} = f(domain(f)(a))
(f::AbsLocalizedRingHom)(a::Integer) = f(domain(f)(a))
(f::AbsLocalizedRingHom)(a::fmpz) = f(domain(f)(a))

@Markdown.doc """
(f::AbsLocalizedRingHom)(I::Ideal)

Returns the ideal generated by the images `f(hᵢ)` of the generators `hᵢ` of `I`.
"""
(f::AbsLocalizedRingHom{RT, RET, DMST, CMST})(I::AbsLocalizedIdeal{RT, RET, DMST}) where {RT, RET, DMST, CMST} = ideal(codomain(f), f.(gens(I)))

(f::AbsLocalizedRingHom{RT, RET, DMST, CMST})(I::Ideal{RET}) where {RT, RET, DMST, CMST} = ideal(codomain(f), f.(domain(f).(gens(I))))

### implementing the Oscar map interface
check_composable(
    f::AbsLocalizedRingHom{RT, RET, MST1, MST2},
    g::AbsLocalizedRingHom{RT, RET, MST2, MST3}
  ) where {RT, RET, MST1, MST2, MST3} = (codomain(f) == domain(g))

function Base.show(io::IO, f::AbsLocalizedRingHom{RT, RET, MST1, MST2}) where {RT, RET, MST1, MST2}
  print(io, "morphism of localized rings from $(domain(f)) to $(codomain(f))")
end
