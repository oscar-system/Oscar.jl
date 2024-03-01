###########################################
# (1) Types for schemes
###########################################

@doc raw"""
    Scheme{BaseRingType<:Ring} 

A scheme over a ring ``𝕜`` of type `BaseRingType`.
"""
abstract type Scheme{BaseRingType} end

@doc raw"""
    AbsAffineScheme{BaseRingType, RingType<:Ring}

An affine scheme ``X = Spec(R)`` with ``R`` of type `RingType` over
a ring ``𝕜`` of type `BaseRingType`.
"""
abstract type AbsAffineScheme{BaseRingType, RingType<:Ring} <: Scheme{BaseRingType} end

@doc raw"""
    AbsCoveredScheme{BaseRingType}

An abstract scheme ``X`` over some `base_ring` ``𝕜`` of type 
`BaseRingType`, given by means of affine charts and their gluings.
"""
abstract type AbsCoveredScheme{BaseRingType} <: Scheme{BaseRingType} end

### Algebraic cycles and divisors 
# CAUTION: This has been moved here from experimental!!!
abstract type AbsAlgebraicCycle{
                                CoveredSchemeType<:AbsCoveredScheme, 
                                CoefficientRingType<:AbstractAlgebra.Ring
                               }
end

abstract type AbsWeilDivisor{CoveredSchemeType, CoefficientRingType} <: AbsAlgebraicCycle{CoveredSchemeType, CoefficientRingType} end

### END Algebraic cycles and divisors

###########################################
# (2) Toric Varieties
###########################################

@attributes mutable struct NormalToricVariety <: AbsCoveredScheme{QQField}
           polymakeNTV::Polymake.BigObject
           NormalToricVariety(polymakeNTV::Polymake.BigObject) = new(polymakeNTV)
end

@attributes mutable struct AffineNormalToricVariety <: AbsAffineScheme{QQField, MPolyQuoRing}
           polymakeNTV::Polymake.BigObject
           AffineNormalToricVariety(polymakeNTV::Polymake.BigObject) = new(polymakeNTV)
end

@attributes mutable struct CyclicQuotientSingularity
    polymakeNTV::Polymake.BigObject
end

const NormalToricVarietyType = Union{NormalToricVariety, AffineNormalToricVariety, CyclicQuotientSingularity}
