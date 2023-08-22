###########################################
# (1) Types for schemes
###########################################

@doc raw"""
    Scheme{BaseRingType<:Ring} 

A scheme over a ring ``𝕜`` of type `BaseRingType`.
"""
abstract type Scheme{BaseRingType} end

@doc raw"""
    AbsSpec{BaseRingType, RingType<:Ring}

An affine scheme ``X = Spec(R)`` with ``R`` of type `RingType` over
a ring ``𝕜`` of type `BaseRingType`.
"""
abstract type AbsSpec{BaseRingType, RingType<:Ring} <: Scheme{BaseRingType} end

@doc raw"""
    AbsCoveredScheme{BaseRingType}

An abstract scheme ``X`` over some `base_ring` ``𝕜`` of type 
`BaseRingType`, given by means of affine charts and their glueings.
"""
abstract type AbsCoveredScheme{BaseRingType} <: Scheme{BaseRingType} end


###########################################
# (2) Toric Varieties
###########################################

@attributes mutable struct NormalToricVariety <: AbsCoveredScheme{QQField}
           polymakeNTV::Polymake.BigObject
           NormalToricVariety(polymakeNTV::Polymake.BigObject) = new(polymakeNTV)
end

@attributes mutable struct AffineNormalToricVariety <: AbsSpec{QQField, MPolyQuoRing}
           polymakeNTV::Polymake.BigObject
           AffineNormalToricVariety(polymakeNTV::Polymake.BigObject) = new(polymakeNTV)
end

@attributes mutable struct CyclicQuotientSingularity
    polymakeNTV::Polymake.BigObject
end

NormalToricVarietyType = Union{NormalToricVariety, AffineNormalToricVariety, CyclicQuotientSingularity}
