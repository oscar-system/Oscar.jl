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

A scheme ``X`` over some `base_ring` ``𝕜`` of type 
`BaseRingType`, given by means of affine charts and their gluings.
"""
abstract type AbsCoveredScheme{BaseRingType} <: Scheme{BaseRingType} end

@doc raw"""
    AbsAlgebraicCycle{CoveredSchemeType<:AbsCoveredScheme, CoefficientRingType<:Ring}
    
An algebraic cycle $D$ on a (locally) Noetherian integral scheme $X$ with coefficients in a ring $R$ 
is a formal linear combination $\sum_i a_i D_i$ with 
$D_i \subseteq X$ integral, closed subschemes and the $a_i \in R$.

Such a cycle is represented non-uniquely as a formal sum $E = \sum_l b_l \mathcal{I}_l$ 
of equidimensional ideal sheaves $\mathcal{I}_l \subseteq \mathcal{O}_X$. 
For an equidimensional ideal sheaf $\mathcal{I}$ its interpretation as a cycle is as follows:
Let $V(\mathcal{I})=E_{1} \cup \dots E_{n}$ be the decomposition of the vanishing locus
of $\mathcal{I}$ into irreducible components $E_i=V(\mathcal{P}_i)$ with $\mathcal{P}_i$ prime. 
Then $E$ corresponds to the cycle $D = \sum_{i=1}^{n} \mathrm{colength}_{\mathcal{P}_i}(\mathcal{I})E_i$.
"""
abstract type AbsAlgebraicCycle{
                                CoveredSchemeType<:AbsCoveredScheme, 
                                CoefficientRingType<:AbstractAlgebra.Ring
                               }
end

@doc raw"""
    AbsWeilDivisor{CoveredSchemeType, CoefficientRingType} <: AbsAlgebraicCycle{CoveredSchemeType, CoefficientRingType}

A Weil divisor with coefficients of type `CoefficientRingType` on a (locally) Noetherian integral scheme ``X``  of type `CoveredSchemeType`.
"""
abstract type AbsWeilDivisor{CoveredSchemeType, CoefficientRingType} <: AbsAlgebraicCycle{CoveredSchemeType, CoefficientRingType} end

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
