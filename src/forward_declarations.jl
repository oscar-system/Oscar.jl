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
Let $V(\mathcal{I})=E_{1} \cup \dots \cup E_{n}$ be the decomposition of the vanishing locus
of $\mathcal{I}$ into irreducible components $E_i=V(\mathcal{P}_i)$ with $\mathcal{P}_i$ prime. 
Then $E$ corresponds to the cycle $D = \sum_{i=1}^{n} \mathrm{colength}_{\mathcal{P}_i}(\mathcal{I})E_i$.

# Examples
```jldoctest
julia> P2 = projective_space(QQ,2); (s0,s1,s2) = homogeneous_coordinates(P2);

julia> I = ideal_sheaf(P2,ideal([s0,s1^2]))
Sheaf of ideals
  on scheme over QQ covered with 3 patches
    1: [(s1//s0), (s2//s0)]   affine 2-space
    2: [(s0//s1), (s2//s1)]   affine 2-space
    3: [(s0//s2), (s1//s2)]   affine 2-space
with restrictions
  1: Ideal (1, (s1//s0)^2)
  2: Ideal ((s0//s1), 1)
  3: Ideal ((s0//s2), (s1//s2)^2)

julia> D = algebraic_cycle(I)
Effective algebraic cycle
  on scheme over QQ covered with 3 patches
with coefficients in integer ring
given as the formal sum of
  1 * sheaf of ideals

julia> irreducible_decomposition(D)
Effective algebraic cycle
  on scheme over QQ covered with 3 patches
with coefficients in integer ring
given as the formal sum of
  2 * sheaf of prime ideals

```
"""
abstract type AbsAlgebraicCycle{
                                CoveredSchemeType<:AbsCoveredScheme, 
                                CoefficientRingType<:AbstractAlgebra.Ring
                               }
end

@doc raw"""
    AbsWeilDivisor{CoveredSchemeType, CoefficientRingType} <: AbsAlgebraicCycle{CoveredSchemeType, CoefficientRingType}

A Weil divisor with coefficients of type `CoefficientRingType` on a (locally) Noetherian integral scheme ``X``  of type `CoveredSchemeType`.

# Examples
```
julia> P2 = projective_space(QQ,2); (s0,s1,s2) = homogeneous_coordinates(P2);

julia> I = ideal((s0*s1)^2);

julia> II = ideal_sheaf(P2, I);

julia> D = weil_divisor(II)
Effective weil divisor
  on scheme over QQ covered with 3 patches
with coefficients in integer ring
given as the formal sum of
  1 * sheaf of ideals

julia> E = irreducible_decomposition(D)
Effective weil divisor
  on scheme over QQ covered with 3 patches
with coefficients in integer ring
given as the formal sum of
  2 * prime ideal sheaf on scheme over QQ covered with 3 patches extended from ideal ((s1//s0)) on affine 2-space
  2 * prime ideal sheaf on scheme over QQ covered with 3 patches extended from ideal ((s0//s1)) on affine 2-space

julia> P = components(E)[1]
Prime ideal sheaf on Scheme over QQ covered with 3 patches extended from Ideal ((s1//s0)) on Affine 2-space

julia> components(D)[1] == II
true

julia> D[II] # to get the coefficient
1

julia> E[P]
2
```
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

