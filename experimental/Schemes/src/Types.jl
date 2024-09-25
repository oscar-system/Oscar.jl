export AbsProjectiveScheme
export EmptyScheme
export ProjectiveScheme
export ProjectiveSchemeMor
export VarietyFunctionField
export VarietyFunctionFieldElem


#####################################################################################################
# Desingularization morphism: birational map between covered schemes with smooth domain
#####################################################################################################
@doc raw"""
    MixedBlowUpSequence{
              DomainType<:AbsCoveredScheme,
              CodomainType<:AbsCoveredScheme
            }<:AbsDesingMor{  DomainType,
                              CodomainType, 
                              MixedBlowUpSequence{DomainType, CodomainType}
                                              }
A datastructure to encode sequences of blow-ups and normalizations of covered schemes
as needed for desingularization of non-embedded schemes by the approaches of Zariski and of
Lipman. 
"""

@attributes mutable struct MixedBlowUpSequence{
                                                DomainType<:AbsCoveredScheme,
                                                CodomainType<:AbsCoveredScheme
                                              }<:AbsDesingMor{  DomainType,
                                                                CodomainType, 
                                                                MixedBlowUpSequence{DomainType, CodomainType}
                                              }
  maps::Vector{Union{<:BlowupMorphism,<:NormalizationMorphism}}       # count right to left:
                                                 # original scheme is codomain of map 1
  # boolean flags
  resolves_sing::Bool                            # domain not smooth yet?
  is_trivial::Bool                               # codomain already smooth?
  is_strong::Bool                                # snc divisors ensured?

  # fields for caching, to be filled during desingularization
  # always carried along to domain(maps[end])) using strict_transform
  ex_div::Vector{AbsIdealSheaf}      # list of exc. divisors arising from individual steps
  dont_meet::Vector{Tuple{Int,Int}}             # mostly for dim=2: intersections which cannot exist
                                                # according to intermediate computations
  caution_multi_charts::Vector{Tuple{Int,Int}}  # only for dim=2: intersection of divisors not
                                                # entirely visible in a single chart

  # keep track of the normalization steps
  normalization_steps::Vector{Int}

  # fields for caching to be filled a posteriori (on demand, only if partial_res==false)
  underlying_morphism::CompositeCoveredSchemeMorphism{DomainType, CodomainType}
  exceptional_divisor::AbsWeilDivisor
  exceptional_locus::AbsAlgebraicCycle

  function MixedBlowUpSequence(maps::Vector{<:AbsCoveredSchemeMorphism})
    n = length(maps)
    for i in 1:n
      @assert all(x->((x isa BlowupMorphism) || (x isa NormalizationMorphism)), maps) "only blow-ups and normalizations allowed"
    end
    for i in 1:n-1
      @assert domain(maps[i]) === codomain(maps[i+1]) "not a sequence of morphisms"
    end
    resi = new{typeof(domain(maps[end])),typeof(codomain(first(maps)))}(maps)
    resi.normalization_steps = [i for i in 1:n if maps[i] isa NormalizationMorphism]
    return resi
  end

end
