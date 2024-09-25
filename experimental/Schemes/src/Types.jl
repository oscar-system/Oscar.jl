export AbsProjectiveScheme
export EmptyScheme
export ProjectiveScheme
export ProjectiveSchemeMor
export VarietyFunctionField
export VarietyFunctionFieldElem

abstract type AbsProjectiveGluing{
                                   GluingType<:AbsGluing,
                                  }
end

@doc raw"""
  LazyProjectiveGluing(
      X::AbsProjectiveScheme,
      Y::AbsProjectiveScheme,
      BG::AbsGluing,
      compute_function::Function,
      gluing_data
    )

Produce a container `pg` to host a non-computed `ProjectiveGluing` of ``X`` with ``Y``.

The arguments consist of

 * the patches ``X`` and ``Y`` to be glued;
 * a gluing `BG` of the `base_scheme`s of ``X`` and ``Y`` over which the
   `ProjectiveGluing` to be computed sits;
 * a function `compute_function` which takes a single argument `gluing_data`
   of arbitrary type and actually carries out the computation;
 * an arbitrary struct `gluing_data` that the user can fill with whatever
   information is needed to properly feed their `compute_function`.

The container `pg` can then be stored as the gluing of ``X`` and ``Y``. As soon
as it is asked about any data on the gluing beyond its two `patches`, it will
invoke the internally stored `compute_function` to actually carry out the computation
of the gluing and then serve the incoming request on the basis of that result.
The latter actual `ProjectiveGluing` will then be cached.
"""
mutable struct LazyProjectiveGluing{
                                     GluingType<:AbsGluing,
                                     GluingDataType
                                    } <: AbsProjectiveGluing{GluingType}
  base_gluing::GluingType
  patches::Tuple{AbsProjectiveScheme, AbsProjectiveScheme}
  compute_function::Function
  gluing_data::GluingDataType
  underlying_gluing::AbsProjectiveGluing

  function LazyProjectiveGluing(
      X::AbsProjectiveScheme,
      Y::AbsProjectiveScheme,
      BG::AbsGluing,
      compute_function::Function,
      gluing_data
    )
    (base_scheme(X), base_scheme(Y)) == patches(BG) || error("gluing is incompatible with provided patches")
    return new{typeof(BG), typeof(gluing_data)}(BG, (X, Y), compute_function, gluing_data)
  end
end

@doc raw"""
    ProjectiveGluing(
        G::GluingType,
        incP::IncType, incQ::IncType,
        f::IsoType, g::IsoType;
        check::Bool=true
      ) where {GluingType<:AbsGluing, IncType<:ProjectiveSchemeMor, IsoType<:ProjectiveSchemeMor}

The `AbsProjectiveSchemeMorphism`s `incP` and `incQ` are open embeddings over open
embeddings of their respective `base_scheme`s.

        PX ↩ PU ≅ QV ↪ QY
      π ↓    ↓    ↓    ↓ π
    G : X  ↩ U  ≅ V  ↪ Y

This creates a gluing of the projective schemes `codomain(incP)` and `codomain(incQ)`
over a gluing `G` of their `base_scheme`s along the morphisms of `AbsProjectiveScheme`s
`f` and `g`, identifying `domain(incP)` and `domain(incQ)`, respectively.
"""
mutable struct ProjectiveGluing{
                                 GluingType<:AbsGluing,
                                 IsoType1<:ProjectiveSchemeMor,
                                 IncType1<:ProjectiveSchemeMor,
                                 IsoType2<:ProjectiveSchemeMor,
                                 IncType2<:ProjectiveSchemeMor,
                                } <: AbsProjectiveGluing{GluingType}
  G::GluingType # the underlying gluing of the base schemes
  inc_to_P::IncType1
  inc_to_Q::IncType2
  f::IsoType1
  g::IsoType2

  ###
  # Given two relative projective schemes and a gluing
  #
  #       PX ↩ PU ≅ QV ↪ QY
  #     π ↓    ↓    ↓    ↓ π
  #   G : X  ↩ U  ≅ V  ↪ Y
  #
  # this constructs the gluing of PX and QY along
  # their open subsets PU and QV, given the two inclusions
  # and isomorphisms over the gluing G in the base schemes.
  function ProjectiveGluing(
      G::GluingType,
      incP::IncType1, incQ::IncType2,
      f::IsoType1, g::IsoType2;
      check::Bool=true
    ) where {GluingType<:AbsGluing, IncType1<:ProjectiveSchemeMor,IncType2<:ProjectiveSchemeMor, IsoType1<:ProjectiveSchemeMor, IsoType2<:ProjectiveSchemeMor}
    (X, Y) = patches(G)
    (U, V) = gluing_domains(G)
    @vprint :Gluing 1 "computing projective gluing\n"
    @vprint :Gluing 2 "$(X), coordinates $(ambient_coordinates(X))\n"
    @vprint :Gluing 2 "and\n"
    @vprint :Gluing 2 "$(Y) coordinates $(ambient_coordinates(X))\n"
    (fb, gb) = gluing_morphisms(G)
    (PX, QY) = (codomain(incP), codomain(incQ))
    (PU, QV) = (domain(incP), domain(incQ))
    (base_scheme(PX) === X && base_scheme(QY) === Y) || error("base gluing is incompatible with the projective schemes")
    domain(f) === codomain(g) === PU && domain(g) === codomain(f) === QV || error("maps are not compatible")
    SPU = homogeneous_coordinate_ring(domain(f))
    SQV = homogeneous_coordinate_ring(codomain(f))
    @check begin
      # check the commutativity of the pullbacks
      all(y->(pullback(f)(SQV(OO(V)(y))) == SPU(pullback(fb)(OO(V)(y)))), gens(base_ring(OO(Y)))) || error("maps do not commute")
      all(x->(pullback(g)(SPU(OO(U)(x))) == SQV(pullback(gb)(OO(U)(x)))), gens(base_ring(OO(X)))) || error("maps do not commute")
      fc = map_on_affine_cones(f, check=false)
      gc = map_on_affine_cones(g, check=false)
      idCPU = compose(fc, gc)
      idCPU == identity_map(domain(fc)) || error("composition of maps is not the identity")
      idCQV = compose(gc, fc)
      idCQV == identity_map(domain(gc)) || error("composition of maps is not the identity")
      # idPU = compose(f, g)
      # all(t->(pullback(idPU)(t) == t), gens(SPU)) || error("composition of maps is not the identity")
      # idQV = compose(g, f)
      # all(t->(pullback(idQV)(t) == t), gens(SQV)) || error("composition of maps is not the identity")
    end
    @vprint :Gluing 1 "done computing the projective gluing\n"
    return new{GluingType, IsoType1, IncType1, IsoType2, IncType2}(G, incP, incQ, f, g)
  end
end

### Proper schemes π : Z → X over a covered base scheme X
#
# When {Uᵢ} is an affine covering of X, the datum stored
# consists of a list of projective schemes
#
#   Zᵢ ⊂ ℙʳ⁽ⁱ⁾(𝒪(Uᵢ)) → Uᵢ
#
# with varying ambient spaces ℙʳ⁽ⁱ⁾(𝒪(Uᵢ)) and a list of
# identifications (transitions)
#
#   Zᵢ ∩ π⁻¹(Uⱼ) ≅ Zⱼ ∩ π⁻¹(Uᵢ)
#
# of projective schemes over Uᵢ∩ Uⱼ for all pairs (i,j).
#
# These structs are designed to accommodate blowups of
# covered schemes along arbitrary centers, as well as
# projective bundles.

@attributes mutable struct CoveredProjectiveScheme{BRT} <: Scheme{BRT}
  Y::AbsCoveredScheme # the base scheme
  BC::Covering # the reference covering of the base scheme
  patches::IdDict{AbsAffineScheme, AbsProjectiveScheme} # the projective spaces over the affine patches in the base covering
  gluings::IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, AbsProjectiveGluing} # the transitions sitting over the affine patches in the gluing domains of the base scheme

  function CoveredProjectiveScheme(
      Y::AbsCoveredScheme,
      C::Covering,
      projective_patches::IdDict{AbsAffineScheme, AbsProjectiveScheme},
      projective_gluings::IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, AbsProjectiveGluing};
      check::Bool=true
    )
    C in coverings(Y) || error("covering not listed")
    for P in values(projective_patches)
      any(x->x===base_scheme(P), patches(C)) || error("base scheme not found in covering")
    end
    for (U, V) in keys(gluings(C))
      (U, V) in keys(projective_gluings) || error("not all projective gluings were provided")
    end
    return new{base_ring_type(Y)}(Y, C, projective_patches, projective_gluings)
  end
end


@doc raw"""
    AbsDesingMor{
                                  DomainType<:AbsCoveredScheme,
                                  CodomainType<:AbsCoveredScheme,
                                  BlowupMorphismType
       } <: AbsCoveredSchemeMorphism{
                                 DomainType,
                                 CodomainType,
                                 Nothing,
                                 BlowupMorphismType
                                }
Abstract type for desingularizations ``f : X -> Y `` of schemes where

  * ``Y`` is the scheme of which the singularities are to be resolved
  * ``f`` is a birational proper map 
          may for instance be BlowUpSequence or Lipman-style combination of blow-ups and normalization
  * ``Y`` is a regular scheme
"""
abstract type AbsDesingMor{
                           DomainType<:AbsCoveredScheme,
                           CodomainType<:AbsCoveredScheme,
                           BlowupMorphismType
                          } <: AbsCoveredSchemeMorphism{
                                                        DomainType,
                                                        CodomainType,
                                                        Nothing,
                                                        BlowupMorphismType
                                                       }
end

########################################################################
# An abstract type for blowup morphisms.
#
# This should also comprise sequences of simple blowups leading
# to a partial or full resolution of singularities. The interface
# is specified below.
########################################################################
abstract type AbsBlowupMorphism{DomainType<:AbsCoveredScheme,
                                  CodomainType<:AbsCoveredScheme,
                                  BlowupMorphismType
                                 } <: AbsDesingMor{DomainType, CodomainType, BlowupMorphismType}
end

########################################################################
# An abstract type for classical blowups of ideal sheaves.
#
# This can either be a BlowupMorphism as below, but also a special
# toric morphism induced by fan subdivisions.
########################################################################
abstract type AbsSimpleBlowupMorphism{DomainType<:AbsCoveredScheme,
                                     CodomainType<:AbsCoveredScheme,
                                     BlowupMorphismType
    } <: AbsBlowupMorphism{
                             DomainType,
                             CodomainType,
                             BlowupMorphismType
                            }
end


########################################################################
# BlowupMorphism
#
# A datastructure to maintain all information necessary to effectively
# handle blowups. This is work in progress and will one day serve as
# a building block for sequences of blowups
########################################################################

@doc raw"""
    BlowupMorphism

A datastructure to encode blowups of covered schemes in some sheaves of ideals.

It is described as a morphism from the new scheme to the blown-up scheme, with
information about its center (i.e. the ideal sheaves blown-up in the bottom
scheme) and its exceptional locus (i.e. the preimage of the center under the
blowup).

# Examples
```jldoctest
julia> R, (x,y,z) = QQ[:x, :y, :z];

julia> A3 = spec(R)
Spectrum
  of multivariate polynomial ring in 3 variables x, y, z
    over rational field

julia> I = ideal(R, [x,y,z])
Ideal generated by
  x
  y
  z

julia> bl = blow_up(A3, I)
Blowup
  of scheme over QQ covered with 1 patch
    1b: [x, y, z]   affine 3-space
  in sheaf of ideals with restriction
    1b: Ideal (x, y, z)
with domain
  scheme over QQ covered with 3 patches
    1a: [(s1//s0), (s2//s0), x]   scheme(0, 0, 0)
    2a: [(s0//s1), (s2//s1), y]   scheme(0, 0, 0)
    3a: [(s0//s2), (s1//s2), z]   scheme(0, 0, 0)
and exceptional divisor
  effective cartier divisor defined by
    sheaf of ideals with restrictions
      1a: Ideal (x)
      2a: Ideal (y)
      3a: Ideal (z)

julia> E = exceptional_divisor(bl)
Effective cartier divisor
  on scheme over QQ covered with 3 patches
    1: [(s1//s0), (s2//s0), x]   scheme(0, 0, 0)
    2: [(s0//s1), (s2//s1), y]   scheme(0, 0, 0)
    3: [(s0//s2), (s1//s2), z]   scheme(0, 0, 0)
defined by
  sheaf of ideals with restrictions
    1: Ideal (x)
    2: Ideal (y)
    3: Ideal (z)

julia> Z = center(bl)
Sheaf of ideals
  on scheme over QQ covered with 1 patch
    1: [x, y, z]   affine 3-space
with restriction
  1: Ideal (x, y, z)
```
"""
@attributes mutable struct BlowupMorphism{
     DomainType<:AbsCoveredScheme, # Not a concrete type in general because this is lazy
     CodomainType<:AbsCoveredScheme,
   } <: AbsSimpleBlowupMorphism{
                                  DomainType,
                                  CodomainType,
                                  BlowupMorphism{DomainType, CodomainType}
                                 }
  projective_bundle::CoveredProjectiveScheme
  codomain::CodomainType   # in general a CoveredScheme
  center::AbsIdealSheaf      # on codomain
  projection::AbsCoveredSchemeMorphism
  domain::AbsCoveredScheme # in general a CoveredScheme
  exceptional_divisor::EffectiveCartierDivisor

  function BlowupMorphism(
      IP::CoveredProjectiveScheme,
      I::AbsIdealSheaf
    )
    X = base_scheme(IP)
    X === scheme(I) || error("ideal sheaf not compatible with blown up variety")
    return new{AbsCoveredScheme, typeof(X)}(IP, X, I)
  end
end

########################################################################
# Resolutions of singularities                                         #
########################################################################

@doc raw"""
    BlowUpSequence{
    DomainType<:AbsCoveredScheme,
    CodomainType<:AbsCoveredScheme
   } <: AbsDesingMor{
                                 DomainType,
                                 CodomainType,
                                }


"""
@attributes mutable struct BlowUpSequence{
                                          DomainType<:AbsCoveredScheme,
                                          CodomainType<:AbsCoveredScheme
                                         }<:AbsBlowupMorphism{
                                                                DomainType, CodomainType, 
                                                                BlowUpSequence{DomainType, CodomainType}
                                                               }
  maps::Vector{<:BlowupMorphism}                 # count right to left:
                                                 # original scheme is codomain of map 1
  
  embeddings::Vector{<:AbsCoveredSchemeMorphism} # if set,
                                                 # assert codomain(maps[i])===codomain(embeddings[i]) 
  # boolean flags
  is_embedded::Bool                              # do not set embeddings, ex_mult, controlled_transform etc
                                                 #     if is_embedded == false
  resolves_sing::Bool                            # domain(maps[end]) smooth?
  is_trivial::Bool                               # codomain already smooth?
  is_strong::Bool                                # snc divisors ensured?
  transform_type::Symbol                         # can be :strict, :weak or :control
                                                 #     only relevant for is_embedded == true

  # fields for caching, may be filled during computation
  ex_div::Vector{<:EffectiveCartierDivisor}      # list of exc. divisors arising from individual steps
                                                 # lives in domain(maps[end])
  control::Int                                   # value of control for controlled transform
  ex_mult::Vector{Int}                           # multiplicities of exceptional divisors removed from
                                                 # total transform, not set for is_embedded == false
                                                 # or transform_type == strict
  controlled_transform::AbsIdealSheaf            # holds weak or controlled transform according to transform_type
  dont_meet::Vector{Tuple{Int,Int}}              # mostly for dim=2: intersections which cannot exist according
                                                 # to intermediate computations
  caution_multi_charts::Vector{Tuple{Int,Int}}   # only for dim=2: intersection of divisors not
                                                 # entirely visible in a single chart


  # fields for caching to be filled a posteriori (on demand, only if partial_res==false)
  underlying_morphism::CompositeCoveredSchemeMorphism{DomainType, CodomainType}
  exceptional_divisor::CartierDivisor            # exceptional divisor of composed_map
  exceptional_locus::WeilDivisor                 # exceptional locus of composed map
  exceptional_divisor_on_X::CartierDivisor          # exceptional divisor of composed_map
                                                 # restricted to domain(embeddings[end])

  function BlowUpSequence(maps::Vector{<:BlowupMorphism})
    n = length(maps)
    for i in 1:n-1
      @assert domain(maps[i]) === codomain(maps[i+1]) "not a sequence of morphisms"
    end
    return new{typeof(domain(maps[end])),typeof(codomain(first(maps)))}(maps)
  end
end

########################################################################
# strict transforms of ideal sheaves                                   #
########################################################################

@attributes mutable struct StrictTransformIdealSheaf{SpaceType, OpenType, OutputType,
                                                     RestrictionType
                                                    } <: AbsIdealSheaf{
                                                                       SpaceType, OpenType,
                                                                       OutputType, RestrictionType
                                                                      }
  morphism::AbsSimpleBlowupMorphism
  orig::AbsIdealSheaf
  underlying_presheaf::AbsPreSheaf

  function StrictTransformIdealSheaf(
      f::AbsSimpleBlowupMorphism,
      J::AbsIdealSheaf
    )
    @assert scheme(J) === codomain(f)
    X = domain(f)
    Ipre = PreSheafOnScheme(X,
                      OpenType=AbsAffineScheme, OutputType=Ideal,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes_without_affine_scheme_open_subscheme(X)
                     )
    I = new{typeof(X), AbsAffineScheme, Ideal, Map}(f, J, Ipre)
    return I
  end
end

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

    
