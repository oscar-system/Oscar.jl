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
