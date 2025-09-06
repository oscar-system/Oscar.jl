@doc raw"""
    morphism(P::AbsProjectiveScheme, Q::AbsProjectiveScheme, f::Map; check::Bool=true )

Given a morphism ``f : T → S`` of the `homogeneous_coordinate_ring`s of `Q` and `P`, respectively, 
construct the associated morphism of projective schemes.
"""
function morphism(P::AbsProjectiveScheme, Q::AbsProjectiveScheme, f::Map; check::Bool=true )
  return ProjectiveSchemeMor(P, Q, f, check=check)
end

@doc raw"""
    morphism(P::AbsProjectiveScheme, Q::AbsProjectiveScheme, f::Map, h::SchemeMor; check::Bool=true )

Suppose ``P ⊂ ℙʳ_A`` and ``Q ⊂ ℙˢ_B`` are projective schemes, ``h : Spec(A) → Spec(B)`` is a 
morphism of their `base_scheme`s, and ``f : T → S`` a morphism of the 
`homogeneous_coordinate_ring`s of `Q` and `P` over ``h^* : B → A``.
This constructs the associated morphism of projective schemes.
"""
function morphism(P::AbsProjectiveScheme, Q::AbsProjectiveScheme, f::Map, h::SchemeMor; check::Bool=true )
  return ProjectiveSchemeMor(P, Q, f, h, check=check)
end

@doc raw"""
    morphism(X::AbsProjectiveScheme, Y::AbsProjectiveScheme, a::Vector{<:RingElem})

Suppose ``X ⊂ ℙʳ`` and ``Y ⊂ ℙˢ`` are projective schemes over the same `base_scheme`.
Construct the morphism of projective schemes associated to the morphism of graded rings 
which takes the generators of the `homogeneous_coordinate_ring` of ``Y`` to the elements 
in `a` of the `homogeneous_coordinate_ring` of ``X``.
"""
function morphism(X::AbsProjectiveScheme, Y::AbsProjectiveScheme, a::Vector{<:RingElem})
  return ProjectiveSchemeMor(X, Y, a)
end

function ProjectiveSchemeMor(
    X::AbsProjectiveScheme, 
    Y::AbsProjectiveScheme, 
    a::Vector{<:RingElem}
  )
  base_ring(X) === base_ring(Y) || error("projective schemes must be defined over the same base ring")
  Q = homogeneous_coordinate_ring(X)
  all(x->parent(x)===Q, a) || return ProjectiveSchemeMor(X, Y, Q.(a))
  P = homogeneous_coordinate_ring(Y)
  return ProjectiveSchemeMor(X, Y, hom(P, Q, a, check=false))
end

function compose(f::AbsProjectiveSchemeMorphism, g::AbsProjectiveSchemeMorphism)
  return ProjectiveSchemeMor(domain(f), codomain(g), compose(pullback(g), pullback(f)))
end

### additional constructors

function fiber_product(f::Map{DomType, CodType}, P::AbsProjectiveScheme{DomType}) where {DomType<:Ring, CodType<:Ring}
  R = base_ring(P) 
  R === domain(f) || error("rings not compatible")
  Rnew = codomain(f)
  S = ambient_coordinate_ring(P)
  Qambient = projective_space(Rnew, symbols(S))
  Snew = ambient_coordinate_ring(Qambient)
  phi = hom(S, Snew, f, gens(Snew), check=false)
  Q = subscheme(Qambient, ideal(homogeneous_coordinate_ring(Qambient), phi.(gens(defining_ideal(P)))))
  return Q, ProjectiveSchemeMor(Q, P, hom(homogeneous_coordinate_ring(P), 
                                          homogeneous_coordinate_ring(Q), 
                                          f, gens(homogeneous_coordinate_ring(Q)),
                                          check=false
                                         ),
                                check=false
                               )
end

function fiber_product(
    f::AbsAffineSchemeMor, 
    P::AbsProjectiveScheme{<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}}
  )
  codomain(f) == base_scheme(P) || error("codomain and base_scheme are incompatible")
  X = domain(f)
  Y = codomain(f)
  Q_ambient = projective_space(X, symbols(ambient_coordinate_ring(P)))
  help_map = hom(ambient_coordinate_ring(P),
                 ambient_coordinate_ring(Q_ambient),
                 pullback(f),
                 gens(ambient_coordinate_ring(Q_ambient)),
                 check=false
                )
  SQ = homogeneous_coordinate_ring(Q_ambient)
  I = ideal(SQ, SQ.(help_map.(gens(defining_ideal(P)))))
  Q = subscheme(Q_ambient, I)
  return Q, ProjectiveSchemeMor(Q, P, 
                                hom(homogeneous_coordinate_ring(P),
                                    homogeneous_coordinate_ring(Q),
                                    pullback(f),
                                    gens(homogeneous_coordinate_ring(Q)), 
                                    check=false
                                   ),
                                check=false
                               )
end

fiber_product(X::AbsAffineScheme, 
              P::AbsProjectiveScheme{<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}}
             ) = fiber_product(inclusion_morphism(X, base_scheme(P)), P)

### canonical map constructors

@doc raw"""
    inclusion_morphism(P::T, Q::T)

Assuming that ``P ⊂ Q`` is a subscheme, both proper over an inclusion of 
their base schemes, this returns the associated `ProjectiveSchemeMor`.
"""
function inclusion_morphism(
    P::AbsProjectiveScheme{<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}},
    Q::AbsProjectiveScheme{<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}}
  )
  X = base_scheme(P)
  Y = base_scheme(Q)
  f = inclusion_morphism(X, Y) # will throw if X and Y are not compatible
  return ProjectiveSchemeMor(P, Q, 
                             hom(homogeneous_coordinate_ring(Q),
                                 homogeneous_coordinate_ring(P),
                                 pullback(f), 
                                 gens(homogeneous_coordinate_ring(P)),
                                 check=false
                                ),
                             check=false
                            )
end

function inclusion_morphism(P::T, Q::T) where {T<:AbsProjectiveScheme{<:AbstractAlgebra.Ring}}
  A = base_ring(Q)
  B = base_ring(P)
  A === B || error("can not compare schemes for non-equal base rings") # TODO: Extend by check for canonical maps, once they are available
  return ProjectiveSchemeMor(P, Q, 
                             hom(homogeneous_coordinate_ring(Q),
                                 homogeneous_coordinate_ring(P),
                                 gens(homogeneous_coordinate_ring(P)),
                                 check=false
                                ),
                             check=false
                            )
end

id_hom(P::AbsProjectiveScheme) = ProjectiveSchemeMor(P, P,
                                                        hom(homogeneous_coordinate_ring(P),
                                                            homogeneous_coordinate_ring(P),
                                                            gens(homogeneous_coordinate_ring(P)),
                                                            check=false
                                                           ),
                                                        check=false
                                                       )

function sub(P::AbsProjectiveScheme, I::Ideal)
  @req base_ring(I) === homogeneous_coordinate_ring(P) "ideal must be defined in the homogeneous coordinate ring of the scheme"
  inc = ProjectiveClosedEmbedding(P, I)
  set_attribute!(domain(inc), :ambient_space, ambient_space(P))
  return domain(inc), inc
end
  
function sub(P::AbsProjectiveScheme, f::RingElem)
  I = ideal(homogeneous_coordinate_ring(P), f)
  return sub(P, I)
end

function sub(P::AbsProjectiveScheme, f::Vector{<:RingElem})
  I = ideal(homogeneous_coordinate_ring(P), f)
  return sub(P, I)
end

function compose(f::ProjectiveClosedEmbedding, g::ProjectiveClosedEmbedding)
  X = domain(f)
  Y = codomain(f)
  Y === domain(g) || error("domain and codomain not compatible")
  Z = codomain(g)
  pb = compose(pullback(g), pullback(f))
  Ig = image_ideal(g)
  SY = homogeneous_coordinate_ring(Y)
  SZ = homogeneous_coordinate_ring(Z)
  If = image_ideal(f)
  push_If = ideal(SZ, [preimage(pullback(g), x) for x in gens(If)])
  J = push_If + Ig
  return ProjectiveClosedEmbedding(compose(underlying_morphism(f), underlying_morphism(g)), J, check=false)
end

function ambient_embedding(X::AbsProjectiveScheme)
  IP = ambient_space(X)
  S = homogeneous_coordinate_ring(IP)
  T = homogeneous_coordinate_ring(X)
  I = defining_ideal(X)
  pb = hom(S, T, gens(T); check=false)
  inc_sub = ProjectiveSchemeMor(X, IP, pb, check=false)
  return ProjectiveClosedEmbedding(inc_sub, I, check=false)
end

function base_change(phi::Any, IP::AbsProjectiveScheme)
  kk = base_ring(IP)
  KK = parent(phi(zero(kk)))
  S = homogeneous_coordinate_ring(IP)
  SS, pb_Phi = change_base_ring(phi, S)
  result = projective_scheme(SS)
  psi = morphism(result, IP, pb_Phi; check=false)
  return result, psi
end

########################################################################
# Rational maps
########################################################################
@doc raw"""
    rational_map(P::AbsProjectiveScheme, Q::AbsProjectiveScheme, f::Map; check::Bool=true )

Given a homomorphism ``f : T → S`` of the `homogeneous_coordinate_ring`s of `Q` and `P`, respectively, 
construct the associated rational map of projective varieties.
"""
function rational_map(P::AbsProjectiveScheme, Q::AbsProjectiveScheme, f::Map; check::Bool=true )
  return RationalMap(P, Q, f, check=check)
end

@doc raw"""
    rational_map(X::AbsProjectiveScheme, Y::AbsProjectiveScheme, a::Vector{<:RingElem})

Suppose ``X ⊂ ℙʳ`` and ``Y ⊂ ℙˢ`` are projective varieties.
Construct the rational map associated to the morphism of graded rings 
which takes the generators of the `homogeneous_coordinate_ring` of ``Y`` to the elements 
in `a` of the `homogeneous_coordinate_ring` of ``X``.
"""
function rational_map(X::AbsProjectiveScheme, Y::AbsProjectiveScheme, a::Vector{<:RingElem}; check::Bool=true)
  SX = homogeneous_coordinate_ring(X)
  SY = homogeneous_coordinate_ring(Y)
  f = hom(SY, SX, a; check)
  return RationalMap(X, Y, f)
end

