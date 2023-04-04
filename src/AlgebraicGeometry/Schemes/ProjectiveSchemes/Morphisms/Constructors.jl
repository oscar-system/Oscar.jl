
function ProjectiveSchemeMor(
    X::AbsProjectiveScheme, 
    Y::AbsProjectiveScheme, 
    a::Vector{<:RingElem}
  )
  base_ring(X) === base_ring(Y) || error("projective schemes must be defined over the same base ring")
  Q = homogeneous_coordinate_ring(X)
  all(x->parent(x)===Q, a) || return ProjectiveSchemeMor(X, Y, Q.(a))
  P = homogeneous_coordinate_ring(Y)
  return ProjectiveSchemeMor(X, Y, hom(P, Q, a))
end

function compose(f::ProjectiveSchemeMor, g::ProjectiveSchemeMor)
  return ProjectiveSchemeMor(domain(f), codomain(g), compose(pullback(g), pullback(f)))
end

### additional constructors

function fiber_product(f::Hecke.Map{DomType, CodType}, P::ProjectiveScheme{DomType}) where {DomType<:Ring, CodType<:Ring}
  R = base_ring(P) 
  R == domain(f) || error("rings not compatible")
  Rnew = codomain(f)
  S = ambient_coordinate_ring(P)
  Qambient = projective_space(Rnew, symbols(S))
  Snew = ambient_coordinate_ring(Qambient)
  phi = hom(S, Snew, f, gens(Snew))
  Q = subscheme(Qambient, ideal(homogeneous_coordinate_ring(Qambient), phi.(gens(defining_ideal(P)))))
  return Q, ProjectiveSchemeMor(Q, P, hom(homogeneous_coordinate_ring(P), 
                                          homogeneous_coordinate_ring(Q), 
                                          f, gens(homogeneous_coordinate_ring(Q))))
end

function fiber_product(
    f::AbsSpecMor, 
    P::ProjectiveScheme{<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}}
  )
  codomain(f) == base_scheme(P) || error("codomain and base_scheme are incompatible")
  X = domain(f)
  Y = codomain(f)
  Q_ambient = projective_space(X, symbols(ambient_coordinate_ring(P)))
  help_map = hom(ambient_coordinate_ring(P),
                 ambient_coordinate_ring(Q_ambient),
                 pullback(f),
                 gens(ambient_coordinate_ring(Q_ambient))
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
                                   )
                               )
end

fiber_product(X::AbsSpec, 
              P::ProjectiveScheme{<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}}
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
                                 gens(homogeneous_coordinate_ring(P))
                                )
                            )
end

function inclusion_morphism(P::T, Q::T) where {T<:AbsProjectiveScheme{<:AbstractAlgebra.Ring}}
  A = base_ring(Q)
  B = base_ring(P)
  A === B || error("can not compare schemes for non-equal base rings") # TODO: Extend by check for canonical maps, once they are available
  return ProjectiveSchemeMor(P, Q, 
                             hom(homogeneous_coordinate_ring(Q),
                                 homogeneous_coordinate_ring(P),
                                 gens(homogeneous_coordinate_ring(P))
                                )
                            )
end

identity_map(P::ProjectiveScheme) = ProjectiveSchemeMor(P, P, 
                                                        hom(homogeneous_coordinate_ring(P),
                                                            homogeneous_coordinate_ring(P),
                                                            gens(homogeneous_coordinate_ring(P))
                                                           )
                                                       )

