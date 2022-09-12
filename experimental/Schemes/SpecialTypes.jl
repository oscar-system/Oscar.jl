export PrincipalOpenSubset
export ambient_scheme, complement_equation, inclusion_morphism

export OpenInclusion
export complement_ideal, complement_scheme

export ClosedEmbedding
export image_ideal

@attributes mutable struct PrincipalOpenSubset{BRT, RT, AmbientType} <: AbsSpec{BRT, RT}
  X::AmbientType
  U::Spec{BRT, RT}
  f::RingElem
  inc::SpecMor

  function PrincipalOpenSubset(X::AbsSpec, f::RingElem)
    parent(f) == ambient_ring(X) || error("element does not belong to the correct ring")
    U = hypersurface_complement(X, [f])
    return new{base_ring_type(X), ring_type(U), typeof(X)}(X, U, f)
  end
  
  function PrincipalOpenSubset(X::AbsSpec, f::Vector{RingElemType}) where {RingElemType<:MPolyElem}
    all(x->(parent(x) == ambient_ring(X)), f) || error("element does not belong to the correct ring")
    U = hypersurface_complement(X, f)
    return new{base_ring_type(X), ring_type(U), typeof(X)}(X, U, prod(f))
  end
end

underlying_scheme(U::PrincipalOpenSubset) = U.U
ambient_scheme(U::PrincipalOpenSubset) = U.X
complement_equation(U::PrincipalOpenSubset) = U.f::elem_type(ambient_ring(U))

function inclusion_morphism(U::PrincipalOpenSubset) 
  if !isdefined(U, :inc)
    X = ambient_scheme(U)
    inc = SpecMor(U, X, hom(OO(X), OO(U), gens(OO(U))))
    U.inc = inc
  end
  return U.inc
end

@attributes mutable struct OpenInclusion{DomainType, CodomainType, PullbackType}<:AbsSpecMor{DomainType, CodomainType, PullbackType, OpenInclusion, Nothing}
  inc::SpecMor{DomainType, CodomainType, PullbackType}
  I::Ideal
  Z::Spec

  function OpenInclusion(f::AbsSpecMor, I::Ideal; check::Bool=true)
    U = domain(f)
    X = codomain(f)
    Z = subscheme(X, I)
    if check
      isempty(preimage(f, Z)) || error("image of the map is not contained in the complement of the vanishing locus of the ideal")
      #TODO: Do checks
    end
    return new{typeof(U), typeof(X), pullback_type(f)}(f, I, Z)
  end
end

underlying_morphism(f::OpenInclusion) = f.inc
complement_ideal(f::OpenInclusion) = f.I
complement_scheme(f::OpenInclusion) = f.Z


@attributes mutable struct ClosedEmbedding{DomainType, 
                                           CodomainType, 
                                           PullbackType
                                          }<:AbsSpecMor{DomainType, 
                                                        CodomainType, 
                                                        PullbackType, 
                                                        ClosedEmbedding,
                                                        Nothing
                                                       }
  inc::SpecMor{DomainType, CodomainType, PullbackType}
  I::Ideal
  U::SpecOpen

  function ClosedEmbedding(X::AbsSpec, I::Ideal)
    base_ring(I) == OO(X) || error("ideal does not belong to the correct ring")
    Y = subscheme(X, I)
    inc = SpecMor(Y, X, hom(OO(X), OO(Y), gens(OO(Y))))
    return new{typeof(Y), typeof(X), pullback_type(inc)}(inc, I)
  end
end

underlying_morphism(f::ClosedEmbedding) = f.inc

image_ideal(f::ClosedEmbedding) = f.I::ideal_type(OO(codomain(f)))

function complement(f::ClosedEmbedding)
  if !isdefined(f, :U)
    U = SpecOpen(codomain(f), image_ideal(f))
    f.U = U
  end
  return f.U
end

ideal_type(::Type{RT}) where {RT<:MPolyRing} = MPolyIdeal{elem_type(RT)}
ideal_type(::Type{RT}) where {PolyType, RT<:MPolyQuo{PolyType}} = MPolyQuoIdeal{PolyType}

export SimpleGlueing

@attributes mutable struct SimpleGlueing{LST<:AbsSpec, 
                                         RST<:AbsSpec, 
                                         LOT<:PrincipalOpenSubset, 
                                         ROT<:PrincipalOpenSubset, 
                                         LMT<:AbsSpecMor, 
                                         RMT<:AbsSpecMor
                                        } <: AbsGlueing{LST, RST, LOT, ROT, LMT, RMT} 
  X::LST
  Y::RST
  U::LOT
  V::ROT
  f::LMT
  g::RMT

  function SimpleGlueing(
      X::AbsSpec, Y::AbsSpec, 
      f::AbsSpecMor{<:PrincipalOpenSubset}, 
      g::AbsSpecMor{<:PrincipalOpenSubset};
      check::Bool=true
    )
    U = domain(f)
    V = domain(g)
    X == ambient_scheme(U) && Y == ambient_scheme(V) || error("schemes are not compatible")
    domain(f) == codomain(g) && domain(g) == codomain(f) || error("maps are not compatible")
    if check
      is_identity_map(compose(f, g)) || error("maps are not inverse to each other")
      is_identity_map(compose(g, f)) || error("maps are not inverse to each other")
    end
    set_attribute!(f, :inverse, g)
    set_attribute!(g, :inverse, f)
    return new{typeof(X), typeof(Y), 
               typeof(U), typeof(V), 
               typeof(f), typeof(g)
              }(X, Y, U, V, f, g)
  end
end

patches(G::SimpleGlueing) = (G.X, G.Y)
glueing_morphisms(G::SimpleGlueing) = (G.f, G.g)
glueing_domains(G::SimpleGlueing) = (G.U, G.V)

@attr SimpleGlueing function inverse(G::SimpleGlueing)
  Ginv = SimpleGlueing(G.Y, G.X, G.g, G.f, check=false)
  set_attribute!(Ginv, :inverse, G)
  return Ginv
end

function compose(G::GT, H::GT) where {GT<:SimpleGlueing}
  if patches(G)[2] == patches(H)[2] 
    return compose(G, inverse(H))
  elseif patches(G)[1] == patches(H)[1]
    return compose(inverse(G), H)
  elseif patches(G)[1] == patches(H)[2]
    return compose(inverse(G), inverse(H))
  end
  X, Y = patches(G)
  Y == patches(H)[1] || error("Glueings not compatible")
  Z = patches(H)[2]
  f, f_inv = glueing_morphisms(G)
  g, g_inv = glueing_morphisms(H)
  U_new = PrincipalOpenSubset(domain(f), pullback(f)(complement_equation(domain(g))))
  W_new = PrincipalOpenSubset(domain(g_inv), pullback(g_inv)(complement_equation(domain(f_inv))))
  V_new = PrincipalOpenSubset(ambient_scheme(domain(g)), 
                              [complement_equation(domain(g)), complement_equation(domain(f_inv))]
                             )
  h = compose(restrict(f, U_new, V_new, check=false), 
              restrict(g, V_new, W_new, check=false))
  h_inv = compose(restrict(g_inv, W_new, V_new, check=false),
                  restrict(f_inv, V_new, U_new, check=false))
  set_attribute!(h, :inverse, h_inv)
  set_attribute!(h_inv, :inverse, h)
  return SimpleGlueing(X, Z, h, h_inv)
end

function restriction(G::SimpleGlueing, X::AbsSpec, Y::AbsSpec; check::Bool=true)
  U, V = glueing_domains(G)
  f, g = glueing_morphisms(G)
  if check
    is_closed_embedding(intersect(X, ambient(U)), ambient(U)) || error("the scheme is not a closed in the ambient scheme of the open set")
    is_closed_embedding(intersect(Y, ambient(V)), ambient(V)) || error("the scheme is not a closed in the ambient scheme of the open set")
  end
  UX = PrincipalOpenSubset(X, complement_equation(U))
  VY = PrincipalOpenSubset(Y, complement_equation(V))
  f_res = restrict(f, UX, VY, check=check)
  g_res = restrict(g, VY, UX, check=check)
  return SimpleGlueing(X, Y, f_res, g_res)
end
