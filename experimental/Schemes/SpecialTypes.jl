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
    parent(f) == OO(X) || error("element does not belong to the correct ring")
    U = hypersurface_complement(X, [f])
    return new{base_ring_type(X), ring_type(U), typeof(X)}(X, U, f)
  end
  
  function PrincipalOpenSubset(X::AbsSpec, f::Vector{RingElemType}) where {RingElemType<:MPolyElem}
    all(x->(parent(x) == OO(X)), f) || error("element does not belong to the correct ring")
    U = hypersurface_complement(X, f)
    return new{base_ring_type(X), ring_type(U), typeof(X)}(X, U, prod(f))
  end
end

underlying_scheme(U::PrincipalOpenSubset) = U.U
ambient_scheme(U::PrincipalOpenSubset) = U.X
complement_equation(U::PrincipalOpenSubset) = U.f::elem_type(OO(ambient_scheme(U)))

### assure compatibility with SpecOpen 
gens(U::PrincipalOpenSubset) = [lifted_numerator(complement_equation(U))]
getindex(U::PrincipalOpenSubset, i::Int) = (i == 1 ? U : error("index out of range"))

function SpecOpen(U::PrincipalOpenSubset) 
  X = ambient_scheme(U)
  h = lifted_numerator(complement_equation(U))
  return SpecOpen(X, [h])
end

function inclusion_morphism(U::PrincipalOpenSubset) 
  if !isdefined(U, :inc)
    X = ambient_scheme(U)
    inc = SpecMor(U, X, hom(OO(X), OO(U), gens(OO(U))))
    U.inc = inc
  end
  return U.inc
end

PrincipalOpenSubset(X::AbsSpec) = PrincipalOpenSubset(X, one(OO(X)))

function preimage(f::AbsSpecMor, U::PrincipalOpenSubset; check::Bool=true) 
  if ambient_scheme(U) != codomain(f) 
    Z = preimage(f, ambient_scheme(U), check=check)
    return PrincipalOpenSubset(Z, OO(Z)(pullback(f)(complement_equation(U)), check=false))
  end
  return PrincipalOpenSubset(domain(f), pullback(f)(complement_equation(U)))
end

@Markdown.doc """
    generic_fraction(a::MPolyLocalizedRingElem, U::PrincipalOpenSubset)

Given a regular function ``a ∈ 𝒪(U)`` on a principal open 
subset ``U ⊂ X`` of an affine scheme ``X``, return a 
fraction ``p/q`` in `Quot(P)` (where ``P`` is the `ambient_ring` of 
the `ambient` scheme ``X`` of ``U``) which represents ``a``
in the sense that the maximal extension of its restriction 
to ``U`` returns ``a``.
"""
function generic_fraction(a::MPolyLocalizedRingElem, U::PrincipalOpenSubset)
  X = ambient_scheme(U)
  parent(a) == OO(U) || error("domains are not compatible")
  return lifted_numerator(a)//lifted_denominator(a)
end

function generic_fraction(a::MPolyQuoLocalizedRingElem, U::PrincipalOpenSubset)
  X = ambient_scheme(U)
  parent(a) == OO(U) || error("domains are not compatible")
  return lifted_numerator(a)//lifted_denominator(a)
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


@Markdown.doc """
    ClosedEmbedding{DomainType, CodomainType, PullbackType}

A closed embedding ``f : X → Y`` of affine schemes ``X = Spec(S)`` 
into ``Y = Spec(R)`` such that ``S ≅ R/I`` via ``f`` for some 
ideal ``I ⊂ R``.
"""
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

@Markdown.doc """
    image_ideal(f::ClosedEmbedding)

For a closed embedding ``f : X → Y`` of affine schemes ``X = Spec(S)`` 
into ``Y = Spec(R)`` such that ``S ≅ R/I`` via ``f`` for some ideal 
``I ⊂ R`` this returns ``I``.
"""
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

########################################################################
# SimpleGlueing is for glueings X ↩ U ≅ V ↪ Y along principal 
# open subsets U ⊂ X and V ⊂ Y along identifications f : U ↔ V : g. 
# For general glueings it can not be guaranteed to have this setup, 
# but it is a situation often encountered and with significant 
# simplification of underlying algorithms in the background. 
# Hence, the special type.
########################################################################
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

function Glueing(
    X::AbsSpec, Y::AbsSpec, 
    f::AbsSpecMor{<:PrincipalOpenSubset}, 
    g::AbsSpecMor{<:PrincipalOpenSubset};
    check::Bool=true)
  return SimpleGlueing(X, Y, f, g, check=check)
end

@attr function is_dense(U::PrincipalOpenSubset)
  f = complement_equation(U)
  return is_non_zero_divisor(f, ambient_scheme(U))
end

### conversion
function Glueing(G::SimpleGlueing)
  X, Y = patches(G)
  f, g = glueing_morphisms(G)
  U, V = glueing_domains(G)
  UU = SpecOpen(U)
  VV = SpecOpen(V)
  ff = SpecOpenMor(UU, VV, [compose(restrict(f, UU[1], V, check=true),
                                    inclusion_map(V, Y))
                           ])
  gg = SpecOpenMor(VV, UU, [compose(restrict(g, VV[1], U, check=true),
                                    inclusion_map(U, X))
                           ])
  return Glueing(X, Y, ff, gg)
end

function compose(G::Glueing, H::SimpleGlueing) 
  return compose(G, Glueing(H))
end

function compose(G::SimpleGlueing, H::Glueing) 
  return compose(Glueing(G), H)
end

