
########################################################################
# Abstract gluings for affine schemes                                 #
########################################################################
@doc raw"""
    AbsGluing

A gluing of two affine schemes ``X`` and ``Y`` (the `patches`) along 
open subsets ``U`` in ``X`` and ``V`` in ``Y`` (the `gluing_domains`)
along mutual isomorphisms ``f : U ↔ V : g`` (the `gluing_morphisms`).
"""
abstract type AbsGluing{LeftAffineSchemeType<:AbsAffineScheme,
                         RightAffineSchemeType<:AbsAffineScheme,
                         LeftOpenType<:Scheme,
                         RightOpenType<:Scheme,
                         LeftMorType<:Map,
                         RightMorType<:Map
                        } end

########################################################################
# Concrete type for general gluings                                   #
########################################################################
@doc raw"""
    Gluing

Concrete instance of an `AbsGluing` for gluings of affine schemes 
``X ↩ U ≅ V ↪ Y`` along open subsets ``U`` and ``V`` of type `AffineSchemeOpenSubscheme`.
"""
@attributes mutable struct Gluing{
                                   LeftAffineSchemeType<:AbsAffineScheme,
                                   RightAffineSchemeType<:AbsAffineScheme,
                                   LeftOpenType<:AffineSchemeOpenSubscheme,
                                   RightOpenType<:AffineSchemeOpenSubscheme,
                                   LeftMorType<:AffineSchemeOpenSubschemeMor,
                                   RightMorType<:AffineSchemeOpenSubschemeMor
                                  } <: AbsGluing{
                                   LeftAffineSchemeType,
                                   RightAffineSchemeType,
                                   LeftOpenType,
                                   RightOpenType,
                                   LeftMorType,
                                   RightMorType
                                  }
  X::LeftAffineSchemeType
  Y::RightAffineSchemeType
  U::LeftOpenType
  V::RightOpenType
  f::LeftMorType # f : U → V
  g::RightMorType

  function Gluing(
      X::AbsAffineScheme, Y::AbsAffineScheme, f::AffineSchemeOpenSubschemeMor, g::AffineSchemeOpenSubschemeMor; check::Bool=true
    )
    ambient_scheme(domain(f)) === X || error("the domain of the gluing morphism is not an open subset of the first argument")
    ambient_scheme(codomain(f)) === Y || error("the codomain of the gluing morphism is not an open subset of the second argument")
    (domain(f) === codomain(g) &&
     domain(g) ===  codomain(f)) || error("maps can not be isomorphisms")
    (OO(domain(f)) === OO(codomain(g)) &&
    OO(domain(g)) ===  OO(codomain(f))) || error("domain/codomain mismatch: maps can not be mutually inverse")
    @check compose(f, g) == identity_map(domain(f)) "gluing maps are not inverse of each other"
    @check compose(g, f) == identity_map(domain(g)) "gluing maps are not inverse of each other"
    return new{typeof(X), typeof(Y),
               typeof(domain(f)), typeof(domain(g)),
               typeof(f), typeof(g)
              }(X, Y, domain(f), domain(g), f, g)
  end
end

########################################################################
# Special type for simple gluings of affine schemes along principal
# open subsets
#
# SimpleGluing is for gluings X ↩ U ≅ V ↪ Y along principal
# open subsets U ⊂ X and V ⊂ Y along identifications f : U ↔ V : g.
# For general gluings it can not be guaranteed to have this setup,
# but it is a situation often encountered and with significant
# simplification of underlying algorithms in the background.
# Hence, the special type.
########################################################################
@doc raw"""
    SimpleGluing

Concrete instance of an `AbsGluing` for gluings of affine schemes 
``X ↩ U ≅ V ↪ Y`` along open subsets ``U`` and ``V`` of type 
`PrincipalOpenSubset`.
"""
@attributes mutable struct SimpleGluing{LST<:AbsAffineScheme,
                                         RST<:AbsAffineScheme,
                                         LOT<:PrincipalOpenSubset,
                                         ROT<:PrincipalOpenSubset,
                                         LMT<:AbsAffineSchemeMor,
                                         RMT<:AbsAffineSchemeMor
                                        } <: AbsGluing{LST, RST, LOT, ROT, LMT, RMT}
  X::LST
  Y::RST
  U::LOT
  V::ROT
  f::LMT
  g::RMT

  function SimpleGluing(
      X::AbsAffineScheme, Y::AbsAffineScheme,
      f::AbsAffineSchemeMor{<:PrincipalOpenSubset},
      g::AbsAffineSchemeMor{<:PrincipalOpenSubset};
      check::Bool=true
    )
    U = domain(f)
    V = domain(g)
    X === ambient_scheme(U) && Y === ambient_scheme(V) || error("schemes are not compatible")
    domain(f) === codomain(g) && domain(g) === codomain(f) || error("maps are not compatible")
    @check is_identity_map(compose(f, g)) "maps are not inverse to each other"
    @check is_identity_map(compose(g, f)) "maps are not inverse to each other"
    set_attribute!(f, :inverse, g)
    set_attribute!(g, :inverse, f)
    return new{typeof(X), typeof(Y),
               typeof(U), typeof(V),
               typeof(f), typeof(g)
              }(X, Y, U, V, f, g)
  end
end

