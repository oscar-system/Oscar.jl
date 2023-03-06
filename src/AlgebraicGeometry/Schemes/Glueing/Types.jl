export AbsGlueing, Glueing, SimpleGlueing

########################################################################
# Abstract glueings for affine schemes                                 #
########################################################################
@Markdown.doc """
    AbsGlueing

A glueing of two affine schemes ``X`` and ``Y`` (the `patches`) along 
open subsets ``U`` in ``X`` and ``V`` in ``Y`` (the `glueing_domains`)
along mutual isomorphisms ``f : U ↔ V : g`` (the `glueing_morphisms`).
"""
abstract type AbsGlueing{LeftSpecType<:AbsSpec,
                         RightSpecType<:AbsSpec,
                         LeftOpenType<:Scheme,
                         RightOpenType<:Scheme,
                         LeftMorType<:Hecke.Map,
                         RightMorType<:Hecke.Map
                        } end

########################################################################
# Concrete type for general glueings                                   #
########################################################################
@Markdown.doc """
    Glueing

Concrete instance of an `AbsGlueing` for glueings of affine schemes 
``X ↩ U ≅ V ↪ Y`` along open subsets ``U`` and ``V`` of type `SpecOpen`.
"""
@attributes mutable struct Glueing{
                                   LeftSpecType<:AbsSpec,
                                   RightSpecType<:AbsSpec,
                                   LeftOpenType<:SpecOpen,
                                   RightOpenType<:SpecOpen,
                                   LeftMorType<:SpecOpenMor,
                                   RightMorType<:SpecOpenMor
                                  } <: AbsGlueing{
                                   LeftSpecType,
                                   RightSpecType,
                                   LeftOpenType,
                                   RightOpenType,
                                   LeftMorType,
                                   RightMorType
                                  }
  X::LeftSpecType
  Y::RightSpecType
  U::LeftOpenType
  V::RightOpenType
  f::LeftMorType # f : U → V
  g::RightMorType

  function Glueing(
      X::AbsSpec, Y::AbsSpec, f::SpecOpenMor, g::SpecOpenMor; check::Bool=true
    )
    ambient_scheme(domain(f)) === X || error("the domain of the glueing morphism is not an open subset of the first argument")
    ambient_scheme(codomain(f)) === Y || error("the codomain of the glueing morphism is not an open subset of the second argument")
    (domain(f) === codomain(g) &&
     domain(g) ===  codomain(f)) || error("maps can not be isomorphisms")
    (OO(domain(f)) === OO(codomain(g)) &&
    OO(domain(g)) ===  OO(codomain(f))) || error("domain/codomain mismatch: maps can not be mutually inverse")
    if check
      compose(f, g) == identity_map(domain(f)) || error("glueing maps are not inverse of each other")
      compose(g, f) == identity_map(domain(g)) || error("glueing maps are not inverse of each other")
    end
    return new{typeof(X), typeof(Y),
               typeof(domain(f)), typeof(domain(g)),
               typeof(f), typeof(g)
              }(X, Y, domain(f), domain(g), f, g)
  end
end

########################################################################
# Special type for simple glueings of affine schemes along principal
# open subsets
#
# SimpleGlueing is for glueings X ↩ U ≅ V ↪ Y along principal
# open subsets U ⊂ X and V ⊂ Y along identifications f : U ↔ V : g.
# For general glueings it can not be guaranteed to have this setup,
# but it is a situation often encountered and with significant
# simplification of underlying algorithms in the background.
# Hence, the special type.
########################################################################
@Markdown.doc """
    SimpleGlueing

Concrete instance of an `AbsGlueing` for glueings of affine schemes 
``X ↩ U ≅ V ↪ Y`` along open subsets ``U`` and ``V`` of type 
`PrincipalOpenSubset`.
"""
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
    X === ambient_scheme(U) && Y === ambient_scheme(V) || error("schemes are not compatible")
    domain(f) === codomain(g) && domain(g) === codomain(f) || error("maps are not compatible")
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

