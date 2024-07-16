
########################################################################
# Morphisms of coverings                                               #
########################################################################
@doc raw"""
    CoveringMorphism

A morphism ``f : C → D`` of two coverings. For every patch ``U`` of ``C`` this
provides a map `f[U']` of type `AffineSchemeMorType` from ``U' ⊂ U`` to
some patch `codomain(f[U])` in `D` for some affine patches ``U'`` covering ``U``.

**Note:** For two affine patches ``U₁, U₂ ⊂ U`` the codomains of `f[U₁]` and `f[U₂]`
do not need to coincide! However, given the gluings in `C` and `D`, all affine maps
have to coincide on their overlaps.
"""
mutable struct CoveringMorphism{DomainType<:Covering,
                                CodomainType<:Covering,
                                MorphismType<:AbsAffineSchemeMor,
                                BaseMorType
                               }
  domain::DomainType
  codomain::CodomainType
  morphisms::IdDict{<:AbsAffineScheme, <:MorphismType} # on a patch X of the domain covering, this
                                         # returns the morphism φ : X → Y to the corresponding
                                         # patch Y of the codomain covering.

  function CoveringMorphism(
      dom::DomainType,
      cod::CodomainType,
      mor::IdDict{<:AbsAffineScheme, MorphismType};
      check::Bool=true
    ) where {
             DomainType<:Covering,
             CodomainType<:Covering,
             MorphismType<:AbsAffineSchemeMor
            }
    # TODO: check domain/codomain compatibility
    # TODO: if check is true, check that all morphisms glue and that the domain patches
    # cover the basic patches of `dom`.
    @check
    for U in keys(mor)
      U in dom || error("patch $U of the map not found in domain")
      codomain(mor[U]) in cod || error("codomain patch not found")
    end
    # check that the whole domain is covered
    for U in basic_patches(dom)
      if !haskey(mor, U)
        !haskey(affine_refinements(dom), U) || error("patch $U of the domain not covered")
        found = false
        for (V, a) in affine_refinements(dom)[U]
          all(x->(haskey(mor, x)), affine_patches(V)) && (found = true)
        end
        !found && error("patch $U of the domain not covered")
      end
    end
    @check begin
      for U in patches(dom)
        for V in patches(dom)
          U === V && continue
          phi_U = mor[U]
          phi_V = mor[V]
          !haskey(gluings(dom), (U, V)) && continue
          glue = dom[U, V]
          f, g = gluing_morphisms(glue)
          !haskey(gluings(cod), (codomain(phi_U), codomain(phi_V))) && error("gluing not found in the codomain")
          U_cod = codomain(phi_U)
          V_cod = codomain(phi_V)
          glue_cod = cod[U_cod, V_cod]
          f_cod, g_cod = gluing_morphisms(glue_cod)
          phi_U_res = restrict(phi_U, domain(f), domain(f_cod))
          phi_V_res = restrict(phi_V, domain(g), domain(g_cod))
          compose(f, phi_V_res) == compose(phi_U_res, f_cod) || error("restrictions do not commute")
          compose(g, phi_U_res) == compose(phi_V_res, g_cod) || error("restrictions do not commute")
        end
      end
    end
    return new{DomainType, CodomainType, MorphismType, Nothing}(dom, cod, mor)
  end
end

