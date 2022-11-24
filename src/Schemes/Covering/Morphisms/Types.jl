export CoveringMorphism

########################################################################
# Morphisms of coverings                                               #
########################################################################
@Markdown.doc """
    CoveringMorphism

A morphism ``f : C → D`` of two coverings. For every patch ``U`` of ``C`` this
provides a map `f[U']` of type `SpecMorType` from ``U' ⊂ U`` to
some patch `codomain(f[U])` in `D` for some affine patches ``U'`` covering ``U``.

**Note:** For two affine patches ``U₁, U₂ ⊂ U`` the codomains of `f[U₁]` and `f[U₂]`
do not need to coincide! However, given the glueings in `C` and `D`, all affine maps
have to coincide on their overlaps.
"""
mutable struct CoveringMorphism{DomainType<:Covering, CodomainType<:Covering, BaseMorType}
  domain::DomainType
  codomain::CodomainType
  morphisms::IdDict{<:AbsSpec, <:AbsSpecMor} # on a patch X of the domain covering, this
                                         # returns the morphism φ : X → Y to the corresponding
                                         # patch Y of the codomain covering.

  function CoveringMorphism(
      dom::DomainType,
      cod::CodomainType,
      mor::IdDict{<:AbsSpec, <:AbsSpecMor};
      check::Bool=true
    ) where {
             DomainType<:Covering,
             CodomainType<:Covering
            }
    # TODO: check domain/codomain compatibility
    # TODO: if check is true, check that all morphisms glue and that the domain patches
    # cover the basic patches of `dom`.
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
    return new{DomainType, CodomainType, Nothing}(dom, cod, mor)
  end
end

