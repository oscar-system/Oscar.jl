export SpecOpenMor

########################################################################
# Morphisms of Zariski-open subsets of affine schemes                  #
########################################################################
@Markdown.doc """
    SpecOpenMor{DomainType<:SpecOpen, CodomainType<:SpecOpen}

Morphisms ``f : U → V`` of open sets ``U ⊂ X`` and ``V ⊂ Y`` of affine schemes.
These are stored as morphisms ``fᵢ: Uᵢ→ Y`` on the affine patches
``Uᵢ`` of ``U``.

The type parameters stand for the following: When ``X = Spec(R)`` and
``Y = Spec(S)`` with ``R = (𝕜[x₁,…,xₘ]/I)[A⁻¹]`` and ``S = (𝕜[y₁,…,yₙ]/J)[B⁻¹]``
then

 * `DomainType` is the type of the domain;
 * `CodomainType` is the type of the codomain;
affine patches of the domain to the affine ambient scheme of the codomain.
"""
mutable struct SpecOpenMor{DomainType<:SpecOpen,
                           CodomainType<:SpecOpen
                          }<:SchemeMor{DomainType, CodomainType, SpecOpenMor, Nothing}
  domain::DomainType
  codomain::CodomainType
  maps_on_patches::Vector{AbsSpecMor}

  # fields used for caching
  inverse::SpecOpenMor
  pullback::Hecke.Map

  function SpecOpenMor(
      U::DomainType,
      V::CodomainType,
      f::Vector{<:AbsSpecMor};
      check::Bool=true
    ) where {DomainType<:SpecOpen, CodomainType<:SpecOpen}
    Y = ambient(V)
    n = length(f)
    n == length(affine_patches(U)) || error("number of patches does not coincide with the number of maps")
    if check
      for i in 1:n
        domain(f[i]) === affine_patches(U)[i] || error("domain of definition of the map does not coincide with the patch")
        codomain(f[i]) === Y || error("codomain is not compatible")
      end
      for i in 1:n-1
	for j in i+1:n
	  A = intersect(domain(f[i]), domain(f[j]))
	  restrict(f[i], A, Y) == restrict(f[j], A, Y) || error("maps don't glue")
	end
      end
      for g in f
        is_empty(subscheme(domain(g), pullback(g).(gens(V)))) || error("image is not contained in the codomain")
      end
    end
    return new{DomainType, CodomainType}(U, V, f)
  end
end

