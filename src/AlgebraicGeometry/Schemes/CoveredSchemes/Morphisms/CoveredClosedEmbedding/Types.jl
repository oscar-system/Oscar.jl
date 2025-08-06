########################################################################
# Closed embeddings                                                    #
########################################################################

@doc raw"""
    CoveredClosedEmbedding <: AbsCoveredSchemeMorphism

Type for closed embeddings of covered schemes.
  
In addition to the closed embedding it stores 
the sheaf of ideals defining the image.
"""
@attributes mutable struct CoveredClosedEmbedding{
    DomainType<:AbsCoveredScheme,
    CodomainType<:AbsCoveredScheme,
    BaseMorphismType
   } <: AbsCoveredSchemeMorphism{
                                 DomainType,
                                 CodomainType,
                                 BaseMorphismType,
                                 CoveredSchemeMorphism
                                }
  f::CoveredSchemeMorphism
  I::AbsIdealSheaf

  function CoveredClosedEmbedding(
      X::DomainType,
      Y::CodomainType,
      f::CoveringMorphism{<:Any, <:Any, MorphismType, BaseMorType};
      check::Bool=true,
      ideal_sheaf::AbsIdealSheaf=IdealSheaf(Y, f, check=check)
    ) where {
             DomainType<:AbsCoveredScheme,
             CodomainType<:AbsCoveredScheme,
             MorphismType<:ClosedEmbedding,
             BaseMorType
            }
    ff = CoveredSchemeMorphism(X, Y, f; check)
    if has_decomposition_info(codomain(f))
      for U in patches(domain(f))
        floc = f[U]
        phi = pullback(floc)
        V = codomain(floc)
        g = Vector{elem_type(OO(V))}(decomposition_info(codomain(f))[V])
        set_decomposition_info!(domain(f), U, Vector{elem_type(OO(U))}(phi.(g)))
      end
    end
    #all(x->(x isa ClosedEmbedding), values(morphisms(f))) || error("the morphisms on affine patches must be `ClosedEmbedding`s")
    return new{DomainType, CodomainType, BaseMorType}(ff, ideal_sheaf)
  end
end
