
########################################################################
# Special Type for closed embeddings of affine schemes                 #
########################################################################
@doc raw"""
    ClosedEmbedding{DomainType, CodomainType, PullbackType}

A closed embedding ``f : X → Y`` of affine schemes ``X = Spec(S)``
into ``Y = Spec(R)`` such that ``S ≅ R/I`` via ``f`` for some
ideal ``I ⊂ R``.
"""
@attributes mutable struct ClosedEmbedding{DomainType<:AbsAffineScheme,
                                           CodomainType<:AbsAffineScheme,
                                           PullbackType<:Map
                                          }<:AbsAffineSchemeMor{DomainType,
                                                        CodomainType,
                                                        PullbackType,
                                                        ClosedEmbedding,
                                                        Nothing
                                                       }
  inc::AffineSchemeMor{DomainType, CodomainType, PullbackType}
  I::Ideal
  U::AffineSchemeOpenSubscheme

  # On an affine scheme X return the subscheme defined by I ⊂ 𝒪(X).
  function ClosedEmbedding(X::AbsAffineScheme, I::Ideal; check::Bool=true)
    base_ring(I) === OO(X) || error("ideal does not belong to the correct ring")
    Y = subscheme(X, I)
    inc = morphism(Y, X, hom(OO(X), OO(Y), gens(OO(Y)); check=false); check=false)
    return new{typeof(Y), typeof(X), pullback_type(inc)}(inc, I)
  end

  # Turn a dummy AffineSchemeMor for a canonical inclusion into a ClosedEmbedding. 
  # This requires specifying an ideal I such that the domain of the 
  # map is naturally isomorphic to the subscheme defined by that ideal.
  function ClosedEmbedding(f::AffineSchemeMor, I::Ideal; check::Bool=true)
    Y = domain(f)
    X = codomain(f)
    base_ring(I) == OO(X) || error("ideal does not belong to the correct ring")
    @check Y == subscheme(X, I) "scheme is not compatible with ideal"
    @check pullback(f).(gens(OO(X))) == gens(OO(Y)) "variables are not preserved by the map"
    return new{typeof(Y), typeof(X), pullback_type(f)}(f, I)
  end
end

