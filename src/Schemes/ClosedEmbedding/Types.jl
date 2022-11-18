export ClosedEmbedding

########################################################################
# Special Type for closed embeddings of affine schemes                 #
########################################################################
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

  # On an affine scheme X return the subscheme defined by I ⊂ 𝒪(X).
  function ClosedEmbedding(X::AbsSpec, I::Ideal)
    base_ring(I) == OO(X) || error("ideal does not belong to the correct ring")
    Y = subscheme(X, I)
    inc = SpecMor(Y, X, hom(OO(X), OO(Y), gens(OO(Y))))
    return new{typeof(Y), typeof(X), pullback_type(inc)}(inc, I)
  end

  # Turn a dummy SpecMor for a canonical inclusion into a ClosedEmbedding. 
  # This requires specifying an ideal I such that the domain of the 
  # map is naturally isomorphic to the subscheme defined by that ideal.
  function ClosedEmbedding(f::SpecMor, I::Ideal; check::Bool=true)
    Y = domain(f)
    X = codomain(f)
    base_ring(I) == OO(X) || error("ideal does not belong to the correct ring")
    if check
      Y == subscheme(X, I)
      pullback(f).(gens(OO(X))) == gens(OO(Y))
    end
    return new{typeof(Y), typeof(X), pullback_type(f)}(f, I)
  end
end

