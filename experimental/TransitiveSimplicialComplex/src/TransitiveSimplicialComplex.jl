# Add your new types, functions, and methods here.

struct TransitiveSimplicialComplex
  K::SimplicialComplex
  dim::Int
  homology::Vector{FinGenAbGroup}
  topological_type::String
  aut_group::PermGroup
end

function transitive_simplicial_complex(facets::Vector{Vector{Int}}, top_type::String)
  K = simplicial_complex(facets)

  return TransitiveSimplicialComplex(
    K, dim(K), [homology(K, i) for i in 1:dim(K)], automorphism_group(K), top_type
  )
end

