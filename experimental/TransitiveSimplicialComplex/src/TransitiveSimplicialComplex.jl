# Add your new types, functions, and methods here.

struct TransitiveSimplicialComplex
  complex::SimplicialComplex
  dim::Int
  homology::Vector{FinGenAbGroup}
  aut_group::PermGroup
  topological_type::String

  TransitiveSimplicialComplex(K::SimplicialComplex, top_type::String) = new(
    K, dim(K), [homology(K, i) for i in 1:dim(K)], automorphism_group(K), top_type
  )

  function TransitiveSimplicialComplex(K::SimplicialComplex,
                                       aut_group::PermGroup,
                                       top_type::String)
    # TODO introduce a check that aut_group is isomorphic to vertex action
    new(K, dim(K), [homology(K, i) for i in 1:dim(K)], aut_group, top_type)
  end
end

function transitive_simplicial_complex(facets::Vector{Vector{Int}}, top_type::String)
  K = simplicial_complex(facets)
  TransitiveSimplicialComplex(K, top_type)
end

simplicial_complex(tsc::TransitiveSimplicialComplex) = tsc.complex
automorphism_group(tsc::TransitiveSimplicialComplex) = tsc.aut_group
dim(tsc::TransitiveSimplicialComplex) = tsc.dim
homology(tsc::TransitiveSimplicialComplex, i::Int) = tsc.homology[i]
topological_type(tsc) = tsc.topological_type



# julia> manifold_2_4_5_1=[[1,2,3],[1,2,4],[1,3,4],[2,3,4]]
# 4-element Vector{Vector{Int64}}:
#  [1, 2, 3]
#  [1, 2, 4]
#  [1, 3, 4]
#  [2, 3, 4]
# 
# julia> tsc = Oscar.transitive_simplicial_complex(manifold_2_4_5_1, "( + ; 0 ) = S^2")
# Oscar.TransitiveSimplicialComplex(Abstract simplicial complex of dimension 2 on 4 vertices, 2, FinGenAbGroup[Z/1, Z], Permutation group of degree 4, "( + ; 0 ) = S^2")
# 
# julia> save("/path/to/file.mrdi", tsc)
#
# julia> load("/path/to/file.mrdi")
# Oscar.TransitiveSimplicialComplex(Abstract simplicial complex of dimension 2 on 4 vertices, 2, FinGenAbGroup[Z/1, Z], Permutation group of degree 4, "( + ; 0 ) = S^2")
# 

@register_serialization_type TransitiveSimplicialComplex

type_params(tsc::TransitiveSimplicialComplex) = TypeParams(
  TransitiveSimplicialComplex,
  automorphism_group(tsc)
)

function save_object(s::SerializerState, tsc::TransitiveSimplicialComplex)
  save_data_dict(s) do 
    save_object(s, simplicial_complex(tsc), :simplicial_complex)
    save_object(s, topological_type(tsc), :topological_type)
  end
end

function load_object(s::DeserializerState, ::Type{TransitiveSimplicialComplex},
                     params::PermGroup)
  load_node(s) do _
    TransitiveSimplicialComplex(
      load_object(s, SimplicialComplex, :simplicial_complex),
      params,
      load_object(s, String, :topological_type)
    )
  end
end
