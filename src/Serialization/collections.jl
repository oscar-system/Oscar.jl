@register_serialization_type LeechPair
type_params(x::LeechPair) = TypeParams(LeechPair, group(x))

function save_object(s::SerializerState, LG::LeechPair)
  save_data_dict(s) do
    save_object(s, Oscar.number(LG), :id)
    save_object(s, Oscar.rank_invariant_lattice(LG), :rank)
    save_object(s, Oscar.order_group(LG), :order)
    save_object(s, Oscar.alpha(LG), :alpha)
    save_object(s, Oscar.index_image_in_Oq_coinvariant(LG), :icoinvbar)
    save_object(s, Oscar.index_image_in_Oq_invariant(LG), :iinvbar)
    save_object(s, Oscar.index_normaliser_modulo_group(LG), :ind)
    save_object(s, Oscar.class_number_invariant_lattice(LG), :hinv)
    save_object(s, Oscar.number_of_niemeier_embeddings(LG), :N)
    # want to avoid using type key here to avoid confusion with file format
    save_object(s, LG.type, :leech_pair_type) 
  end
end

function load_object(s::DeserializerState, ::Type{LeechPair}, G::MatrixGroup)
  # dev until we move collection to offical db
  db = Oscar.Oscardb.get_db(;dev=true)
  leech = Oscar.Oscardb.find_one(db["zzlattices"], Dict("_id" => "leech"))

  LeechPair(
    leech,
    G,
    load_object(s, Int, :id),
    load_object(s, Int, :rank),
    load_object(s, ZZRingElem, :order),
    load_object(s, Int, :alpha),
    load_object(s, Int, :icoinvbar),
    load_object(s, Int, :iinvbar),
    load_object(s, Int, :ind),
    load_object(s, Int, :hinv),
    load_object(s, Int, :N), 
    load_object(s, String, :leech_pair_type),
  )
end

################################################################################
# TransitiveSimplicialComplex

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
