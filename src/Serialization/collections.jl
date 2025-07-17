################################################################################
# Leech Pair
struct LeechPair
  leech::ZZLat
  G::MatrixGroup{QQFieldElem, QQMatrix}
  _id::Int
  rank::Int 
  order::ZZRingElem
  alpha::Int
  icoinvbar::Int
  iinvbar::Int
  ind::Int
  hinv::Int
  N::Int
  type::String

  function LeechPair(id::Int, leech::ZZLat, gg::Vector{QQMatrix})
    l = readlines("/homes/combi/vecchia/local/repositories/FinGrpTorelliVar/HKData/LeechPairs/query.csv")[id]
    vl = split(l, "|")
    return new(
      leech,
      matrix_group(gg),
      id,
      parse(Int, vl[2]),
      Hecke.parse(ZZRingElem, vl[3]),
      parse(Int, vl[4]),
      parse(Int, vl[5]),
      parse(Int, vl[6]),
      parse(Int, vl[7]),
      parse(Int, vl[8]),
      parse(Int, vl[9]),
      vl[10]
    )
  end
end

###############################################################################
#
# Accessors and methods
#
###############################################################################

Oscar.lattice(LG::LeechPair) = LG.leech

Oscar.group(LG::LeechPair) = LG.G

group_description(LG::LeechPair) = describe(LG.G)

number(LG::LeechPair) = LG._id

rank_invariant_lattice(LG::LeechPair) = LG.rank

order_group(LG::LeechPair) = LG.order

alpha(LG::LeechPair) = LG.alpha

index_image_in_Oq_coinvariant(LG::LeechPair) = LG.icoinvbar

index_image_in_Oq_invariant(LG::LeechPair) = LG.iinvbar

index_normaliser_modulo_group(LG::LeechPair) = LG.ind

class_number_invariant_lattice(LG::LeechPair) = LG.hinv

number_of_niemeier_embeddings(LG::LeechPair) = LG.N

Oscar.type(LG::LeechPair) = LG.type

Oscar.invariant_lattice(LG::LeechPair) = invariant_lattice(lattice(LG), group(LG))

Oscar.coinvariant_lattice(LG::LeechPair) = coinvariant_lattice(lattice(LG), group(LG))

Oscar.invariant_coinvariant_pair(LG::LeechPair) = invariant_coinvariant_pair(lattice(LG), group(LG))

###############################################################################
#
# Printing
#
###############################################################################

function Base.show(io::IO, LG::LeechPair)
  print(io, "Leech pair no. $(number(LG))")
end

function Base.show(io::IO, ::MIME"text/plain", LG::LeechPair)
  if get(io, :supercompact, false)
    print(io, "Leech pair")
  else
    println(io, "Leech pair no. $(number(LG))")
    println(io, "Rank:  $(rank_invariant_lattice(LG))")
    println(io, "Order: $(order_group(LG))")
    print(io, "Type:  $(type(LG))")
  end
end

@register_serialization_type LeechPair
type_params(x::LeechPair) = TypeParams(LeechPair, group(x))

function save_object(s::SerializerState, ::LeechPair)
  save_data_dict(s) do
    save_object(s, group(), :group)
    save_object(s, number(), :_id)
    save_object(s, rank_invariant_lattice(), :rank)
    save_object(s, index_image_in_Oq_coinvariant(), :icoinvbar)
    save_object(s, index_image_in_Oq_invariant(), :iinvbar)
    save_object(s, index_normaliser_modulo_group(), :hinv)
    save_object(s, number_of_niemeier_embeddings(), :N)
    # want to avoid using type key here to avoid confusion with file format
    save_object(s, number_of_niemeier_embeddings(), :leech_pair_type) 
  end
end

struct LeechPair
  leech::ZZLat
  G::MatrixGroup{QQFieldElem, QQMatrix}
  _id::Int
  rank::Int 
  order::ZZRingElem
  alpha::Int
  icoinvbar::Int
  iinvbar::Int
  ind::Int
  hinv::Int
  N::Int
  type::String
end

###############################################################################
#
# Accessors and methods
#
###############################################################################

Oscar.lattice(LG::LeechPair) = LG.leech

Oscar.group(LG::LeechPair) = LG.G

group_description(LG::LeechPair) = describe(LG.G)

number(LG::LeechPair) = LG._id

rank_invariant_lattice(LG::LeechPair) = LG.rank

order_group(LG::LeechPair) = LG.order

alpha(LG::LeechPair) = LG.alpha

index_image_in_Oq_coinvariant(LG::LeechPair) = LG.icoinvbar

index_image_in_Oq_invariant(LG::LeechPair) = LG.iinvbar

index_normaliser_modulo_group(LG::LeechPair) = LG.ind

class_number_invariant_lattice(LG::LeechPair) = LG.hinv

number_of_niemeier_embeddings(LG::LeechPair) = LG.N

Oscar.type(LG::LeechPair) = LG.type

Oscar.invariant_lattice(LG::LeechPair) = invariant_lattice(lattice(LG), group(LG))

Oscar.coinvariant_lattice(LG::LeechPair) = coinvariant_lattice(lattice(LG), group(LG))

Oscar.invariant_coinvariant_pair(LG::LeechPair) = invariant_coinvariant_pair(lattice(LG), group(LG))

###############################################################################
#
# Printing
#
###############################################################################

function Base.show(io::IO, LG::LeechPair)
  print(io, "Leech pair no. $(number(LG))")
end

function Base.show(io::IO, ::MIME"text/plain", LG::LeechPair)
  if get(io, :supercompact, false)
    print(io, "Leech pair")
  else
    println(io, "Leech pair no. $(number(LG))")
    println(io, "Rank:  $(rank_invariant_lattice(LG))")
    println(io, "Order: $(order_group(LG))")
    print(io, "Type:  $(type(LG))")
  end
end

@register_serialization_type LeechPair
type_params(x::LeechPair) = TypeParams(LeechPair, group(x))

function save_object(s::SerializerState, LG::LeechPair)
  save_data_dict(s) do
    save_object(s, group(LG), :group)
    save_object(s, number(LG), :_id)
    save_object(s, rank_invariant_lattice(LG), :rank)
    save_object(s, index_image_in_Oq_coinvariant(LG), :icoinvbar)
    save_object(s, index_image_in_Oq_invariant(LG), :iinvbar)
    save_object(s, index_normaliser_modulo_group(LG), :hinv)
    save_object(s, number_of_niemeier_embeddings(LG), :N)
    # want to avoid using type key here to avoid confusion with file format
    save_object(s, number_of_niemeier_embeddings(LG), :leech_pair_type) 
  end
end

function load_object(s::DeserializerState, ::Type{LeechPair}, G::MatrixGroup)
  # need to be able to load from oscar-dev 
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
