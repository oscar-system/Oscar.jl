################################################################################
# TransitiveSimplicialComplex

# Add your new types, functions, and methods here.

struct TransitiveSimplicialComplex
  _id::String
  complex::SimplicialComplex
  dim::Int
  n_vertices::Int
  f_vector::Vector{Int}
  betti_numbers::Vector{Int}
  aut_group::PermGroup
  topological_type::String # this may not be there, allow Nothing?

  TransitiveSimplicialComplex(args...) = new(args...)
  
  TransitiveSimplicialComplex(
    name::String,
    K::SimplicialComplex,
    top_type::String
  ) = new(
      name,
      K,
      dim(K),
      n_vertices(K),
      f_vector(K),
      betti_numbers(K),
      automorphism_group(K),
      top_type
  )
  
  TransitiveSimplicialComplex(
    name::String,
    K::SimplicialComplex,
    aut_group::PermGroup,
    top_type::String
  ) = new(
    name,
    K,
    dim(K),
    n_vertices(K),
    f_vector(K),
    betti_numbers(K),
    aut_group,
    top_type
  )
end

function transitive_simplicial_complex(name::String,
                                       facets::Vector{Vector{Int}},
                                       top_type::String)
  K = simplicial_complex(facets)
  TransitiveSimplicialComplex(name, K, top_type)
end

# these functions should be described in a more detail documentation to come
simplicial_complex(tsc::TransitiveSimplicialComplex) = tsc.complex
automorphism_group(tsc::TransitiveSimplicialComplex) = tsc.aut_group
dim(tsc::TransitiveSimplicialComplex) = tsc.dim
n_vertices(tsc::TransitiveSimplicialComplex) = tsc.n_vertices
f_vector(tsc::TransitiveSimplicialComplex) = tsc.f_vector
betti_numbers(tsc::TransitiveSimplicialComplex) = tsc.betti_numbers
homology(tsc::TransitiveSimplicialComplex, i::Int) = homology(simplicial_complex(tsc), i)
topological_type(tsc::TransitiveSimplicialComplex) = tsc.topological_type
name(tsc::TransitiveSimplicialComplex) = tsc._id

function Base.show(io::IO, tsc::TransitiveSimplicialComplex)
  println(io, "TransitiveSimplicialComplex: $(name(tsc))")
end
