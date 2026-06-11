struct UniformHypergraph
  n_vertices::Int
  k::Int
  faces::Vector{Vector{Int}}

  function UniformHypergraph(n_vertices::Int, k::Int, faces::Vector{Vector{Int}})
    @req all(length(f) == k && 1 <= minimum(f) && maximum(f) <= n_vertices for f in faces) "Parameters don't define a uniform hypergraph."
    return new(n_vertices, k, faces)
  end
end

const CollectionTypes = Union{Vector{Set{Int}}, Vector{Vector{Int}}, Vector{Combination{Int}}, Set{Set{Int}}}
uniform_hypergraph(faces::CollectionTypes, n::Int, k::Int) = UniformHypergraph(n, k, sort(unique(sort.(unique.(faces)))))
uniform_hypergraph(faces::CollectionTypes, n::Int) = uniform_hypergraph(faces, n, only(unique(length.(faces))))

function uniform_hypergraph(faces::CollectionTypes) 
  isempty(faces) && return UniformHypergraph(0, 0, Vector{Int}[])
  return uniform_hypergraph(faces, maximum(maximum.(faces)))
end

@doc raw"""
     uniform_hypergraph(faces, n::Int, k::Int)
     uniform_hypergraph(faces, n::Int)
     uniform_hypergraph(faces)
     uniform_hypergraph(K::SimplicialComplex, k::Int)

Create a uniform hypergraph using `faces`, which should be an iterable of size-k iterables of integers at most n.
One can also create a `UniformHypergraph` for the `k-1`-faces of a `SimplicialComplex` `K`.

#Examples
```jldoctest
julia> U = uniform_hypergraph([[1, 2], [2, 3]], 4)
UniformHypergraph(4, 2, [[1, 2], [2, 3]])

julia> U = uniform_hypergraph([[1, 2], [2, 3]])
UniformHypergraph(3, 2, [[1, 2], [2, 3]])

julia> U = uniform_hypergraph(simplicial_complex([[1 ,2, 3], [1, 3, 4]]), 2)
UniformHypergraph(4, 2, [[1, 2], [1, 3], [1, 4], [2, 3], [3, 4]])

julia> face_size(U)
2

julia> faces(U)
5-element Vector{Vector{Int64}}:
 [1, 2]
 [1, 3]
 [1, 4]
 [2, 3]
 [3, 4]
```
"""
uniform_hypergraph(K::SimplicialComplex, k::Int) = uniform_hypergraph(faces(K, k-1), n_vertices(K), k)
uniform_hypergraph(G::Graph{Undirected}) = uniform_hypergraph([[src(e), dst(e)] for e in edges(G)], n_vertices(G))
simplicial_complex(K::UniformHypergraph) = simplicial_complex([[[i] for i in 1:n_vertices(K)]; faces(K)])
n_vertices(K::UniformHypergraph) = K.n_vertices

faces(K::UniformHypergraph) = K.faces
# added for covenience when writting functions for Simplicial complex and Uniform Hypergraph
facets(K::UniformHypergraph) = Set.(K.faces)
face_size(K::UniformHypergraph) = K.k
Base.isempty(K::UniformHypergraph) = isempty(faces(K))

Base.hash(K :: UniformHypergraph, u :: UInt) = hash(K.n_vertices, hash(K.k, hash(K.faces, u)))
Base.:(==)(K :: UniformHypergraph, L :: UniformHypergraph) = K.n_vertices == L.n_vertices && K.k == L.k && K.faces == L.faces
@doc raw"""
     alexander_dual(K::UniformHypergraph)

Given a `UniformHypergraph` return the Alexander dual, seen as bijection ``\binom{[n]}{k} \to \binom{[n]}{n-k}``

# Examples
```jldoctest
julia> K = uniform_hypergraph([[1, 2], [2, 3], [3, 4]])
UniformHypergraph(4, 2, [[1, 2], [2, 3], [3, 4]])

julia> alexander_dual(K)
UniformHypergraph(4, 2, [[1, 3], [2, 3], [2, 4]])
```
"""
alexander_dual(K::UniformHypergraph) = uniform_hypergraph(alexander_dual(simplicial_complex(K)), n_vertices(K) - face_size(K))
is_shifted(K::UniformHypergraph) = is_shifted(simplicial_complex(K))
