# this might not be necessary if we can use polymake hasse diagram functionality?
function complex_faces(K :: SimplicialComplex, d :: Int) :: Vector{Vector{Int}}
  return sort(union(subsets.(sort.(collect.(facets(K))), d+1)...))
  # return union(subsets.(facets(K), d+1)...)
end

struct UniformHypergraph
  n_vertices::Int
  k::Int
  faces::Vector{Vector{Int}}

  function UniformHypergraph(n_vertices::Int, k::Int, faces::Vector{Vector{Int}})
    @req all(length(f) == k && 1 <= minimum(f) && maximum(f) <= n_vertices for f in faces) "Parameters don't define a uniform hypergraph."
    return new(n_vertices, k, faces)
  end
end

function uniform_hypergraph(faces::Vector{Vector{Int}}, n::Int, k::Int)
  faces = sort(collect(Set(sort.(collect.(Set.(faces))))))
  return UniformHypergraph(n, k, faces)
end

function uniform_hypergraph(faces::Vector{Vector{Int}}, n::Int)
  return uniform_hypergraph(faces, n, only(Set(length.(faces))))
end

function uniform_hypergraph(faces::Vector{Vector{Int}})
  return uniform_hypergraph(faces, maximum(maximum.(faces)))
end

@doc raw"""
     uniform_hypergraph(faces::Vector{Vector{Int}}, n::Int, k::Int)
     uniform_hypergraph(faces::Vector{Vector{Int}}, n::Int)
     uniform_hypergraph(faces::Vector{Vector{Int}})
     uniform_hypergraph(K::SimplicialComplex, k::Int)

Create a uniform hypergraph using `faces`, the size of each face should be `k` and all faces should be subsets of $[n]$.
One can also create a `UniformHypergraph` for the `k` faces of a `SimplicialComplex` `K`.

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
function uniform_hypergraph(K::SimplicialComplex, k::Int)
  return uniform_hypergraph(complex_faces(K, k-1), n_vertices(K), k)
end

function simplicial_complex(K::UniformHypergraph)
  return simplicial_complex([[[i] for i in 1:n_vertices(K)]; faces(K)])
end

n_vertices(K::UniformHypergraph) = K.n_vertices

faces(K::UniformHypergraph) = K.faces
# added for covenience when writting functions for Simplicial complex and Uniform Hypergraph
facets(K::UniformHypergraph) = Set.(K.faces)
face_size(K::UniformHypergraph) = K.k

function Base.hash(K :: UniformHypergraph, u :: UInt)
  return hash(K.n_vertices, hash(K.k, hash(K.faces, u)))
end

function Base.:(==)(K :: UniformHypergraph, L :: UniformHypergraph)
  return K.n_vertices == L.n_vertices && K.k == L.k && K.faces == L.faces
end

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
function alexander_dual(K::UniformHypergraph)
  return uniform_hypergraph(alexander_dual(simplicial_complex([[[i] for i in 1:K.n_vertices]; K.faces])), K.n_vertices - K.k)
end

is_shifted(K::UniformHypergraph) = is_shifted(simplicial_complex(K))
