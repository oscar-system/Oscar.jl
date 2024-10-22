# this might not be necessary if we can use polymake hasse diagram functionality?
function complex_faces(K :: SimplicialComplex, d :: Int) :: Vector{Vector{Int}}
  return sort(union(subsets.(sort.(collect.(facets(K))), d+1)...))
  # return union(subsets.(facets(K), d+1)...)
end

struct UniformHypergraph
  n_vertices :: Int
  k :: Int
  faces :: Vector{Vector{Int}}
end

function uniform_hypergraph(faces :: Vector{Vector{Int}}, n :: Int, k :: Int)
  faces = sort(collect(Set(sort.(collect.(Set.(faces))))))
  @req all(length(f) == k && 1 <= minimum(f) && maximum(f) <= n for f in faces) "Parameters don't define a uniform hypergraph."
  return UniformHypergraph(n, k, faces)
end

function uniform_hypergraph(faces :: Vector{Vector{Int}}, n :: Int)
  return uniform_hypergraph(faces, n, only(Set(length.(faces))))
end

function uniform_hypergraph(faces :: Vector{Vector{Int}})
  return uniform_hypergraph(faces, maximum(maximum.(faces)))
end

function uniform_hypergraph(K :: SimplicialComplex, k :: Int)
  return uniform_hypergraph(complex_faces(K, k-1), n_vertices(K), k)
end

import Oscar: n_vertices, simplicial_complex

function simplicial_complex(K :: UniformHypergraph)
  return simplicial_complex([[[i] for i in 1:n_vertices(K)]; faces(K)])
end

function n_vertices(K :: UniformHypergraph)
  return K.n_vertices  
end

faces(K::UniformHypergraph) = K.faces
face_size(K::UniformHypergraph) = K.k

function Base.hash(K :: UniformHypergraph, u :: UInt)
  return hash(K.n_vertices, hash(K.k, hash(K.faces, u)))
end

function Base.:(==)(K :: UniformHypergraph, L :: UniformHypergraph)
  return K.n_vertices == L.n_vertices && K.k == L.k && K.faces == L.faces
end

raw""" Alexander dual, seen as bijection ``\binom{[n]}{k} \to \binom{[n]}{n-k}`` """
function alexander_dual(K :: UniformHypergraph)
  return uniform_hypergraph(alexander_dual(simplicial_complex([[[i] for i in 1:K.n_vertices]; K.faces])), K.n_vertices - K.k)
end
