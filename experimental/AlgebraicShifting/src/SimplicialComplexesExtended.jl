""" complex_faces_by_dimension(K::SimplicialComplex, mindim::Int=0)

  Returns a Vector{Set{Set{Int}}} fs such that fs[d+1] contains all d-simplices of K.
"""
function complex_faces_by_dimension(K :: SimplicialComplex; mindim::Int=0)
  po = face_poset(K)
  faces_by_dim = [Set{Set{Int}}() for _ in 1:rank(po)]
  HD = K.pm_simplicialcomplex.HASSE_DIAGRAM
  for i in 1:length(po)
    d = HD.DECORATION[i]
    Polymake.decoration_rank(d) < mindim && continue
    Polymake.decoration_face(d) == Set([-1]) && continue
    push!(faces_by_dim[Polymake.decoration_rank(d)],
          Polymake.to_one_based_indexing(Polymake.decoration_face(d)))
  end
  return faces_by_dim
end

@doc raw"""
complex_faces(K::SimplicialComplex)

Returns a Set{Set{Int}} containing all simplices of K.

complex_faces(K::SimplicialComplex, dim::Int)

Returns a Set{Set{Int}} containing all d-simplices of K.

# Examples
```jldoctest
julia> K = simplicial_complex([[1,2,3],[1,2,5],[3,4]])
Abstract simplicial complex of dimension 2 on 5 vertices

julia> complex_faces(K)
Set{Set{Int64}} with 13 elements:
  Set([4, 3])
  Set([5, 1])
  Set([5, 2, 1])
  Set([2, 3])
  Set([4])
  Set([3])
  Set([5])
  Set([1])
  Set([2, 3, 1])
  Set([2])
  Set([3, 1])
  Set([2, 1])
  Set([5, 2])

julia> complex_faces(K, 1)
Set{Set{Int64}} with 6 elements:
  Set([4, 3])
  Set([5, 1])
  Set([3, 1])
  Set([2, 3])
  Set([2, 1])
  Set([5, 2])```
"""
complex_faces(K :: SimplicialComplex) = union(complex_faces_by_dimension(K)...) :: Set{Set{Int}}
complex_faces(K :: SimplicialComplex, dim::Int) = complex_faces_by_dimension(K, mindim=dim)[dim+1] :: Set{Set{Int}}

