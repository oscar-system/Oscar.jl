""" complex_faces_by_dimension(K::SimplicialComplex, mindim::Int=0)

  Returns a Vector{Set{Set{Int}}} fs such that fs[d+1] contains all d-simplices of K.
"""
function complex_faces_by_dimension(K :: SimplicialComplex; mindim::Int=0)
  po = face_poset(K)
  return map(k -> Set(Set{Int}.(
    data.(
      elements_of_rank(po, k))
  )), 1:rank(po))
end

@doc raw"""
complex_faces(K::SimplicialComplex)

Returns a Set{Set{Int}} containing all simplices of K.

complex_faces(K::SimplicialComplex, dim::Int)

Returns a Set{Set{Int}} containing all d-simplices of K.
"""
complex_faces(K :: SimplicialComplex) = union(complex_faces_by_dimension(K)...) :: Set{Set{Int}}
complex_faces(K :: SimplicialComplex, dim::Int) = complex_faces_by_dimension(K, mindim=dim)[dim+1] :: Set{Set{Int}}

