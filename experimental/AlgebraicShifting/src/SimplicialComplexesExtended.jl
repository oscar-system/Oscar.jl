using DataStructures
""" groupby(f, xs)

  Returns a dictionary Dict(i => {x ∈ xs | f(x) = i}).
"""
function groupby(f, xs)
  mergewith(
    union,
    Dict(),
    (Dict(f(x) => Set([x])) for x in xs)...
  )
end

""" complex_faces_by_dimension(K::SimplicialComplex, mindim::Int=0)

  Returns a Vector{Set{Set{Int}}} fs such that fs[d+1] contains all d-simplices of K.
"""
function complex_faces_by_dimension(K :: SimplicialComplex; mindim::Int=0)
  fcts = DefaultDict(Set{Set{Int}}(), groupby(length, facets(K))) # group facets by dimension+1
  fs = Vector{Set{Set{Int}}}(undef, dim(K)+1) # faces by dimension+1
  fs[dim(K)+1] = fcts[dim(K)+1]
  for d in dim(K)-1:-1:mindim # obtain d-simplices
    # println(fcts[d+1])
    fs[d+1] = union(fcts[d+1], subsets.(fs[d+2], d+1)...)
  end
  return fs
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
DefaultDict{Int64, Set{Set{Int64}}, Set{Set{Int64}}}(2 => Set([Set([4, 3])]), 3 => Set([Set([2, 3, 1]), Set([5, 2, 1])]))
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
DefaultDict{Int64, Set{Set{Int64}}, Set{Set{Int64}}}(2 => Set([Set([4, 3])]), 3 => Set([Set([2, 3, 1]), Set([5, 2, 1])]))
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

