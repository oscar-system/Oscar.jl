@doc raw"""
    quantum_symmetric_group(n::Int)

Return the ideal that defines the quantum symmetric group on `n` elements.
It is comprised of `2*n + n^2 + 2*n*n*(n-1)` many generators.

The relations are:

  - row and column sum relations: `2*n` relations
  - idempotent relations: `n^2` relations
  - relations of type `u[i,j]*u[i,k]` and `u[j,i]*u[k,i]` for `k != j`: `2*n*n*(n-1)` relations

# Examples

```jldoctest
julia> S4 = quantum_symmetric_group(4);

julia> length(gens(S4))
120
```
"""
function quantum_symmetric_group(n::Int)
  A, u = free_associative_algebra(QQ, :u => (1:n, 1:n))

  relations = elem_type(A)[]

  #Idempotent relations
  for i in 1:n, j in 1:n
    new_relation = u[i,j] * u[i,j] - u[i,j]
    push!(relations, new_relation)
    for k in 1:n
      if k != j
        new_relation = u[i,j] * u[i,k]
        push!(relations, new_relation)
        new_relation = u[j,i] * u[k,i]
        push!(relations, new_relation)
      end
    end
  end

  #row and column sum relations
  for i in 1:n
    push!(relations, sum(u[i,:]) - 1)
    push!(relations, sum(u[:,i]) - 1)

  end
  return ideal(relations)
end

@doc raw"""
    _quantum_automorphism_group_indices(M::Matroid, structure::Symbol=:bases)

Return the indices of the relations that define the quantum automorphism group of a matroid for a given structure.

# Examples

```jldoctest
julia> M = uniform_matroid(3,4);

julia> idx = Oscar._quantum_automorphism_group_indices(M,:bases);

julia> length(idx)
1920
```
"""
function _quantum_automorphism_group_indices(M::Matroid, structure::Symbol=:bases)
  @req structure in [:bases, :circuits, :flats] "structure must be one of :bases, :circuits, :flats"

  matroid_groundset(M) isa Vector{Int} ? nothing : M=Matroid(Oscar.pm_object(M))  

  n = length(M)
  grdSet = matroid_groundset(M)

  rels = Vector{Tuple{Int,Int}}[]

  sets  = getproperty(Oscar, structure)(M)
  sizes = unique!(map(length, sets))


  #Computing the sizehint for the relations
  setSizes = map(length, sets)

  sets_count = [0 for _ in 1:maximum(setSizes)]
  foreach(size->sets_count[size] += 1, setSizes)
  for i in 1:length(sets_count)
    if sets_count[i] == 0 || sets_count[i] == n^i
      sets_count[i] = 0
      continue
    end
    nonsets_count = n^i - (sets_count[i] * factorial(i))
    sets_count[i] = sets_count[i] * factorial(i) * nonsets_count * 2
  end

  vector_size = sum(sets_count)
  sizehint!(rels, vector_size)

  for size in sizes 
    size == 0 && continue

    S = symmetric_group(size)
    permutedSets = [[set[g(i)] for i in 1:size] for g in S for set in sets if length(set) == size]

    nonSets = MultiPartitionIterator{Int}(grdSet, size; toskip=permutedSets)

    for nonset in nonSets
      for set in permutedSets
        push!(rels, [(set[i],nonset[i]) for i in 1:size])
        push!(rels, [(nonset[i],set[i]) for i in 1:size])
      end
    end
  end
  @assert length(rels) == vector_size "The number of relations is not correct $(length(rels)) != $(vector_size)"

  return Vector{Vector{Tuple{Int,Int}}}(rels)
end

@doc raw"""
    quantum_automorphism_group(M::Matroid, structure::Symbol=:bases)

Return the ideal that defines the quantum automorphism group of a matroid for a given structure.

# Examples

```jldoctest
julia> G = complete_graph(4)
Undirected graph with 4 nodes and the following edges:
(2, 1)(3, 1)(3, 2)(4, 1)(4, 2)(4, 3)

julia> M = cycle_matroid(G);

julia> qAut = quantum_automorphism_group(M,:bases);

julia> length(gens(qAut))
23448
```
"""
function quantum_automorphism_group(
  M::Matroid,
  structure::Symbol=:bases)

  relation_indices = _quantum_automorphism_group_indices(M, structure)
  n = length(M)

  SN = quantum_symmetric_group(n)
  A = base_ring(SN)
  u = permutedims(reshape(gens(A),(n, n)),[2, 1])

  new_relations = elem_type(A)[prod(gen -> u[gen[1], gen[2]], relation; init=one(A)) for relation in relation_indices]
  return ideal(vcat(gens(SN), new_relations))

end

@doc raw"""
    quantum_automorphism_group(G::Graph{Undirected})

Return the ideal that defines the quantum automorphism group of the undirected graph `G`.

# Examples
```jldoctest
julia> G = graph_from_edges([[1, 2], [2, 4]]);

julia> qAut = quantum_automorphism_group(G);

julia> length(gens(qAut))
184
```
"""
function quantum_automorphism_group(G::Graph{Undirected})
  n = nv(G)

  SN = quantum_symmetric_group(n)
  A = base_ring(SN)
  u = permutedims(reshape(gens(A),(n, n)),[2, 1])

  nonedges = Tuple{Int,Int}[(i, j) for i in 1:n for j in i:n if !has_edge(G,i,j)]

  edgs = map(edg -> (src(edg), dst(edg)), edges(G)) # This should be |E| many edges

  new_relations = copy(gens(SN))
  sizehint!(new_relations, ngens(SN) + length(edgs) * length(nonedges)*8)

  function _addrelations(edge,nonedg)
    r1 = u[edge[1], nonedg[1]] * u[edge[2], nonedg[2]]
    r2 = u[edge[2], nonedg[1]] * u[edge[1], nonedg[2]]
    r3 = u[edge[1], nonedg[2]] * u[edge[2], nonedg[1]] 
    r4 = u[edge[2], nonedg[2]] * u[edge[1], nonedg[1]]
    push!(new_relations, r1, r2, r3, r4)
  end
  for edge in edgs, nonedg in nonedges
    _addrelations(edge,nonedg)
    _addrelations(nonedg,edge)
  end
  unique!(new_relations)
  return ideal(new_relations)
end

struct MultiPartitionIterator{T}
  elements::Vector{T}
  size::Int
  toskip::Vector{Vector{T}}
  function MultiPartitionIterator{T}(elements::Vector{T}, size::Int; toskip::Vector{Vector{T}}=Vector{Vector{T}}()) where T
    @req map(length, toskip) == fill(size, length(toskip)) "The size of the subsets to skip must be the same as the size of the subsets"
    @req all(x->x in elements, reduce(vcat, toskip)) "The elements to skip must be in the elements"
    new(elements, size, toskip)
  end
end

function Base.iterate(S::MultiPartitionIterator{T}, state=0) where T
  n = length(S.elements)
  while true
    if state > n^S.size-1
      return nothing
    else
      subset = [S.elements[div(state, n^(i-1)) % n + 1] for i in 1:S.size]
      if subset in S.toskip
        state += 1
        continue
      else
        return (subset, state + 1)
      end
    end
  end
end

Base.length(S::MultiPartitionIterator{T}) where T = length(S.elements)^S.size - length(S.toskip)
Base.IteratorSize(::Type{MultiPartitionIterator{T}}) where T = Base.HasLength()
Base.IteratorEltype(::Type{MultiPartitionIterator{T}}) where T = Base.HasEltype()
Base.eltype(::Type{MultiPartitionIterator{T}}) where T = Vector{T}
