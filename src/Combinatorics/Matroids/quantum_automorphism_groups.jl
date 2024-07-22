@doc raw"""
    quantum_symmetric_group(n::Int)

Return the ideal that defines the quantum symmetric group on `n` elements.
It is comprised of `2*n + n^2 + 2*n*n*(n-1)` many generators.

The relations are:

  - row and column sum relations: `2*n` relations
  - idempotent relations: `n^2` relations
  - relations of type `u[i,j]*u[i,k]` and `u[j,i]*u[k,i]` for `k != j`: `2*n*n*(n-1)` relations

# Example

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

  n = length(M)
  grdSet = matroid_groundset(M)

  b    = [Vector{Int}[] for _ in 1:n]
  nb   = [Vector{Int}[] for _ in 1:n]
  rels = [Vector{Tuple{Int,Int}}[] for _ in 1:n]

  sets  = getproperty(Oscar, structure)(M)
  sizes = unique!(map(length, sets))
  tempGrdSet = reduce(vcat,[grdSet for i in 1:n])

  for size in sizes 
    size == 0 && continue
    powerSet = unique(sort.(subsets(tempGrdSet,size)))

    setsOfSize = filter(x->length(x)==size,sets)
    nonSets = setdiff(powerSet,setsOfSize) 

    S = symmetric_group(size)

    for set in setsOfSize
      append!(b[size],[[set[g(i)] for i in 1:size] for g in S])
    end
    for nonset in nonSets
      append!(nb[size],[[nonset[g(i)] for i in 1:size] for g in S])
    end
    for set in b[size]
      for nonset in nb[size]
        push!(rels[size], [(set[i],nonset[i]) for i in 1:size])
        push!(rels[size], [(nonset[i],set[i]) for i in 1:size])
      end
    end
  end
  rels = unique!.(rels)
  return Vector{Vector{Tuple{Int,Int}}}(reduce(vcat,rels))
end

@doc raw"""
    quantum_automorphism_group(M::Matroid, structure::Symbol=:bases)

Return the Ideal that defines the quantum automorphism group of a matroid for a given structure.

# Examples

```jldoctest
julia> M = uniform_matroid(3,4);

julia> qAut = quantum_automorphism_group(M,:bases);

julia> length(gens(qAut))
2040
```
"""
function quantum_automorphism_group(
  M::Matroid,
  structure::Symbol=:bases)

  @req matroid_groundset(M) isa Vector{Int} "Matroid must have integer groundset"

  relation_indices = _quantum_automorphism_group_indices(M,structure)
  n = length(M)

  SN = quantum_symmetric_group(n)
  A = base_ring(SN)
  u = permutedims(reshape(gens(A),(n,n)),[2,1])

  new_relations = elem_type(A)[prod(gen -> u[gen[1], gen[2]], relation; init=one(A)) for relation in relation_indices]
  return ideal([gens(SN)..., new_relations...])

end

@doc raw"""
    quantum_automorphism_group(G::Graph)

Return the Ideal that defines the quantum automorphism group of a graph.

# Examples
```jldoctest
julia> G = graph_from_edges([[1,2],[2,4]]);

julia> qAut = quantum_automorphism_group(G);

julia> length(gens(qAut))
168
```
"""
function quantum_automorphism_group(G::Graph)
  @req G isa Graph{Undirected} "Graph must be undirected"

  n = nv(G)

  SN = quantum_symmetric_group(n)
  A = base_ring(SN)
  u = permutedims(reshape(gens(A),(n,n)),[2,1])

  nonedges = Tuple{Int,Int}[(i, j) for i in 1:n for j in i:n if !has_edge(G,i,j)] # This should be n^2 - |E| many nonedges

  edgs = map(edg -> (src(edg), dst(edg)), edges(G)) # This should be |E| many edges

  function _addrelations(edge,nonedg)
    r1 = u[edge[1],nonedg[1]] * u[edge[2],nonedg[2]]
    r2 = u[edge[2],nonedg[1]] * u[edge[1],nonedg[2]]
    r3 = u[edge[1],nonedg[2]] * u[edge[2],nonedg[1]] 
    r4 = u[edge[2],nonedg[2]] * u[edge[1],nonedg[1]]
    return [r1,r2,r3,r4]
  end
  new_relations = reduce(vcat,[_addrelations(edg, nonedg) for edg in edgs for nonedg in nonedges])
  return ideal([gens(SN)..., unique(new_relations)...])
end
