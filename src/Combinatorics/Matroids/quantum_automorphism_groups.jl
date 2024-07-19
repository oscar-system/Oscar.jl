@doc raw"""
    quantum_symmetric_group(n::Int)

Get the Ideal that defines the quantum symmetric group on `n` elements. It is comprised of `2*n + n^2 + 2*n*n*(n-1)` many generators.

# The relations are:

  - row and column sum relations: `2*n` relations
  - idempotent relations: `n^2` relations
  - relations of type `u[i,j]*u[i,k]` and `u[j,i]*u[k,i]` for `k != j`: `2*n*n*(n-1)` relations

# Output


# Example

```jldoctest
S4 = quantum_symmetric_group(4)
length(gens(S4))
# output

120
```


"""
function quantum_symmetric_group(n::Int)
  generator_strings = String[]

  for i in 1:n, j in 1:n
    push!(generator_strings, "u[$i,$j]")
  end

  A, g = free_associative_algebra(Oscar.QQ, generator_strings)
  u = Matrix{elem_type(A)}(undef, n, n)
  for i in 1:n, j in 1:n
    u[i, j] = g[(i-1)*n + j]
  end

  relations = elem_type(A)[]

  #Idempotent relations
  for i in 1:n, j in 1:n
    new_relation = u[i, j] * u[i, j] - u[i, j]
    push!(relations, new_relation)
    for k in 1:n
      if k != j
        new_relation = u[i,j] * u[i, k]
        push!(relations, new_relation)
        new_relation = u[j, i]*u[k, i]
        push!(relations, new_relation)
      end
    end
  end

  #row and column sum relations
  for i in 1:n
    new_relation_row = -1
    new_relation_col = -1
    for k in 1:n
      new_relation_row += u[i,k]
      new_relation_col += u[k,i]
    end
    push!(relations, new_relation_row)
    push!(relations, new_relation_col)

  end
  return ideal(relations)
end

@doc raw"""
    _quantum_automorphism_group_indices(M::Matroid, structure::Symbol=:bases)

Get the indices of the relations that define the quantum automorphism group of a matroid for a given structure.

# Examples
```jldoctest
M = uniform_matroid(3,4)
bases(M)
idx = Oscar._quantum_automorphism_group_indices(M,:bases)
length(idx)

# output

1920
"""
function _quantum_automorphism_group_indices(M::Matroid, structure::Symbol=:bases)
  @req structure in [:bases, :circuits, :flats] "structure must be one of :bases, :circuits, :flats"

  n = length(M)
  grdSet = matroid_groundset(M)

  b    = [[] for _ in 1:n]
  nb   = [[] for _ in 1:n]
  rels = [[] for _ in 1:n]

  sets  = eval(structure)(M)
  sizes = unique(map(x -> length(x), sets))
  tempGrdSet = reduce(vcat,[grdSet for i in 1:n])

  for size in sizes 
    size == 0 && continue
    powerSet = unique(sort.(Oscar.subsets(tempGrdSet,size)))

    setsOfSize = filter(x->length(x)==size,sets)
    nonSets = setdiff(powerSet,setsOfSize) 

    for set in setsOfSize
      append!(b[size],collect(Oscar.permutations(set)))
    end
    for nonset in nonSets
      append!(nb[size],collect(Oscar.permutations(nonset)))
    end
    for set in b[size]
      for nonset in nb[size]
        rel = []
        for i in 1:size
          push!(rel,(set[i],nonset[i]))
        end
        push!(rels[size],rel)
        rel = []
        for i in 1:size
          push!(rel,(nonset[i],set[i]))
        end
        push!(rels[size],rel)
      end
    end
  end

  rels = unique.(rels)
  return Vector{Vector{Tuple{Int,Int}}}(reduce(vcat,rels))
end

@doc raw"""
    quantum_automorphism_group(M::Matroid, structure::Symbol=:bases)

Get the Ideal that defines the quantum automorphism group of a matroid for a given structure.

# Examples

```jldoctest
M = uniform_matroid(3,4);
qAut = quantum_automorphism_group(M,:bases);
length(gens(qAut))
# output

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

  new_relations = []

  for relation in relation_indices
    temp = one(A)
    for gen in relation
      temp = temp * u[gen[1], gen[2]]
    end
    push!(new_relations,temp)
  end
  return ideal([gens(SN)..., new_relations...])

end

@doc raw"""
    quantum_automorphism_group(G::Graph)

Get the Ideal that defines the quantum automorphism group of a graph.

# Examples
```jldoctest
G = graph_from_edges([[1,2],[2,4]])
qAut = quantum_automorphism_group(G)
gens(qAut)

# output

168
```
"""
function quantum_automorphism_group(G::Graph)
  @req typeof(G) == Graph{Undirected} "Graph must be undirected"

  n = nv(G)

  SN = quantum_symmetric_group(n)
  A = base_ring(SN)
  u = permutedims(reshape(gens(A),(n,n)),[2,1])

  nonedges = Tuple{Int,Int}[] # This should be n^2 - |E| many nonedges
  for i in 1:n 
    for j in i:n
      if !has_edge(G,i,j)
        push!(nonedges,(i,j))
      end
    end
  end


  edgs = map(edg -> (src(edg), dst(edg)), collect(edges(G))) # This should be |E| many edges

  new_relations = []
  function _addrelations(edge,nonedg)
    r1 = u[edge[1],nonedg[1]] * u[edge[2],nonedg[2]]
    r2 = u[edge[2],nonedg[1]] * u[edge[1],nonedg[2]]
    r3 = u[edge[1],nonedg[2]] * u[edge[2],nonedg[1]] 
    r4 = u[edge[2],nonedg[2]] * u[edge[1],nonedg[1]]
    append!(new_relations,[r1,r2,r3,r4])
    return;
  end
  for edg in edgs
    for nonedg in nonedges
      _addrelations(edg,nonedg)
    end
  end
  return ideal([gens(SN)..., unique(new_relations)...])
end
