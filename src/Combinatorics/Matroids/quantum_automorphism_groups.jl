

function _index_number(i::Int)
  dgs = reverse(digits(i))
  return join(["₀₁₂₃₄₅₆₇₈₉"[3*d+1] for d in dgs], "") 
end

function _magic_unitary_symbols(n::Int=4)
  u = Matrix{String}(undef, n, n)
  for i in 1:n
    for j in 1:n
      if i > 9 || j > 9
        u[i, j] = "u$(_index_number(i))₋$(_index_number(j))"
      else
        u[i, j] = "u$(_index_number(i))$(_index_number(j))"
      end
    end
  end
  return u
end

function _quantum_symmetric_group_groebner_basis(n::Int; u::Matrix{Generic.FreeAssociativeAlgebraElem{T}}) where T
  function rwel(k::Int, j::Int, h::Int=3, v::Int=3)
    return sum([u[2, k] * u[s, j] for s in h:n]; init=zero(T)) - sum([u[s, k] * u[1, j] for s in v:n]) + u[1, j] - u[2, k]
  end

  function rinj(k::Int, j::Int, h::Int=3, v::Int=3)
    n = size(u)[1]
    return sum([u[k, 2] * u[j, s] for s in h:n]) - sum([u[k, s] * u[j, 1] for s in v:n]) + u[j, 1] - u[k, 2]
  end

  function wel(i::Int, j::Int, k::Int) 
    @assert 1 <= i <= size(u)[1] && 1 <= j <= size(u)[2] && 1 <= k <= size(u)[2] "Indices out of bounds"
    return u[i, j] * u[i, k]
  end

  function inj(i::Int, j::Int, k::Int)
    @assert 1 <= i <= size(u)[1] && 1 <= j <= size(u)[2] && 1 <= k <= size(u)[1] "Indices out of bounds"
    return u[i, j] * u[k, j]
  end

  function ip(i::Int, j::Int)
    @assert 1 <= i <= size(u)[1] && 1 <= j <= size(u)[2] "Indices out of bounds"
    return u[i, j] * u[i, j] - u[i, j]
  end
  function row_sum(i::Int)
    return sum([u[i, x] for x in 1:size(u)[1]]) - one(parent(u[1, 1]))
  end
  function col_sum(i::Int)
    return sum([u[x, i] for x in 1:size(u)[1]]) - one(parent(u[1, 1]))
  end

  function bg(z::Int, k::Int, j::Int, i::Int)
    z == 2  && return u[k, 2] * inj(j, 3, i, u=u) - rinj(k, j, u=u) * u[i, 3] 
    z == 8  && return u[2, k] * wel(3, j, i, u=u) - rwel(k, j, u=u) * u[3, i]
  end

  cs = [col_sum(i) for i in 1:n]
  rs = [row_sum(i) for i in 2:n]

  ips = [ip(i,j) for i in 2:n for j in 2:n]
  wels = [wel(i,j,k) for i in 2:n for j in 2:n for k in 2:n if j != k]
  injs = [inj(i,j,k) for i in 2:n for j in 2:n for k in 2:n if i != k]
  rwels = [rwel(k, j) for j in 2:n for k in 2:n if  k != j &&!(j == 3 && k == 2)]
  rinjs = [rinj(k, j) for j in 2:n for k in 2:n if k != j]

  e1 = [bg(2,k,j,i) for k=3:n for j=3:n for i=2:n if i!=j && j!=k]
  e1_s = [bg(8,k,j,i) for k=3:n for j=3:n for i=2:n if i!=j && j!=k]
  e2 = [bg(2,k,2,i) for k=3:n for i=4:n]

  e2_s = [bg(8,k,2,i) for k=3:n for i=4:n]
  e3 = [bg(2,2,j,i) for j in 5:n for i in 2:n if j!=i]
  e3_s = [bg(8,2,j,i) for j in 5:n for i in 2:n if j!=i]
  e4 = [bg(2,2,4,i) for i in 2:n if 4!=i && 3!=i]
  e4_s = [bg(8,2,4,i) for i in 2:n if 4!=i && 3!=i]
  e5 = [bg(2,2,4,3)]

  return vcat(cs, rs, ips, wels, injs, rwels, rinjs, e1, e1_s, e2, e2_s, e3, e3_s, e4, e4_s, e5)

end



@doc raw"""
    quantum_symmetric_group(n::Int; reduced_gb::Bool=true)

Return the ideal that defines the quantum symmetric group on `n` elements.
It is comprised of `2*n + n^2 + 2*n*n*(n-1)` many generators.
For `n >= 5`, a Gröbner basis is provided following the construction in [SW25](@cite).
If `reduced_gb` is set to `false`, the original generators are used to extend the Gröbner basis, which may improve reduction speed.

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
function quantum_symmetric_group(n::Int; reduced_gb::Bool=true)
  if Oscar.is_unicode_allowed()
    A, u = free_associative_algebra(QQ, _magic_unitary_symbols(n))
  else
    A, u = free_associative_algebra(QQ, :u => (1:n, 1:n))
  end

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
  
  I = ideal(relations)

  if n >= 5
    gb = _quantum_symmetric_group_groebner_basis(n; u=u)
    if !reduced_gb
      gb = vcat(gens(I), gb)
    end
    I.gb = IdealGens(gb)
  end

  return I
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

  if structure == :bases
    sets = bases(M)
  elseif structure == :circuits
    sets = circuits(M)
  elseif structure == :flats
    sets = flats(M)
  else
    error("unreachable")
  end

  #Computing the sizehint for the relations
  setSizes = map(length, sets)
  sizes = unique(setSizes)

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
  u = permutedims(reshape(gens(A),(n, n)))

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
  u = permutedims(reshape(gens(A),(n, n)))

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
