# -*- Undirected graphs -*-

@doc raw"""
    pairwise_markov(G::Graph{Undirected})

Return all *elementary* CI statements in the pairwise Markov property of `G`,
i.e., all `CIStmt([i],[j],K)` such that `i` and `j` are not adjacent and
`K = setdiff(vertices(G), [i,j])`.

## Examples

```jldoctest
julia> G = graph_from_edges(Undirected, [[1,2],[2,3],[1,4],[4,5]])
Undirected graph with 5 nodes and the following edges:
(2, 1)(3, 2)(4, 1)(5, 4)

julia> pairwise_markov(G)
6-element Vector{CIStmt}:
 [1 _||_ 3 | {2, 4, 5}]
 [1 _||_ 5 | {2, 3, 4}]
 [2 _||_ 4 | {1, 3, 5}]
 [2 _||_ 5 | {1, 3, 4}]
 [3 _||_ 4 | {1, 2, 5}]
 [3 _||_ 5 | {1, 2, 4}]
```
"""
function pairwise_markov(G::Graph{Undirected})::Vector{CIStmt}
  n = n_vertices(G)
  [ci_stmt(i, j, setdiff(1:n, [i,j])) for i in 1:n for j in (i+1):n
    if !has_edge(G, i, j)]
end

@doc raw"""
    local_markov(G::Graph{Undirected})

Return all *elementary* CI statements in the local Markov property of `G`,
i.e., all `CIStmt([i],[j],K)` such that `i` and `j` are not adjacent and
`K` is the set of neighbors of `i`.

## Examples

```jldoctest
julia> G = graph_from_edges(Undirected, [[1,2],[2,3],[1,4],[4,5]])
Undirected graph with 5 nodes and the following edges:
(2, 1)(3, 2)(4, 1)(5, 4)

julia> local_markov(G)
12-element Vector{CIStmt}:
 [1 _||_ 3 | {2, 4}]
 [1 _||_ 5 | {2, 4}]
 [2 _||_ 4 | {1, 3}]
 [2 _||_ 5 | {1, 3}]
 [1 _||_ 3 | 2]
 [3 _||_ 4 | 2]
 [3 _||_ 5 | 2]
 [2 _||_ 4 | {1, 5}]
 [3 _||_ 4 | {1, 5}]
 [1 _||_ 5 | 4]
 [2 _||_ 5 | 4]
 [3 _||_ 5 | 4]
```
"""
function local_markov(G::Graph{Undirected})::Vector{CIStmt}
  n = n_vertices(G)
  unique([ci_stmt(i, j, neighbors(G, i)) for i in 1:n for j in 1:n
          if i != j && !has_edge(G, i, j)])
end

@doc raw"""
    are_separated(G::Graph{Undirected}, i::Int, j::Int, K::Vector{Int})

Check if `i` and `j` are separated by `K` in the graph `G`, i.e., whether
every path connecting `i` and `j` contains a node from `K`.

## Examples

```jldoctest
julia> G = graph_from_edges(Undirected, [[1,2],[2,3],[2,4],[3,4]])
Undirected graph with 4 nodes and the following edges:
(2, 1)(3, 2)(4, 2)(4, 3)

julia> are_separated(G, 1, 4, [2])
true

julia> are_separated(G, 1, 4, [3])
false
```
"""
function are_separated(G::Graph{Undirected}, i::Int, j::Int, K::Vector{Int})::Bool
  @req i != j && !(i in K) && !(j in K) "i and j must be distinct and not in K"
  V = vertices(G)
  todo = [ [setdiff(V, i), [i]] ]
  while !isempty(todo)
    cur = pop!(todo)
    W = cur[1]
    p = cur[2]
    v = p[end]
    if v == j
      if length(intersect(p, K)) == 0
        return false
      end
    end
    append!(todo, [ [setdiff(W, w), [p..., w]] for w in neighbors(G, v) if w in W])
  end
  return true
end

@doc raw"""
    global_markov(G::Graph{Undirected})

Return all *elementary* CI statements in the global Markov property of `G`,
i.e., all `CIStmt([i],[j],K)` such that `K` separates `i` and `j` in `G`
(see `are_separated`).

## Examples

```jldoctest
julia> G = graph_from_edges(Undirected, [[1,2],[2,3],[1,4],[4,5]])
Undirected graph with 5 nodes and the following edges:
(2, 1)(3, 2)(4, 1)(5, 4)

julia> global_markov(G)
31-element Vector{CIStmt}:
 [1 _||_ 3 | 2]
 [1 _||_ 3 | {2, 4}]
 [1 _||_ 3 | {2, 5}]
 [1 _||_ 3 | {2, 4, 5}]
 [1 _||_ 5 | 4]
 [1 _||_ 5 | {2, 4}]
 [1 _||_ 5 | {3, 4}]
 [1 _||_ 5 | {2, 3, 4}]
 [2 _||_ 4 | 1]
 [2 _||_ 4 | {1, 3}]
 ⋮
 [3 _||_ 4 | {2, 5}]
 [3 _||_ 4 | {1, 2, 5}]
 [3 _||_ 5 | 1]
 [3 _||_ 5 | 2]
 [3 _||_ 5 | 4]
 [3 _||_ 5 | {1, 2}]
 [3 _||_ 5 | {1, 4}]
 [3 _||_ 5 | {2, 4}]
 [3 _||_ 5 | {1, 2, 4}]
```
"""
function global_markov(G::Graph{Undirected})::Vector{CIStmt}
  [stmt for stmt in ci_statements(collect(vertices(G)))
   if are_separated(G, stmt.I[1], stmt.J[1], stmt.K)]
end

# -*- Directed acyclic graphs -*-

parents(G::Graph{Directed}, i::Int)  =  inneighbors(G, i)
children(G::Graph{Directed}, i::Int) = outneighbors(G, i)

# TODO: These helpers are either already in polymake or should be optimized
# and put somewhere deeper in the Graph code.
# TODO: There is some amount of code duplication because all of these
# algorithms are variants of depth-first search.

function topological_sort(G::Graph{Directed})::Vector{Int}
  @req is_acyclic(G) "the digraph must be acyclic"
  degrees = indegree(G)
  todo = findall(i -> degrees[i] == 0, vertices(G))
  ordering = Int[]
  while !isempty(todo)
    v = popfirst!(todo)
    push!(ordering, v)
    append!(todo, [c for c in children(G, v) if all(w -> w in ordering, parents(G, c))])
  end
  return ordering
end

function descendants(G::Graph{Directed}, i::Int)::Set{Int}
  V = vertices(G)
  d = Set{Int}(i)
  todo = [ [setdiff(V, i), i] ]
  while !isempty(todo)
    cur = pop!(todo)
    W = cur[1]
    v = cur[2]
    C = children(G, v)
    union!(d, C)
    append!(todo, [ [setdiff(W, w), w] for w in C if w in W])
  end
  return d
end

nondescendants(G::Graph{Directed}, i::Int) = setdiff(vertices(G), descendants(G, i))

@doc raw"""
    pairwise_markov(G::Graph{Directed})

Return all *elementary* CI statements in the pairwise Markov property of `G`,
i.e., all `CIStmt([i],[j],K)` such that `j` is a non-parent as well as non-descendant
of `i` and `K = setdiff(nondescendants(G, i), j)`.

## Examples

```jldoctest
julia> G = graph_from_edges(Directed, [[1,2],[1,3],[2,3],[3,4]])
Directed graph with 4 nodes and the following edges:
(1, 2)(1, 3)(2, 3)(3, 4)

julia> pairwise_markov(G)
2-element Vector{CIStmt}:
 [1 _||_ 4 | {2, 3}]
 [2 _||_ 4 | {1, 3}]
```
"""
function pairwise_markov(G::Graph{Directed})::Vector{CIStmt}
  V = vertices(G)
  stmts = CIStmt[]
  for i in V
    nd = nondescendants(G, i)
    B = setdiff(nd, parents(G, i))
    append!(stmts, [ci_stmt(i, j, setdiff(nd, j)) for j in B])
  end
  unique(stmts)
end

@doc raw"""
    local_markov(G::Graph{Directed})

Return all *elementary* CI statements in the local Markov property of `G`,
i.e., all `CIStmt([i],[j],K)` such that `j` is a non-parent as well as non-descendent
of `i` and `K = parents(G, i)`.

## Examples

```jldoctest
julia> G = graph_from_edges(Directed, [[1,2],[1,3],[2,3],[3,4]])
Directed graph with 4 nodes and the following edges:
(1, 2)(1, 3)(2, 3)(3, 4)

julia> local_markov(G)
2-element Vector{CIStmt}:
 [1 _||_ 4 | 3]
 [2 _||_ 4 | 3]
```
"""
function local_markov(G::Graph{Directed})::Vector{CIStmt}
  V = vertices(G)
  stmts = CIStmt[]
  for i in V
    pa = parents(G, i)
    B = setdiff(nondescendants(G, i), pa)
    append!(stmts, [ci_stmt(i, j, pa) for j in B])
  end
  unique(stmts)
end

function is_ancestral(G::Graph{Directed}, A::Vector{Int})::Bool
  for v in setdiff(vertices(G), A)
    for w in A
      if has_edge(G, v, w) && !(v in A)
        return false
      end
    end
  end
  return true
end

function ancestral_closure(G::Graph{Directed}, A::Vector{Int})::Vector{Int}
  for v in setdiff(vertices(G), A)
    for w in A
      if has_edge(G, v, w) && !(v in A)
        # Add v and restart
        return ancestral_closure(G, [A..., v])
      end
    end
  end
  return A # if we get here, A was ancestral
end

function moralization(G::Graph{Directed}, A::Vector{Int})
  @req is_ancestral(G, A) "the moralization set must be ancestral"
  V = vertices(G)
  # Since Oscar's graphs can only have vertex sets of the form 1:m,
  # we have to very carefully relabel the edges so that vertex i in
  # the moralization corresponds to vertex A[i] in the original graph.
  M = Graph{Undirected}(length(A))
  for i in 1:length(A)
    for j in 1:length(A)
      if i == j
        continue
      end
      if has_edge(G, A[i], A[j]) || any(w -> has_edge(G, A[i], w) && has_edge(G, A[j], w), A)
        add_edge!(M, i, j)
      end
    end
  end
  return M
end

@doc raw"""
    are_separated(G::Graph{Directed}, i::Int, j::Int, K::Vector{Int})

Check if `i` and `j` are d-separated by `K` in the graph `G`, i.e., whether
`i` and `j` are separated by `K` in the moral graph `moralization(G, ancestral_closure(G, [i,j,K...]))`.

## Examples

```jldoctest
julia> G = graph_from_edges(Directed, [[1,2],[2,3],[2,4],[3,4]])
Directed graph with 4 nodes and the following edges:
(1, 2)(2, 3)(2, 4)(3, 4)

julia> are_separated(G, 1, 4, [2])
true

julia> are_separated(G, 1, 4, [3])
false
```
"""
function are_separated(G::Graph{Directed}, i::Int, j::Int, K::Vector{Int})::Bool
  @req i != j && !(i in K) && !(j in K) "i and j must be distinct and not in K"
  A = ancestral_closure(G, [i, j, K...])
  U = moralization(G, A)
  # The moralization has vertices indexed by 1:m which correspond to the
  # elements of A in order. We must test whether 1 is separated by 2 given
  # the elements 3:length(K)+2
  return are_separated(U, 1, 2, collect(3:(length(K)+2)))
end

@doc raw"""
    global_markov(G::Graph{Directed})

Return all *elementary* CI statements in the global Markov property of `G`,
i.e., all `CIStmt([i],[j],K)` such that `K` d-separates `i` and `j` in `G`
(see `are_separated`).

## Examples

```jldoctest
julia> G = graph_from_edges(Directed, [[1,2],[1,3],[2,3],[3,4]])
Directed graph with 4 nodes and the following edges:
(1, 2)(1, 3)(2, 3)(3, 4)

julia> global_markov(G)
4-element Vector{CIStmt}:
 [1 _||_ 4 | 3]
 [1 _||_ 4 | {2, 3}]
 [2 _||_ 4 | 3]
 [2 _||_ 4 | {1, 3}]
```

```jldoctest
julia> G = graph_from_edges(Directed, [[1,3],[2,3]])
Directed graph with 3 nodes and the following edges:
(1, 3)(2, 3)

julia> global_markov(G)
1-element Vector{CIStmt}:
 [1 _||_ 2 | {}]
```
"""
# TODO: Bayes Ball to get the CI statements more quickly?
# TODO: For now, this is precisely the same algorithm as for the
# undirected case because all the work is done in `are_separated`.
function global_markov(G::Graph{Directed})::Vector{CIStmt}
  [stmt for stmt in ci_statements(collect(vertices(G)))
   if are_separated(G, stmt.I[1], stmt.J[1], stmt.K)]
end
