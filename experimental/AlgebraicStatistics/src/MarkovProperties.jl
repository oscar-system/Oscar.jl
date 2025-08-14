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
  V = vertices(G)
  [ci_stmt([i], [j], setdiff(V, [i,j])) for i in 1:length(V) for j in (i+1):length(V)
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
  V = vertices(G)
  unique([ci_stmt([i], [j], neighbors(G, i)) for i in 1:length(V) for j in 1:length(V)
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
 [1 _||_ 3 | {2, 5}]
 [1 _||_ 3 | {2, 4}]
 [1 _||_ 3 | {2, 4, 5}]
 [1 _||_ 5 | 4]
 [1 _||_ 5 | {2, 4}]
 [1 _||_ 5 | {3, 4}]
 [1 _||_ 5 | {2, 3, 4}]
 [2 _||_ 4 | 1]
 [2 _||_ 4 | {1, 5}]
 [2 _||_ 4 | {1, 3}]
 [2 _||_ 4 | {1, 3, 5}]
 [2 _||_ 5 | 1]
 [2 _||_ 5 | 4]
 [2 _||_ 5 | {1, 4}]
 [2 _||_ 5 | {3, 4}]
 [2 _||_ 5 | {1, 3}]
 [2 _||_ 5 | {1, 3, 4}]
 [3 _||_ 4 | 1]
 [3 _||_ 4 | 2]
 [3 _||_ 4 | {1, 5}]
 [3 _||_ 4 | {2, 5}]
 [3 _||_ 4 | {1, 2}]
 [3 _||_ 4 | {1, 2, 5}]
 [3 _||_ 5 | 1]
 [3 _||_ 5 | 2]
 [3 _||_ 5 | 4]
 [3 _||_ 5 | {1, 4}]
 [3 _||_ 5 | {2, 4}]
 [3 _||_ 5 | {1, 2}]
 [3 _||_ 5 | {1, 2, 4}]
```
"""
function global_markov(G::Graph{Undirected})::Vector{CIStmt}
  V = vertices(G)
  stmts = CIStmt[]
  for i in 1:length(V)
    for j in 1:length(V)
      if i == j
        continue
      end
      M = setdiff(V, [i,j])
      for k in 0:length(M)
        for K in subsets(M, k)
          if are_separated(G, i, j, K)
            push!(stmts, ci_stmt([i], [j], K))
          end
        end
      end
    end
  end
  return unique(stmts)
end
