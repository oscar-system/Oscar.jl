@doc raw"""
    circuits(M::Matroid, S:: GroundsetType)

Return the list of circuits of the matroid contained in `S`. 

# Example
```jldoctest
julia> M = uniform_matroid(3,7);

julia> S = [1,3,4,6,7];

julia> circuits(M,S)
5-element Vector{Vector{Int64}}:
 [1, 3, 4, 6]
 [1, 3, 4, 7]
 [1, 3, 6, 7]
 [1, 4, 6, 7]
 [3, 4, 6, 7]
```
"""
function circuits(M::Matroid, S::T) where T<:GroundsetType
  @req issubset(S, matroid_groundset(M)) "The restriction set has to be a subset of the matroid's ground set"
  return circuits(restriction(M,S))
end

@doc raw"""
    cocircuits(M::Matroid, S::T) where T<:GroundsetType

Return the list of cocircuits of matroid contained in `S`.

# Example
```jldoctest
julia> M = uniform_matroid(3,7);

julia> S = [1,3,4,6,7];

julia> cocircuits(M,S)
1-element Vector{Vector{Int64}}:
 [1, 3, 4, 6, 7]
```
"""
function cocircuits(M::Matroid, S::T) where T<:GroundsetType
  MD = dual_matroid(M)
  return circuits(MD,S)
end

@doc raw"""
    tutte_group(M::Matroid; char::Int=0)

Computes the Tutte group of a matroid `M` over a ring with characteristic `char`.
It should be noted, that the `char` only matters if it is two.
For more details see [DW89](@cite).

# Example
```jldoctest
julia> T = tutte_group(fano_matroid());

julia> ngens(T)
29
```
"""
function tutte_group(M::Matroid; char::Int=0)
  B = bases(Int, M)
  idx = Dict{BitSet, Int}(BitSet(b) => i for (i,b) in enumerate(B))
  v = zeros(Int, length(B)+1)
  if char == 2
    v[end] = 1 #this is for the epsilon
  else
    v[end] = 2 #this is for the epsilon
  end
  relations = [v]
  all_circuits = circuits(Int, M)
  all_cocircuits = cocircuits(Int, M)
  all_circ_bs = BitSet.(all_circuits)
  all_cocirc_bs = BitSet.(all_cocircuits)
  n = length(M)
  for X in nonbases(Int, M)
    Xbs = BitSet(X)
    ci = findfirst(c -> issubset(c, Xbs), all_circ_bs)
    ci === nothing && continue
    findnext(c -> issubset(c, Xbs), all_circ_bs, ci+1) !== nothing && continue
    Ybs = BitSet(1:n); setdiff!(Ybs, Xbs)
    di = findfirst(c -> issubset(c, Ybs), all_cocirc_bs)
    di === nothing && continue
    C = copy(all_circuits[ci])
    D = copy(all_cocircuits[di])
    e = popfirst!(C)
    f = popfirst!(D)
    for g in C
      for h in D
        v = zeros(Int, length(B)+1)
        I = setdiff(Xbs, (e, g))
        v[idx[union(I, (e, f))]] = 1
        v[idx[union(I, (g, h))]] = 1
        v[idx[union(I, (e, h))]] = -1
        v[idx[union(I, (g, f))]] = -1
        v[end] = count([e<f, f<g, g<h, h<e]) #this is the index for the epsilon
        push!(relations, v)
      end
    end
  end
  relations_matrix = matrix(ZZ, relations)
  return abelian_group(relations_matrix)
end

@doc raw"""
    is_tutte_realizable(M::Matroid)

Returns whether the matroid satisfies the Tutte realizability condition.
If `false`, this implies that the matroid cannot be realized over a field with characteristic other than `2`.
If `true`, we don't have a conclusive answer on realizability, since the Tutte
group only yields a necessary (and no sufficient) criterion for realizability
of `M`; see Corollary 1 in Section 3 of [DW89](@cite).

# Example
```jldoctest
julia> is_tutte_realizable(uniform_matroid(2,4))
true

julia> is_tutte_realizable(fano_matroid())
false
```
"""
function is_tutte_realizable(M::Matroid)
  T = tutte_group(M)
  return is_tutte_realizable(T)
end

function is_tutte_realizable(G::FinGenAbGroup)
  n=ngens(G)
  return !is_one(G[n])
end
