@doc raw"""
    circuits(M::Matroid, S::T) where T<:GroundsetType

Return the list of circuits of the matroid contained in `S`. 
"""
function circuits(M::Matroid, S::T) where T<:GroundsetType
  @req issubset(S, matroid_groundset(M)) "The restriction set has to be a subset of the matroid's ground set"
  return circuits(restriction(M,S))
end

function _circuit(M::Matroid, S::T) where T<:GroundsetType
  C = circuits(M, S)
  @req !isempty(C) "S does not contain a circuit of the matroid"
  return C[1] 
end

@doc raw"""
    cocircuits(M::Matroid, S::T) where T<:GroundsetType

Return the list of cocircuits of matroid contained in `S`.
"""
function cocircuits(M::Matroid, S::T) where T<:GroundsetType
  MD = dual_matroid(M)
  return circuits(MD,S)
end

function _cocircuit(M::Matroid, S::T) where T<:GroundsetType
  @req !iszero(rank(M)) "The rank of the matroid should be smaller than the cardinality of the ground set"
  return cocircuits(M,S)[1] 
end

@doc raw"""
    tutte_group(M::Matroid; char::Int=0)

Computes the Tutte group of a matroid `M` over a ring with characteristic `char`.
It should be noted, that the `char` only matters if it is two.
For more details; see [DW89](@cite).
"""
function tutte_group(M::Matroid; char::Int=0)
  B = bases(M)
  idx = Dict{Set{Int}, Int}(Set(k) => i for (i,k) in enumerate(B))
  gs = matroid_groundset(M)
  v = zeros(Int, length(B)+1)
  if char == 2
    v[end] = 1 #this is for the epsilon
  else
    v[end] = 2 #this is for the epsilon
  end
  relations = [v]
  for X in nonbases(M)
    if rank(M,X) == rank(M)-1
      C = _circuit(M,X)
      D = _cocircuit(M,setdiff(gs,X))
      e = popfirst!(C)
      f = popfirst!(D)
      #push!(ret,D)
      for g in C
        for h in D
          v = zeros(Int, length(B)+1)
          I = setdiff(Set(X), [e,g])
          v[idx[union(I, [e,f])]] = 1
          v[idx[union(I, [g,h])]] = 1
          v[idx[union(I, [e,h])]] = -1
          v[idx[union(I, [g,f])]] = -1
          v[end] = count([e<g, g<f, f<h, h<e]) #this is the index for the epsilon
          push!(relations, v)
       end
     end
    end
  end
  relations_matrix = matrix(ZZ, relations)
  return abelian_group(relations_matrix) #this implicitely computes the Smith normal form
end

@doc raw"""
    is_tutte_realizable(M::Matroid)

Returns whether the matroid fulfills the Tutte realizability condition.
If `false`, this implies that the matroid cannot be realized over a field with characteristic other than `2`.
If `true`, we don't have a conclusive answer on realizability, since the Tutte
group only yields a necessary (and no sufficient) criterion for realizability
of M; see Corollary 1 in Section 3 of [DW89] (@cite).

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
