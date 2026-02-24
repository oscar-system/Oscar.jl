"""
  circuits(M::Matroid, S::Union{AbstractVector{T}, AbstractSet{T}}) where T<:GroundsetType

  Return the list of circuits of the matroid contained in `S`. 
"""
function circuits(M::Matroid, S::Union{AbstractVector{T}, AbstractSet{T}}) where T<:GroundsetType
  @req all(s in matroid_groundset(M) for s in S) "The restriction set has to be a subset of the matroid's ground set"
  return circuits(restriction(M,S))
end

function _circuit(M::Matroid, S::Union{AbstractVector{T}, AbstractSet{T}}) where T<:GroundsetType
  C = circuits(M, S)
    @req len(C) > 0 "S does not contain a circuit of the matroid"
  return C[1] 
end

"""
  cocircuits(M::Matroid, S::Union{AbstractVector{T}, AbstractSet{T}}) where T<:GroundsetType

  Return the list of cocircuits of matroid contained in S.
"""
function cocircuits(M::Matroid, S::Union{AbstractVector{T}, AbstractSet{T}}) where T<:GroundsetType
  MD = dual_matroid(M)
  return circuits(MD,S)
end

function _cocircuit(M::Matroid, S::Union{AbstractVector{T}, AbstractSet{T}}) where T<:GroundsetType
  
  @req !iszero(rank(M)) "The rank of the matroid should be smaller than the cardinality of the ground set"
  return cocircuits(M,S)[1] 
end

"""
  Computes the Tutte group of a matroid M for a ring with characteristic `char`.
  For more details; see [DW89](@cite).
"""
function tutte_group(R::Ring, M::Matroid)
  B = bases(M)
  idx = Dict{Set{Int}, Int}(Set(k) => i for (i,k) in enumerate(B))
  #idx[Set{Int}()] = length(B) + 1 #this is to index the epsilon
  gs = matroid_groundset(M)
  v = zeros(Int, length(B)+1)
  if char == 2
    v[end] = 1 #this is for the epsilon
  else
    v[end] = 2
  end
  relations = [v]
  for X in nonbases(M)
    if rank(M,X) == rank(M)-1
      C = _circuit(M,X)
      D = cocircuit(M,setdiff(gs,X))
      e = popfirst!(C)
      f = popfirst!(D)
      #push!(ret,D)
      for g in C
        for h in D
          v = zeros(Int, length(B)+1)
          I = setdiff(Set(X), [e,g])
          v[idx[union(I, [e,f])]] = 1
          v[idx[union(I, [e,h])]] = 1
          v[idx[union(I, [g,f])]] = -1
          v[idx[union(I, [g,h])]] = -1
          v[end] = sum([e<g, g<f, f<h, h<e]) #this is the index for the epsilon
          push!(relations, v)
       end
     end
    end
  end
  relations_matrix = matrix(ZZ, relations)
  return abelian_group(relations_matrix)
end

"""
  False if matroid M is not realizable over R.
  True if inconclusive. 
  Note that the Tutte group only yields a necessary (and no sufficient) criterion for realizability of M over R;
  see Corollary 1 in Section 3 of [DW89](@cite).
```julia
julia> tutte_realizable(QQ,uniform_matroid(2,4));
false
julia> tutte_realizable(QQ,fano_matroid())
true
```
"""
function tutte_realizable(R::Ring, M::Matroid)
  T = tutte_group(R,M);
  return tutte_realizable(T)
end

function tutte_realizable(G::FinGenAbGroup)
  return !is_one(G[end])
end

export tutte_group, circuits, cocircuit, tutte_realizable
