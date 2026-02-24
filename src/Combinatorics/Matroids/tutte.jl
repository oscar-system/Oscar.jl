

function circuits(M::Matroid, S)
  return circuits(restriction(M,S))
end

function circuit(M::Matroid, S)
  #Check for independence
  @req rank(M) != length(matroid_groundset(M)) "The rank of the matroid should be smaller than the cardinality of the ground set"
  return circuits(M,S)[1] 
end

function cocircuits(M::Matroid, S)
  MD = dual_matroid(M)
  return circuits(MD,S)
end

function cocircuit(M::Matroid, S)
  #Check for independence
  @req !iszero(rank(M)) "The rank of the matroid should be smaller than the cardinality of the ground set"
  return cocircuits(M,S)[1] 
end




function tutte_group(R::Ring, M::Matroid)
  B = bases(M)
  idx = Dict{Set{Int}, Int}(Set(k) => i for (i,k) in enumerate(B))
  #idx[Set{Int}()] = length(B) + 1 #this is to index the epsilon
  gs = matroid_groundset(M)
  if characteristic(R)==2
            v = zeros(Int, length(B)+1)
            v[end] = 1 #this is for the epsilon
            relations = [v]
  else
            v = zeros(Int, length(B)+1)
            v[end] = 2
            relations = [v]
  end
  for X in nonbases(M)
    if rank(M,X) == rank(M)-1
      C = circuit(M,X)
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


export tutte_group, circuits, circuit, circuit

struct TutteGroup
  group::FinGenAbGroup
end

