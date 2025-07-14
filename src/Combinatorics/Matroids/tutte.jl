
export tutte_group, circuits, circuit, circuit

struct TutteGroup
  group::FinGenAbGroup
end

#=
    def __init__(self, M, char2 = False):
        
        B = list(M.bases())
        idx = {B[i] : i for i in range(len(B))}
        idx[frozenset()] = len(B)
        self._M = M
        self._bases = B
        self._idx = idx
        
        if char2:
            R = [{frozenset():1}]
        else:
            R = [{frozenset():2}]
        
        
        for X in M.nonbases():
            if M.rank(X)==M.full_rank()-1: 
                C=set(M.circuit(X))
                D=set(M.cocircuit(M.groundset()-X))
                e = C.pop()
                f = D.pop()
                for g in C:
                    for h in D:
                        R.append(self.cross_ratio(X-set([e,g]),e,g,f,h))
        self._TT = matrix(ZZ, 1,len(B)+1, sparse = True)
        
        self.add_relations(R)
        

=#
#=
M0 = matroid_from_nonbases([[3,4]],4)
tutte_group(M0)
M = pappus_matroid()
tutte_group(M)
circuit(M0,[3,4])
M1 = matroid_from_circuits([[1,2,3,4],[5,6,7,2],[1,5,6,7,3,4]],7)
tutte_group(M1)

NB = nonbases(M1)
circuits.(Ref(M1),NB)
circuits(M1,matroid_groundset(M1))
circuits(M1)
=#

function circuits(M::Matroid, S)
  return circuits(restriction(M,S))
end

function circuit(M::Matroid, S)
  #Check for independence
  return circuits(M,S)[1] 
end

function cocircuits(M::Matroid, S)
  cocircuits(M,
end



function tutte_group(M::Matroid)
  B = bases(M)
  gs = matroid_groundset(M)
  ret = []
  for X in nonbases(M)
    if rank(M,X) == rank(M)-1
      C = circuit(M,X)
      D = setdiff(gs,X)
      push!(ret,D)
    else
      println("hello there")
    end
  end
  return ret
end

