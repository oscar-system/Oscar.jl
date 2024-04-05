#get all exponent vectors of a polynomial f
exponent_vectors = f->exponent_vector.(monomials(f), Ref(1))

#=
returns the set of bounding vectors of a marked Gröbner basis G, with markings given by Lm 

INPUT: A Gröbner basis G, a set of monomials Lm (the leading monomials of G)
OUTPUT: A set of n-dimensional integer vectors of the form a-b,
where a is the exponent vector of a leading monomoial of some g in G,
and b is the exponent vector of some other term
=# 
function difference_lead_tail(
  G::Oscar.IdealGens, Lm::Vector{L}
) where {L<:MPolyRingElem}
  lead_exp = Lm .|> exponent_vectors .|> first
  
  v = zip(lead_exp, exponent_vectors.(G)) .|> splat((l, t) -> Ref(l).-t)

  return unique!(reduce(vcat, v))[2:end] #temporary solution: the first element is always the zero vector? 
end


function difference_lead_tail(I::Oscar.IdealGens)
  lead_exp = leading_term.(I; ordering=ordering(I)) .|> exponent_vectors .|> first
  tail_exps = tail.(I; ordering=ordering(I)) .|> exponent_vectors
  
  v = zip(lead_exp, tail_exps) .|> splat((l, t) -> Ref(l).-t)

  return unique!(reduce(vcat, v))
end





#returns 'true' if Mv <_{lex} 0 , 'false' otherwise 
# <_{lex} is the lexicographic ordering on Q^n
#with the notation from my thesis, this is equivalent to v <_M 0
function new_less_than_zero(M::ZZMatrix, v::Vector{ZZRingElem})
  if is_zero(v)
      return false
  else
  i = 1
  end
  while dot(M[i, :], v) == 0
      i += 1
  end
  return dot(M[i, :], v) < 0
end

#returns all v in V with 0 <_S v  and v <_T 0 
#when V are the bounding vectors of a cone, this is the set of all "candidates" for the next facet normal  
function new_filter_by_ordering(start::MonomialOrdering, target::MonomialOrdering, V::Vector{Vector{ZZRingElem}})
    pred = v->(
      new_less_than_zero(canonical_matrix(target), ZZ.(v)) && !new_less_than_zero(canonical_matrix(start), ZZ.(v))
    )
    return unique!(filter(pred, V))
  end
  

#returns true if v <_M w , false otherwise  
#i.e elements of the vectors Mv and Mw are compared until a tie is broken

function matrix_less_than(M::ZZMatrix, v::Vector{ZZRingElem}, w::Vector{ZZRingElem})
    i = 1
    while dot(M[i, :], v) == dot(M[i, :], w) && i != number_of_rows(M) 
        i += 1 
    end
    return dot(M[i, :], v) < dot(M[i, :], w)
end

#returns true if u < v w.r.t the Facet preorder (cf. "The generic Gröbner walk" (Fukuda et a;. 2007), pg. 10)
#Comment: It may make more sense to have the monomial orderings as inputs.
function facet_less_than(S::ZZMatrix, T::ZZMatrix, u::Vector{ZZRingElem}, v::Vector{ZZRingElem})
    i = 1
    while dot(T[i,:], u)v == dot(T[i,:], v)u && i != number_of_rows(T)
        i += 1
    end
    return matrix_less_than(S, dot(T[i,:], u)v, dot(T[i,:], v)u) 
end

#returns all elements of V smaller than w w.r.t the facet preorder
function new_filter_lf(w::Vector{ZZRingElem}, start::MonomialOrdering, target::MonomialOrdering, V::Vector{Vector{ZZRingElem}})
    btz = Vector{Vector{ZZRingElem}}()
    for v in V
      if facet_less_than(canonical_matrix(start),canonical_matrix(target),w,v) &&  !(v in btz)
        push!(btz, v)
      end
    end
    return btz
end

#=
given a Gröbner basis G with markings given by Lm and the previous normal vector w, find the next one 
This method is compute_last_w from Fukuda 2007, pg. 12
Initialization is with w = ZZ.([0])
termination condition is V = []
=#  

function new_next_gamma(
    G::Oscar.IdealGens, Lm::Vector{L}, w::Vector{ZZRingElem}, start::MonomialOrdering, target::MonomialOrdering
  ) where {L<:MPolyRingElem}
    V = new_filter_by_ordering(start, target, [ZZ.(v) for v  in difference_lead_tail(G, Lm)])
    if w != ZZ.([0]) 
        V = new_filter_lf(w, start, target, V)
    end
    if isempty(V)
        return V
    end
    minV = first(V)
    for v in V 
        if facet_less_than(canonical_matrix(start), canonical_matrix(target),v, minV)
            minV = v
        end
    end
    return minV 
end





#=
function difference_lead_tail(
  G::Oscar.IdealGens, Lm::Vector{L}, T::MonomialOrdering
) where {L<:MPolyRingElem}
  lead_exp = leading_term.(Lm, ordering=T) .|> exponent_vectors .|> first
  
  v = zip(lead_exp, exponent_vectors.(G)) .|> splat((l, t) -> Ref(l).-t)

  return unique!(reduce(vcat, v))[2:end] #TEMPORARY SOLUTION 
end
=#