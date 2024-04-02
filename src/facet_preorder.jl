#returns 'true' if Mv <_{lex} 0 , 'false' otherwise 
# <_{lex} is the lexicographic ordering on Q^n
#with the notation from my thesis, this is equivalent to v <_M 0
function new_less_than_zero(M::ZZMatrix, v::Vector{ZZRingElem})
    i = 1
    while dot(M[i, :], v) == 0
        i += 1
    end
    return dot(M[i, :], v) < 0
end

#returns all v in V with 0 <_S v  and v <_T 0 
#when V are the bounding vectors of a cone, this is the set of all "candidates" for the next facet normal  
function new_filter_by_ordering(S::MonomialOrdering, T::MonomialOrdering, V::Vector{Vector{ZZRingElem}})
    pred = v->(
      new_less_than_zero(canonical_matrix(T), ZZ.(v)) && !new_less_than_zero(canonical_matrix(S), ZZ.(v))
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

#returns true if v < w w.r.t the Facet preorder (cf. "The generic Gröbner walk" (Fukuda et a;. 2007), pg. 10)
function facet_less_than(S::ZZMatrix, T::ZZMatrix, u::Vector{ZZRingElem}, v::Vector{ZZRingElem})
    i = 1
    while dot(T[i,:], u)v == dot(T[i,:], v)u && i != number_of_rows(T)
        i += 1
    end
    return matrix_less_than(S, dot(T[i,:], u)v, dot(T[i,:], v)u) 
end

#returns all elements of V smaller than w w.r.t the facet preorder
#QUESTION: Is there any particular reason to use sets here? 
function new_filter_lf(w::Vector{ZZRingElem}, S::ZZMatrix, T::ZZMatrix, V::Set{Vector{ZZRingElem}})
    btz = Set{Vector{ZZRingELem}}()
    for v in V
      if facet_less_than(S,T,v,w)
        push!(btz, v)
      end
    end
    return btz
end


#given a Gröbner basis G and the previous normal vector w, find the next one 
#This method is compute_last_w from Fukuda 2007, pg. 12
#Initialization is with w = ZZ.([0])
#termination condition is V = [] 
function new_next_gamma(
    G::Oscar.IdealGens, w::Vector{ZZRingElem}, start::MonomialOrdering, target::MonomialOrdering
  ) where {L<:MPolyRingElem}
    V = new_filter_by_ordering(start, target, [ZZ.(v) for v  in difference_lead_tail(G)])
    if w != ZZ.([0]) 
        V = new_filter_lf(w, S, T, V)
    end
    if isempty(V)
        return v
    end
    minV = first(V)
    for v in V 
        if facet_less_than(canonical_matrix(start), canonical_matrix(target),v, minV)
            minV = v
        end
    end
    return minV 
end

#= tests


R, (x,y,z) = polynomial_ring(QQ, ["x","y","z"])

o_s = lex(R)

o_t= weight_ordering([1,3,0], deglex(R))

S = canonical_matrix(o_s)
T = canonical_matrix(o_t)

V = [ZZ.([0, 3, -3]), ZZ.([1, -2, 1]), ZZ.([1,1,-2]), ZZ.([2,-1,-1])]

u = ZZ.([1,-2,1])
v = ZZ.([2,-1,-1])
filter_by_ordering(o_s, o_t, V)



matrix_less_than(S, ZZ.([3,2,1]), ZZ.([3,3,3])) 


I = ideal([x^2 + y*z, x*y + z^2])

Gfinal = groebner_walk(I)

=#