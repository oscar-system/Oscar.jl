using Oscar

#include("new_next_gamma.jl")



#TODO: This function is asymmetric; For v = 0 and any u, new_is_parallel(u,v) = true but new_is_parallel(v,u) is false
#This is not a problem, as both u and v are non-zero in our setting
#Nevertheless, this should be fixed
function new_is_parallel(u::Vector{ZZRingElem}, v::Vector{ZZRingElem})
    count = 1
    x = 0
    for i = 1:length(u)
        if u[i] == 0
            if v[count] == 0
                count += +1
            else
                return false
            end
        else
            x = v[count] // u[i]
            count += 1
            break
        end
    end
    if count > length(v)
        return true
    end
    for i = count:length(v)
        @inbounds if v[i] != x * u[i]
            return false
        end
    end
    return true

end


#INPUT: A Gröbner basis G, the set lm of leading monomials of elements of G, a vector v 
#OUTPUT: The set of initial forms inw(G), where w is a vector lying on the facet defined by v
function new_facet_initials(G::Oscar.IdealGens, lm::Vector{T}, v::Vector{ZZRingElem}
  ) where {T<:MPolyRingElem}
inwG = copy(lm)
generators = gens(G)
  for i in 1:length(lm)
    a = first(exponent_vector.(monomials(lm[i]), Ref(1))) 
    for (b, coeff) in zip(exponent_vectors(generators[i]),
       [leading_coefficient(term) for term in terms(generators[i])])
      if new_is_parallel(ZZ.(a-b), v)
        inwG[i] += coeff*monomial(R,b)
      end
    end
  end
  return inwG
end


#= The lifting step, as described in "The generic Gröbner walk" (Fukuda et al 2007), pg.8 

    Given a set of initial forms inwG of the marked Gröbner basis G, 
    convert inwG to a Gröbner basis M of InwI w.r.t "target", and
    "lift" M to a Gröbner basis of I by subtracting from each m in M 
    its normal form w.r.t the starting basis G.

=# 

function new_lift_generic(G::Vector{T}, Lm::Vector{T}, inwG::Vector{T}, target::MonomialOrdering
) where {T<:MPolyRingElem}
  M = gens(groebner_basis(ideal(inwG), ordering = target, complete_reduction = true, algorithm=:buchberger))
  leading_newGB = Vector{MPolyRingElem}()
  newGB = Vector{MPolyRingElem}()
  for m in M
    push!(newGB, m - reduce_walk(m, G, Lm, target))
    push!(leading_newGB, leading_term(m, ordering = target))
  end
  return newGB, leading_newGB
end



#TODO: I still need to have a proper look at these functions. Improvements may be possible 

#INPUT: A polynomial p, the generators of a marked G.B G, the leading monomials Lm, a monomial order ord
#OUTPUT: The normal form of p w.r.t the marked Gröbner basis G 
#QUESTION: Where do I need ord? I think things all work without it
function reduce_walk(
  p::MPolyRingElem, G::Vector{T}, Lm::Vector{T}, ord::MonomialOrdering
) where {T<:MPolyRingElem}
  for i in 1:length(G)
    (q, b) = divides_walk(p, Lm[i], parent(p))
    if b
      return reduce_walk(p - (q * G[i]), G, Lm, ord)
    end
  end
  return p
end


function divides_walk(p::MPolyRingElem, lm::MPolyRingElem, S::MPolyRing)
  div = false
  newpoly = MPolyBuildCtx(S)
  for term in terms(p)
    (b, c) = divides(term, lm)
    if b
      push_term!(
        newpoly, first(coefficients(c)), first(AbstractAlgebra.exponent_vectors(c))
      )
      div = true
    end
  end
  return finish(newpoly), div
end



function new_generic_step(
  G::Oscar.IdealGens, Lm::Vector{T}, v::Vector{ZZRingElem}, ord::MonomialOrdering
) where {T<:MPolyRingElem}
  inwG = new_facet_initials(G, Lm, v)
 
  G, Lm = new_lift_generic(gens(G), Lm, inwG, ord)
  G = interreduce(G, Lm, ord)
  return G, Lm
end

function interreduce(
  G::Vector{T}, Lm::Vector{T}, ord::MonomialOrdering
) where {T<:MPolyRingElem}
  for i in 1:length(G)
    G[i] = reduce_walk(G[i], G[1:end .!= i], Lm[1:end .!= i], ord)
  end
  return G
end



