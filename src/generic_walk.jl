function generic_walk(G::Oscar.IdealGens, start::MonomialOrdering, target::MonomialOrdering)
  @vprintln :groebner_walk "Results for generic_walk"
  @vprintln :groebner_walk "Facets crossed for: "

  Lm = leading_term.(G, ordering = start)
  MG = markedGB(gens(G), Lm)
  v = markedGB_next_gamma(MG, ZZ.([0]), start, target)

  while !isempty(v)
    MG = markedGB_generic_step(MG, v, target)
    v = markedGB_next_gamma(MG, v, start, target)

    # TODO: increase step_counter here
    @vprintln :groebner_walk v
    @vprintln :groebner_walk 2 G
  end
  return Oscar.IdealGens(MG.gens, target; isGB = true)
end


function markedGB_generic_walk(G::Oscar.IdealGens, start::MonomialOrdering, target::MonomialOrdering)
  Lm = leading_term.(G, ordering = start)
  MG = markedGB(gens(G), Lm)
  v = markedGB_next_gamma(MG, ZZ.([0]), start, target)

  while !isempty(v)
    MG = markedGB_generic_step(MG, v, target)
    v = markedGB_next_gamma(MG, v, start, target)
  end
  return Oscar.IdealGens(MG.gens, target; isGB = true)
end

#Given the "old markedGB" GB and the newly computed facet normal v 
#compute the next markedGB by taking G.B of initial forms H w.r.t less 
#and lifting it with markedGB_lift_generic. Subsequently reduce 

function markedGB_generic_step(MG::markedGB, v::Vector{ZZRingElem}, ord::MonomialOrdering)
facet_Generators = markedGB_facet_initials(MG, v)
H = groebner_basis(
  ideal(facet_Generators); ordering=ord, complete_reduction=true, algorithm=:buchberger
)
newGB = markedGB_lift_generic(MG, gens(H), ord)
newGB = reductionalg(newGB)
return newGB
end


#= 

----- thesis example (chap3) 
R, (x,y,z) = polynomial_ring(QQ, ["x","y","z"])

o_s = lex(R)

o_t= weight_ordering([1,3,0], deglex(R))
S = canonical_matrix(o_s)
T = canonical_matrix(o_t)

I = ideal([x^2 + y*z, x*y + z^2])
G = groebner_basis(I, ordering = o_s, complete_reduction = true)

markedGB_generic_walk(G, o_s, o_t)


=#


#Auxiliary functions for the generic Gröbner walk, ordered by subroutine
#----------------------------------


#------next_gamma (Goal: Get next facet normal along generic path)

exponent_vectors = f->exponent_vector.(monomials(f), Ref(1)) #returns exponent vectors of all terms in f 



#returns a list of integer vectors of the form a - b
# (where a is a leading exponent and b is in the tail of some g in MG) 
function markedGB_difference_lead_tail(MG::markedGB)
    (G,Lm) = MG.gens, MG.markings 
  lead_exp = Lm .|> exponent_vectors .|> first
  
  v = zip(lead_exp, exponent_vectors.(G)) .|> splat((l, t) -> Ref(l).-t)

  return [ZZ.(v) for v in unique!(reduce(vcat, v))[2:end]] #temporary solution: the first element is always the zero vector? 
end

#given a marked GB, the previous weight vector w and monomial orderings
#returns the "next" facet normal
#i.e. the bounding vector v fulfilling w<v w.r.t facet preorder which is minimal (w.r.t facet preorder)

function markedGB_next_gamma(
    MG::markedGB, w::Vector{ZZRingElem}, start::MonomialOrdering, target::MonomialOrdering
  )
    V = new_filter_by_ordering(start, target, markedGB_difference_lead_tail(MG))
    if w != ZZ.([0]) 
        V = new_filter_lf(w, start, target, V)
    end
    if isempty(V)
        return V
    end
    minV = first(V)
    for v in V 
        if new_facet_less_than(canonical_matrix(start), canonical_matrix(target),v, minV)
            minV = v
        end
    end
    return minV 
end

#----------------------------

#------- facet initials 

#Given a markedGB MG and a facet normal v 
#Return corresponding G.B of initial forms by truncating to all bounding vectors parallel to v
function markedGB_facet_initials(MG::markedGB, v::Vector{ZZRingElem})
    inwG = copy(MG.markings)
    gens = copy(MG.gens)
    R = parent(first(gens))
    for i in 1:length(MG.markings)
      a = first(exponent_vector.(monomials(MG.markings[i]), Ref(1))) 
      for (b, coeff) in zip(exponent_vectors(gens[i]),
         [leading_coefficient(term) for term in terms(gens[i])])
        if new_is_parallel(ZZ.(a-b), v)
          inwG[i] += coeff*monomial(R,b)
        end
      end
    end
    return inwG
  end


#----------------------------

#------ lifting step 

#Given a markedGB MG, a reduced GB of initial forms H w.r.t ord, and a monomial order 
#"lift" H to a markedGB of I (ordering unknown!) by subtracting initial forms according to Fukuda, 2007
function markedGB_lift_generic(MG::markedGB, H::Vector{<:MPolyRingElem}, ord::MonomialOrdering)
    R = parent(first(MG.gens))
    Newlm = Array{elem_type(R),1}(undef, 0)
    liftPolys = Array{elem_type(R),1}(undef, 0)
    for g in H
      push!(Newlm, leading_term(g; ordering=ord))
      push!(liftPolys, g - normal_form(g, MG))
    end
    return markedGB(liftPolys, Newlm)
  end
  


#--------------------------




#helper functions 


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
  

#given two monomial orderings start and target and a vector of integer vectors V
#return all elements v of V with 0 <_T v and v <_S 0  
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
function new_facet_less_than(S::ZZMatrix, T::ZZMatrix, u::Vector{ZZRingElem}, v::Vector{ZZRingElem})
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
      if new_facet_less_than(canonical_matrix(start),canonical_matrix(target),w,v) &&  !(v in btz)
        push!(btz, v)
      end
    end
    return btz
end


#returns true if u is a non-zero integer multiple of v

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

#TODO: add tests 