using Oscar
include("markedGB.jl")
include("markedGB_helpers.jl")






function markedGB_generic_walk(G::Oscar.IdealGens, start::MonomialOrdering, target::MonomialOrdering)
    Lm = leading_term.(G, ordering = start)
    MG = markedGB(gens(G), Lm)
    v = markedGB_next_gamma(MG, ZZ.([0]), start, target)
  
    while !isempty(v)
      MG = markedGB_generic_step(MG, v, target)
      v = markedGB_next_gamma(MG, v, start, target)
    end
    return Oscar.IdealGens(MG.gens, o_t; isGB = true)
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


#= 