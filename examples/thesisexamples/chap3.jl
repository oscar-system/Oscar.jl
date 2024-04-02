#Testing welpj's Groebner walk with my examples 

#Examples from chapter 3

using Oscar 
using Oscar.GroebnerWalk


#We perform computations on the same ideal throughout

R, (x,y,z) = polynomial_ring(QQ, ["x","y","z"])

set_verbosity_level(:groebner_walk, 1)


I = ideal([x^2 + y*z, x*y + z^2])

#First example: conversion from degrevlex to lex 

Gfinal = groebner_walk(I)

#The algorithm terminates after 1 standard step.
#Check the subroutines: 
G = groebner_basis(I; ordering=degrevlex(R), complete_reduction=true)

w = next_weight(G, ZZ.([1,1,1]), ZZ.([1,0,0]))

# I agree with the output of next_weight
# However, I think infoLevel should say "Crossed Cones in: [1,1,1]" (and not [1,0,0])

#First iteration 
inwI = ideal(initial_forms(G,w))
H = groebner_basis(inwI; ordering = weight_ordering(w, lex(R)), complete_reduction = true ) 
Gnew = lift2(G, degrevlex(R), H, weight_ordering(w, lex(R)))

#w != tau, so we go again, updating to Gnew, etc. 

w2 = next_weight(Gnew, ZZ.([1,1,1]), ZZ.([1,0,0]))

inwI2 = ideal(initial_forms(Gnew, w2))

H2 = groebner_basis(inwI2; ordering = weight_ordering(w2, lex(R)), complete_reduction = true ) 
G2new = lift2(Gnew, ordering(Gnew), H2, weight_ordering(w2, lex(R)))

#G2new is Glex, as w = tau 




#Second example (same conversion with different choice of weight vectors) 

o2 = weight_ordering([1,1,0], degrevlex(R))

o3 = weight_ordering([3,1,0], lex(R))

Gfinal2 = groebner_walk(I, o2, o3)

#gens(Gfinal) == gens(Gfinal2) 

#the conversion works but like before I don't like that it says "Crossed cones in [3,1,0]"
#[3,1,0] lies in the interior of a cone.


#Third example: convert from lex to deglex, refined by [1,3,0] 
o4 = weight_ordering([1,3,0], deglex(R))

Gfinal3 = groebner_walk(I, lex(R), o4)

#two cones are crossed. 
#Integer weight vectors is the standard choice. 


#leading_term.(gens(Gfinal3), ordering = ordering(Gfinal3))



#= 
groebner_walk(I, degrevlex(R), lex(R))

M = [1 3 0; 1 1 1; 1 0 0; 0 1 0; 0 0 1]

o2 = matrix_ordering([x,y,z], M)

result2 = groebner_walk(I, lex(R), o2)

gens(result2)

o3 = weight_ordering([1,3,0], deglex(R))

canonical_matrix(o2)

canonical_matrix(o3)

#leading_term.(gens(Gnew))
#tail.(gens(Gnew)))

=#