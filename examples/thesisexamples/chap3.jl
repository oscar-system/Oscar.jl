#Testing welpj's Groebner walk with my examples 

#Examples from chapter 3

using Oscar 


#We perform computations on the same ideal throughout

R, (x,y,z) = polynomial_ring(QQ, ["x","y","z"])

set_verbosity_level(:groebner_walk, 1)


I = ideal([x^2 + y*z, x*y + z^2])
os = lex(R)
ot = weight_ordering([1,3,0], deglex(R))
Glex = groebner_basis(I; ordering = ot, complete_reduction = true)

G1 = groebner_walk(I, ot, os, algorithm = :generic)

G2 = groebner_walk(I, ot, os, algorithm = :standard)

