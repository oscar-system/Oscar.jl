using Oscar 
using Oscar.GroebnerWalk


#The generic walk currently doesn't work (type conflicts)
#This example tests "new_next_gamma"
#We perform computations on the same ideal throughout

R, (x,y,z) = polynomial_ring(QQ, ["x","y","z"])

set_verbosity_level(:groebner_walk, 1)
    

I = ideal([x^2 + y*z, x*y + z^2])

o_s = lex(R)

o_t = weight_ordering([1,3,0], deglex(R))

G = groebner_basis(I, ordering = o_s) 
S = canonical_matrix(o_s)
T = canonical_matrix(o_t)

new_next_gamma(G, ZZ.([0]), o_s, o_t)

Gfinal = groebner_walk(I, lex(R), o_t; walk_type = :generic) #throws an error




#= 

R2, (x,y) = polynomial_ring(QQ, ["x","y"])

I2 = ideal([x^2-y^3, x^3 -y^2 - x])

G2 = groebner_basis(I2)
#= Fukuda Jensen example 

