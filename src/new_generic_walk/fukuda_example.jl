#The example from section 5 of "The generic Gröbner Walk" (Fukuda et al. 2007)

using Oscar


include("new_generic_walk.jl")

R, (x,y) = polynomial_ring(QQ, ["x","y"])

I = ideal([x^2 - y^3, x^3 - y^2 - x])

start = degrevlex(R)
target = lex(R)
G = groebner_basis(I, ordering = start)
new_generic_walk(G, start, target) #success! 



#= Go through the example step by step
Lm = leading_term.(G, ordering = start)

w = new_next_gamma(G, Lm, [ZZ.(0)], start, target)

inwG = new_facet_initials(G, Lm, w)

(G2, Lm2) = new_generic_step(G, Lm, w, target)

G2 = Oscar.IdealGens(G2)

w2 = new_next_gamma(G2, Lm2, w , start, target)

inw2G = new_facet_initials(G2, Lm2, w2)

new_lift_generic(gens(G2), Lm2, inw2G, target) #success! This is correct 

(G, Lm) = new_generic_step(G, Lm, w2, target)

G = Oscar.IdealGens(G)

w3 = new_next_gamma(G, Lm, w2, start, target)

inwG = new_facet_initials(G, Lm, w3)

(G,Lm) = new_generic_step(G,Lm,w3,target)

G = Oscar.IdealGens(G)

w4 = new_next_gamma(G,Lm,w3, start, target)

terminate
=#