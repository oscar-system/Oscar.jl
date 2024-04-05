using Oscar

include("new_generic_walk.jl")

R, (x,y,z) = polynomial_ring(QQ, ["x","y","z"])

o_s = lex(R)

o_t= weight_ordering([1,3,0], deglex(R))
S = canonical_matrix(o_s)
T = canonical_matrix(o_t)

I = ideal([x^2 + y*z, x*y + z^2])
G = groebner_basis(I, ordering = o_s)

new_generic_walk(G, o_s, o_t)

#=
Tests
Lm = leading_term.(G, ordering = o_s)
w = new_next_gamma(G, Lm, ZZ.([0]), o_s, o_t)
(newG, newLm) =  new_generic_step(G, Lm, w, o_t)


newG = Oscar.IdealGens(newG)

difference_lead_tail(newG, Lm)

w2 = new_next_gamma(newG, newLm, w, o_s, o_t)

(newnewG, newnewLm) = new_generic_step(newG, newLm, w2, o_t)

newnewG = Oscar.IdealGens(newnewG)

w3 = new_next_gamma(newnewG, newnewLm, w2,o_s, o_t) #throws an error. It should give an empty list (or some other signal that we're done) 

(finalG, finalLm) = new_generic_step(newnewG, newnewLm, w2, o_t)
finalG = Oscar.IdealGens(finalG) 

finalG, newnewG
new_filter_lf(w2, S, T, Vector(Vector{ZZRingElem}[]))

generic_walk(G, o_s, o_t)

difference_lead_tail(newnewG, newnewLm, o_t) #confusing output 

=# 