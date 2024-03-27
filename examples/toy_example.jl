using Oscar
using Oscar.GroebnerWalk

R,(x,y) = polynomial_ring(QQ, ["x","y"])
I = ideal([y^4+ x^3-x^2+x,x^4])

set_verbosity_level(:groebner_walk, 1)

groebner_walk(I)

S = degrevlex(R)
T = lex(R)

s = canonical_matrix(S)[1,:]
t = canonical_matrix(T)[1,:]

G = groebner_basis(I; ordering=degrevlex(R), complete_reduction=true)

G1 = standard_step(G, s, T)

next_weight(G1, ZZ.([1,1]), ZZ.([1,0]))

groebner_walk(I; walk_type=:generic)