IP2 = projective_space(QQ, 2)
X = covered_scheme(IP2)
S = ambient_ring(IP2)
(u,v,w) = gens(S)
Ihom = ideal(S, u^2 - v*w)
I = IdealSheaf(IP2, Ihom)
U = patches(default_covering(X))
@test I(U[1]) isa Ideal 
@test I(U[2]) isa Ideal 
@test I(U[3]) isa Ideal 
V = PrincipalOpenSubset(U[1], gens(OO(U[1]))[1])
rho = I(U[1], V)
@test I(V) == ideal(OO(V), rho.(gens(I(U[1]))))
