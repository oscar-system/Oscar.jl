using Oscar

R,(x,y) = polynomial_ring(QQ, [:x,:y])
I = ideal([y^4+ x^3-x^2+x,x^4])

set_verbosity_level(:groebner_walk, 1)

groebner_walk(I)
groebner_walk(I; algorithm=:generic)

