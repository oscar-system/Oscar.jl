#=
    Implicitization of Bezier surface. Example taken from Amrhein, Gloor, Küchlin. "On the walk" (2004)
=#
using Oscar

R, (x,y,z,u,v) = polynomial_ring(QQ, [:x,:y,:z,:u,:v])

o1 = matrix_ordering(R, [1 1 1 0 0; 0 0 0 1 1; 0 0 0 1 0; 1 1 0 0 0; 1 0 0 0 0])
o2 = matrix_ordering(R, [0 0 0 1 1; 1 1 1 0 0; 1 1 0 0 0; 1 0 0 0 0; 0 0 0 1 0])

I = ideal([
    u + u^2 - 2*v - 2*u^2*v + 2*u*v^2 - x,
    -6*u + 2*v + v^2 - 5*v^3 + 2*u*v^2 - 4*u^2*v^2 - y,
    -2 + 2*u^2 + 6*v - 3*u^2*v^2 - z
])
set_verbosity_level(:groebner_walk, 1)

G = groebner_walk(I, o2, o1; algorithm=:standard)

